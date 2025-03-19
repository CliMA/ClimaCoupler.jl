using ClimaAnalysis, Dates
import ClimaCalibrate
import ClimaCoupler
import JLD2
import EnsembleKalmanProcesses as EKP
include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/observation_utils.jl"))

function ClimaCalibrate.observation_map(iteration)
    observation_vec = JLD2.load_object(observation_path)
    single_member_dims = length(EKP.get_obs(first(observation_vec))) * batch_size
    G_ensemble = Array{Float64}(undef, single_member_dims, ensemble_size)
    for m in 1:ensemble_size
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, "model_config/output_active")
        @info "Processing member $m: $simdir_path"
        try
            G_ensemble[:, m] .= process_member_data(SimDir(simdir_path))

        catch e
            @error "Error processing member $m, filling observation map entry with NaNs" exception = e
            bt = catch_backtrace()
            println("Stacktrace:")
            display(stacktrace(bt))
            G_ensemble[:, m] .= NaN
        end
    end
    return G_ensemble
end

# Process a single ensemble member's data into a vector
function process_member_data(simdir::SimDir)
    rsut = process_outputvar(simdir, "rsut")
    rlut = process_outputvar(simdir, "rlut")
    rsutcs = process_outputvar(simdir, "rsutcs")
    rlutcs = process_outputvar(simdir, "rlutcs")
    rsdt = process_outputvar(simdir, "rsdt")
    # This needs to be averaged over lat, lon, and seasons
    net_rad = rlut + rsut - rsdt
    cre = rsut + rlut - rsutcs - rlutcs

    pr = process_outputvar(simdir, "pr")
    shf = process_outputvar(simdir, "shf")
    ts = process_outputvar(simdir, "ts")

    ta = process_outputvar(simdir, "ta")
    hur = process_outputvar(simdir, "hur")
    hus = process_outputvar(simdir, "hus")

    ql = process_outputvar(simdir, "ql")
    qi = process_outputvar(simdir, "qi")

    # Map over each year
    year_observations = map(1:4:length(rsut)) do year_start
        year_end = min(year_start + 3, length(rsut))
        yr_ind = year_start:year_end

        net_rad_yr = mean(mean.(getproperty.(net_rad[yr_ind], :data)))

        rsut_yr = downsample_and_vectorize(rsut[yr_ind])
        rlut_yr = downsample_and_vectorize(rlut[yr_ind])
        cre_yr = downsample_and_vectorize(cre[yr_ind])
        pr_yr = downsample_and_vectorize(pr[yr_ind])
        shf_yr = downsample_and_vectorize(shf[yr_ind])
        ts_yr = downsample_and_vectorize(ts[yr_ind])

        ta_yr = downsample_and_vectorize(ta[yr_ind])
        hur_yr = downsample_and_vectorize(hur[yr_ind])
        hus_yr = downsample_and_vectorize(hus[yr_ind])
        ql_yr = downsample_and_vectorize(ql[yr_ind])
        qi_yr = downsample_and_vectorize(qi[yr_ind])

        vcat(net_rad_yr, rsut_yr, rlut_yr, cre_yr, pr_yr, shf_yr, ts_yr, ta_yr, hur_yr, hus_yr, ql_yr, qi_yr)
    end
    return vcat(year_observations...)
end

function downsample_and_vectorize(seasonal_avgs)
    downsampled_seasonal_avg_arrays = downsample.(seasonal_avgs, 3)
    return vcat(vec.(downsampled_seasonal_avg_arrays)...)
end

# Process an outputvar into a vector of seasonal averages
function process_outputvar(simdir, name)
    monthly_avgs = preprocess_monthly_averages(simdir, name)
    seasons = split_monthly_averages_into_seasons(monthly_avgs)
    # Ensure each season has three months
    @assert all(map(x -> length(times(x)) == 3, seasons))
    seasonal_avgs = average_time.(seasons)
    return seasonal_avgs
end

# Preprocess monthly averages to the right dimensions and dates, remove NaNs
function preprocess_monthly_averages(simdir, name)
    monthly_avgs = get_monthly_averages(simdir, name)
    # Interpolate to pressure coordinates to match observations
    if has_altitude(monthly_avgs)
        pressure = get_monthly_averages(simdir, "pfull")
        monthly_avgs = ClimaAnalysis.Atmos.to_pressure_coordinates(monthly_avgs, pressure)
        monthly_avgs = limit_pressure_dim_to_era5_range(monthly_avgs)
    end
    # Line up dates for monthly averages
    monthly_avgs = ClimaAnalysis.shift_to_start_of_previous_month(monthly_avgs)
    # Remove spinup time
    monthly_avgs = window(monthly_avgs, "time"; left = spinup_time)
    global_mean = monthly_avgs |> average_lat |> average_lon |> average_time
    # Replace NaNs with global mean
    monthly_avgs = ClimaAnalysis.replace(monthly_avgs, NaN => global_mean)
    return monthly_avgs
end

function split_monthly_averages_into_seasons(monthly_avgs)
    all_times = times(monthly_avgs)
    @assert length(all_times) % 3 == 0

    # Window over 3 months at a time to create seasonal outputvars
    split_by_seasons = map(all_times[1:3:end]) do t
        window(monthly_avgs, "time"; left = t, right = t + 2months)
    end
    return split_by_seasons
end

function rebalance_months_per_season!(seasons)
    months_per_season = 3
    # Keep track of which seasons have been modified
    modified_indices = Set{Int}()

    for i in 1:(length(seasons) - 1)
        s1, s2 = seasons[i], seasons[i + 1]
        s1_times = s1 |> times
        s2_times = s2 |> times
        s1_months = length(s1_times)
        s2_months = length(s2_times)

        s1_months == s2_months == months_per_season && continue

        if s1_months == 4 && s2_months == 2
            # Pop last month of season 1, prepend to season 2
            # Copy the data so it's not a subarray
            first_three_months = copy(window(s1, "time"; left = s1_times[1], right = s1_times[3]).data)
            s1_new_dims = s1.dims
            s1_new_dims[time_name(s1)] = s1_times[1:3]
            s1_remake = remake(s1; dims = s1_new_dims, data = first_three_months)

            first_month_s2 = copy(window(s1, "time"; left = s1_times[end]).data)
            s2_new_dims = s2.dims
            s2_new_dims[time_name(s2)] = [s1_times[end], s2_times...]
            s2_new_data = vcat(first_month_s2, s2.data)
            s2_remake = remake(s2; dims = s2_new_dims, data = s2_new_data)

            seasons[i] = s1_remake
            seasons[i + 1] = s2_remake

            push!(modified_indices, i, i + 1)

        elseif s1_months == 2 && s2_months == 4
            # Pop first month of season 2, append to season 1
            first_month_s2 = copy(window(s2, "time"; left = s2_times[1], right = s2_times[1]).data)
            s1_new_dims = s1.dims
            s1_new_dims[time_name(s1)] = [s1_times..., s2_times[1]]
            s1_new_data = vcat(s1.data, first_month_s2)
            s1_remake = remake(s1; dims = s1_new_dims, data = s1_new_data)

            last_three_months = copy(window(s2, "time"; left = s2_times[2], right = s2_times[end]).data)

            s2_new_dims = s2.dims
            s2_new_dims[time_name(s2)] = s2_times[2:end]
            s2_remake = remake(s2; dims = s2_new_dims, data = last_three_months)

            seasons[i] = s1_remake
            seasons[i + 1] = s2_remake

            push!(modified_indices, i, i + 1)

        else
            @info "Only one month has an imbalance. The next pair ($(i+1), $(i+2)) should address this"
        end
    end

    # Log which seasons were modified
    if !isempty(modified_indices)
        @info "Modified seasons at indices: $(sort(collect(modified_indices)))"
    else
        @info "No seasons needed modification"
    end

end
