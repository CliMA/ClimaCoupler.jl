using ClimaAnalysis, Dates
import ClimaCalibrate
import ClimaCoupler
import JLD2
import EnsembleKalmanProcesses as EKP
import NaNStatistics
import CairoMakie
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
            G_ensemble[:, m] .= NaN
        end
    end
    return G_ensemble
end

function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    plot_constrained_params_and_errors(output_dir, ekp, priors)
end

function plot_constrained_params_and_errors(output_dir, ekp, priors)
    dim_size = sum(length.(EKP.batch(prior)))
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500))
    for i in 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size + 1], ekp)
    CairoMakie.save(joinpath(output_dir, "constrained_params_and_error.png"), fig)
    return nothing
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
    ts = process_outputvar(simdir, "ts")
    lwp = process_outputvar(simdir, "lwp")
    lwp = window.(lwp, "latitude"; left = -60, right = 60)

    # Map over each year
    year_observations = map(1:4:length(rsut)) do year_start
        year_end = min(year_start + 3, length(rsut))
        yr_ind = year_start:year_end

        net_rad_yr = mean(mean.(getproperty.(net_rad[yr_ind], :data)))

        rsut_yr = vectorize(rsut[yr_ind])
        rlut_yr = vectorize(rlut[yr_ind])
        cre_yr = vectorize(cre[yr_ind])
        pr_yr = vectorize(pr[yr_ind])
        ts_yr = vectorize(ts[yr_ind])
        lwp_yr = vectorize(lwp[yr_ind])

        vcat(net_rad_yr, cre_yr, rsut_yr, rlut_yr, pr_yr, ts_yr, lwp_yr)
    end
    return vcat(year_observations...)
end

function vectorize(seasonal_avgs)
    return vcat(vec.(getproperty.(seasonal_avgs, :data))...)
end

# Process an outputvar into a vector of seasonal averages
function process_outputvar(simdir, name)
    monthly_avgs = preprocess_monthly_averages(simdir, name)
    seasons = split_monthly_averages_into_seasons(monthly_avgs) # replace with split_by_season after testing
    # Ensure each season has three months
    @assert all(map(x -> length(times(x)) == 3, seasons))
    seasonal_avgs = average_time.(seasons)
    return seasonal_avgs
end

# Preprocess monthly averages to the right dimensions and dates, remove NaNs
function preprocess_monthly_averages(simdir, name)
    monthly_avgs = get_monthly_averages(simdir, name)
    # Interpolate to pressure coordinates to match observations
    pressure = get_monthly_averages(simdir, "pfull")
    if has_altitude(monthly_avgs)
        # This fails sometimes
        monthly_avgs = ClimaAnalysis.Atmos.to_pressure_coordinates(monthly_avgs, pressure)
        monthly_avgs = limit_pressure_dim_to_era5_range(monthly_avgs)
    end
    monthly_avgs.attributes["start_date"] = pressure.attributes["start_date"]
    # Line up dates for monthly averages
    monthly_avgs = ClimaAnalysis.shift_to_start_of_previous_month(monthly_avgs)
    # Remove spinup time
    monthly_avgs = window(monthly_avgs, "time"; left = spinup_time)
    global_mean = monthly_avgs |> average_lat |> average_lon |> average_time
    FT = monthly_avgs.data |> eltype
    # Replace NaNs with global mean
    monthly_avgs = ClimaAnalysis.replace(monthly_avgs, NaN => FT(NaNStatistics.nanmean(global_mean.data)))
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
