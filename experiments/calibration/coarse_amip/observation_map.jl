using ClimaAnalysis, Dates
import ClimaCalibrate
import ClimaCoupler
import JLD2
import EnsembleKalmanProcesses as EKP
import NaNStatistics
import CairoMakie
include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/observation_utils.jl"))
include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/leaderboard/leaderboard.jl"))

function ClimaCalibrate.observation_map(iteration)
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))
    current_minibatch = EKP.get_current_minibatch(ekp)
    single_obs_len = sum(length(EKP.get_obs(ekp))) 
    single_member_len = single_obs_len * length(current_minibatch)
    ensemble_size = EKP.get_N_ens(ekp)
    G_ensemble = Array{Float64}(undef, single_member_len, ensemble_size)
    for m in 1:ensemble_size
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, "model_config/output_active")
        @info "Processing member $m: $simdir_path"
        try
            G_ensemble[:, m] .= process_member_data(SimDir(simdir_path), current_minibatch)

        catch e
            @error "Error processing member $m, filling observation map entry with NaNs" exception = e
            G_ensemble[:, m] .= NaN
        end
    end
    return G_ensemble
end


function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    try
        plot_output_path = ClimaCalibrate.path_to_iteration(output_dir, iteration)
        plot_constrained_params_and_errors(plot_output_path, ekp, prior)
        
        for m in 1:EKP.get_N_ens(ekp)
            output_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
            diagnostics_folder_path = joinpath(output_path, "model_config", "output_active")
            compute_leaderboard(output_path, diagnostics_folder_path, spinup_months)
        end
    catch e
        @error "Error in `analyze_iteration`" exception = catch_backtrace()
    end
end

function plot_constrained_params_and_errors(output_dir, ekp, prior)
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
function process_member_data(simdir::SimDir, current_minibatch)
    rsdt = preprocess_diagnostic_monthly_averages(simdir, "rsdt")
    rsut = preprocess_diagnostic_monthly_averages(simdir, "rsut")
    rlut = preprocess_diagnostic_monthly_averages(simdir, "rlut")
    rsutcs = preprocess_diagnostic_monthly_averages(simdir, "rsutcs")
    rlutcs = preprocess_diagnostic_monthly_averages(simdir, "rlutcs")

    net_rad = rlut + rsut - rsdt  |> average_lat |> average_lon |> get_yearly_averages
    pr = preprocess_diagnostic_monthly_averages(simdir, "pr")
    ts = preprocess_diagnostic_monthly_averages(simdir, "ts")
    lwp = preprocess_diagnostic_monthly_averages(simdir, "lwp")
    lwp = window(lwp, "latitude"; left = -60, right = 60)

    year_observations = map(1:length(current_minibatch)) do yr
        rsut_yr = year_of_seasonal_averages(rsut, yr)
        rlut_yr = year_of_seasonal_averages(rlut, yr)
        # This needs to be averaged over lat, lon, and seasons
        cre_yr = year_of_seasonal_averages(rsut + rlut - rsutcs - rlutcs, yr)

        pr_yr = year_of_seasonal_averages(pr, yr)
        ts_yr = year_of_seasonal_averages(ts, yr)
        lwp_yr = year_of_seasonal_averages(lwp, yr)
        return vcat([vec(net_rad[yr]), cre_yr, rsut_yr, rlut_yr, pr_yr, ts_yr, lwp_yr]...)
    end
    return vcat(year_observations...)
end

# Preprocess monthly averages to the right dimensions and dates, remove NaNs
function preprocess_diagnostic_monthly_averages(simdir, name)
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
