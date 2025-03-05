using ClimaAnalysis, Dates
import ClimaCalibrate
import ClimaCoupler
import JLD2
import EnsembleKalmanProcesses as EKP
import NaNStatistics
# import CairoMakie
include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/observation_utils.jl"))
include(joinpath(pkgdir(ClimaCoupler), "experiments/ClimaEarth/leaderboard/leaderboard.jl"))

function ClimaCalibrate.observation_map(iteration)
    ekp = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir, iteration))
    current_minibatch = EKP.get_current_minibatch(ekp)
    obs = EKP.get_obs(ekp)
    single_obs_len = sum(length(obs))
    single_member_len = single_obs_len * length(current_minibatch)
    ensemble_size = EKP.get_N_ens(ekp)
    short_names = split(obs_series.observations[1].names[1], ";") # This relies on the naming convention

    G_ensemble = Array{Float64}(undef, single_member_len, ensemble_size)
    for m in 1:ensemble_size
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, "model_config")
        @info "Processing member $m: $simdir_path"
        try
            G_ensemble[:, m] .= process_member_data(SimDir(simdir_path), short_names, current_minibatch)

        catch e
            @error "Error processing member $m, filling observation map entry with NaNs" exception = e
            G_ensemble[:, m] .= NaN
        end
    end
    total_elements = length(G_ensemble)
    nan_count = count(isnan, G_ensemble)
    # Check for 50% nans
    @assert nan_count < total_elements / 2
    # TODO: use nanmean
    @info "Mean bias y - G, averaged across the ensemble" bias = mean(G_ensemble, dims = 2) - obs |> mean
    return G_ensemble
end

function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    plot_output_path = ClimaCalibrate.path_to_iteration(output_dir, iteration)
    plot_constrained_params_and_errors(plot_output_path, ekp, prior)

    for m in 1:EKP.get_N_ens(ekp)
        output_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        diagnostics_folder_path = joinpath(output_path, "model_config")
        try
            compute_leaderboard(output_path, diagnostics_folder_path, spinup_months)
        catch e
            @error "Error in `analyze_iteration`" error = e
        end
    end
end

function plot_constrained_params_and_errors(output_dir, ekp, prior)
    dim_size = sum(length.(EKP.batch(prior)))
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500))
    for i in 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size + 1], ekp)
    EKP.Visualize.plot_error_over_time(fig[1, dim_size + 2], ekp)
    CairoMakie.save(joinpath(output_dir, "constrained_params_and_error.png"), fig)
    return nothing
end

# Process a single ensemble member's data into a vector
function process_member_data(simdir::SimDir, short_names, current_minibatch)
    # Define standard diagnostic fields to preprocess
    diagnostic_var_names = ["rsdt", "rsut", "rlut", "rsutcs", "rlutcs", "pr", "ts", "lwp", "clivi"]

    # Preprocess all diagnostic fields
    processed_data = Dict{String, Any}()
    for name in diagnostic_var_names
        processed_data[name] = preprocess_diagnostic_monthly_averages(simdir, name)
    end

    # Calculate derived fields
    processed_data["sw_cre"] = processed_data["rsut"] - processed_data["rsutcs"]
    processed_data["lw_cre"] = processed_data["rlut"] - processed_data["rlutcs"]
    processed_data["net_rad"] =
        processed_data["rlut"] + processed_data["rsut"] - processed_data["rsdt"] |> average_lat |> average_lon


    # Apply latitude window and rsdt weighting to IWP/LWP
    rsdt_mask = remake(processed_data["rsdt"], data = processed_data["rsdt"].data .> 0)
    rsdt_mask = window(rsdt_mask, "latitude"; left = -60, right = 60)
    processed_data["iwp"] = processed_data["clivi"]
    for field in ["lwp", "iwp"]
        processed_data[field] = window(processed_data[field], "latitude"; left = -60, right = 60)
        processed_data[field] = processed_data[field] * rsdt_mask
    end

    start_year = minimum(current_minibatch) + 2002
    year_range = (start_year):(start_year + length(current_minibatch) - 1)

    year_observations = map(year_range) do yr
        # Process seasonal data consistently
        seasonal_data = map(short_names) do short_name
            # TODO: Replace this with ClimaAnalysis/ClimaAnalysisExt
            year_averages = year_of_seasonal_averages(processed_data[short_name], yr)
            @show short_name, length(ClimaAnalysis.flatten(year_averages).data )
            ClimaAnalysis.flatten(year_averages).data 
        end
        return vcat(seasonal_data...)
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
    start_date = DateTime(monthly_avgs.attributes["start_date"], dateformat"yyyy-mm-ddTHH:MM:SS")
    monthly_avgs = window(monthly_avgs, "time"; left = DateTime(Year(start_date).value, 12, 1))
    global_mean = monthly_avgs |> average_lat |> average_lon |> average_time
    FT = monthly_avgs.data |> eltype

    # Replace NaNs with global mean
    monthly_avgs = ClimaAnalysis.replace(monthly_avgs, NaN => FT(NaNStatistics.nanmean(global_mean.data)))
    return monthly_avgs
end
