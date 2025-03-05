using ClimaAnalysis, Dates
import ClimaCalibrate
import ClimaCoupler
import JLD2
import EnsembleKalmanProcesses as EKP
include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/observation_utils.jl"))

function ClimaCalibrate.observation_map(iteration)
    observation_vec = JLD2.load_object(observation_path)
    single_member_dims = length(EKP.get_obs(first(observation_vec)))
    G_ensemble = Array{Float64}(undef, single_member_dims, ensemble_size)
    for m in 1:ensemble_size
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, "model_config/output_active")
        try
            G_ensemble[:, m] .= process_member_data(SimDir(simdir_path))

        catch e
            @warn "Error processing member $m, filling observation map entry with NaNs" exception = e
            G_ensemble[:, m] .= NaN
        end
    end
    return G_ensemble
end

function process_member_data(simdir::SimDir)
    pressure = get_monthly_averages(simdir, "pfull")

    rsdt_full = get_monthly_averages(simdir, "rsdt")
    rsut_full = get_monthly_averages(simdir, "rsut")
    rlut_full = get_monthly_averages(simdir, "rlut")
    
    year_net_radiation = (rlut_full + rsut_full - rsdt_full) |> average_lat |> average_lon |> average_time

    rsut = process_outputvar(simdir, "rsut")
    rlut = process_outputvar(simdir, "rlut")
    rsutcs = process_outputvar(simdir, "rsutcs")
    rlutcs = process_outputvar(simdir, "rlutcs")
    cre = rsut + rlut - rsutcs - rlutcs

    pr = process_outputvar(simdir, "pr")
    # shf = process_outputvar(simdir, "shf")
    ts = process_outputvar(simdir, "ts")

    ta = process_outputvar(simdir, "ta")
    hur = process_outputvar(simdir, "hur")
    hus = process_outputvar(simdir, "hus")
    # clw = get_seasonal_averages(simdir, "clw")
    # cli = get_seasonal_averages(simdir, "cli")

    return vcat(year_net_radiation.data, rsut, rlut, cre, pr, ts)#, ta, hur, hus)
end

function process_outputvar(simdir, name)
    days = 86_400

    monthly_avgs = get_monthly_averages(simdir, name)
    # Preprocess to match observations
    if has_altitude(monthly_avgs)
        pressure = get_monthly_averages(simdir, "pfull")
        monthly_avgs = ClimaAnalysis.Atmos.to_pressure_coordinates(monthly_avgs, pressure)
        monthly_avgs = limit_pressure_dim_to_era5_range(monthly_avgs)
    end
    monthly_avgs = ClimaAnalysis.replace(monthly_avgs, missing => 0.0, NaN => 0.0)
    monthly_avgs = ClimaAnalysis.shift_to_start_of_previous_month(monthly_avgs)

    # Cut off first 3 months
    single_year = window(monthly_avgs, "time"; left = 92days)
    seasons = split_by_season_across_time(single_year)
    # Ensure each season has three months
    @assert all(map(x -> length(times(x)) == 3, seasons))
    seasonal_avgs = average_time.(seasons)

    downsampled_seasonal_avg_arrays = downsample.(seasonal_avgs, 3)
    return vcat(vec.(downsampled_seasonal_avg_arrays)...)
end
