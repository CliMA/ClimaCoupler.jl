using ClimaAnalysis, Dates
import ClimaCalibrate
import ClimaCoupler
import JLD2
import EnsembleKalmanProcesses as EKP
obs = JLD2.load_object("experiments/calibration/coarse_amip/observations.jld2")
const single_member_dims = length(EKP.get_obs(first(obs)))
include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/observation_utils.jl"))

function ClimaCalibrate.observation_map(iteration)
    G_ensemble = Array{Float64}(undef, single_member_dims, ensemble_size)
    for m in 1:ensemble_size
        @show m
        member_path = ClimaCalibrate.path_to_ensemble_member(output_dir, iteration, m)
        simdir_path = joinpath(member_path, "model_config/output_active")
        if isdir(simdir_path)
            simdir = SimDir(simdir_path)
            G_ensemble[:, m] .= process_member_data(simdir)
        else
            @info "No data found for member $m."
            G_ensemble[:, m] .= NaN
        end
    end
    return G_ensemble
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
    # Ensure we are splitting evenly across seasons
    @assert all(map(x -> length(times(x)) == 3, seasons))
    seasonal_avgs = average_time.(seasons)

    downsampled_seasonal_avg_arrays = downsample.(seasonal_avgs, 3)
    return vcat(vec.(downsampled_seasonal_avg_arrays)...)
    # return vectorize_nyears_of_seasonal_outputvars(seasonal_avgs, 1)
end

function process_member_data(simdir::SimDir)
    isempty(simdir) && return fill!(zeros(single_member_length), NaN)

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
