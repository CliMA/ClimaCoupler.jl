using ClimaAnalysis, Dates
import ClimaCalibrate
import ClimaCoupler
import JLD2
import EnsembleKalmanProcesses as EKP

function ClimaCalibrate.observation_map(iteration)
    observation_vec = JLD2.load_object(observation_path)
    single_member_dims = length(EKP.get_obs(observations))
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
    rsut = preprocess_monthly_averages(simdir, "rsut")
    rsutcs = preprocess_monthly_averages(simdir, "rsutcs")
    cre = rsutcs - rsut
    return vcat(vec(rsut.data), vec(cre.data))
end

# Preprocess monthly averages to the right dimensions and dates, remove NaNs
days = 86_400
spinup_time = 93days
function preprocess_monthly_averages(simdir, name)
    monthly_avgs = get_monthly_averages(simdir, name)
    monthly_avgs = ClimaAnalysis.shift_to_start_of_previous_month(monthly_avgs)
    # Remove spinup time
    monthly_avgs = window(monthly_avgs, "time"; left = spinup_time)
    global_mean = monthly_avgs |> average_lat |> average_lon |> average_time # average_time(average_lon(average_lat(monthly_avgs)))
    monthly_avgs = ClimaAnalysis.replace(monthly_avgs, NaN => first(global_mean.data))
    return monthly_avgs
end

get_monthly_averages(simdir, var_name) = get(simdir; short_name = var_name, reduction = "average", period = "1M")
