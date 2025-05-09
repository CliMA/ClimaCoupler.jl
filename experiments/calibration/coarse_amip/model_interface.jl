import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
import JLD2
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))

function ClimaCalibrate.forward_model(iter, member)

    config_file = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "coarse_amip", "model_config.yml")
    config_dict = get_coupler_config_dict(config_file)

    output_dir_root = config_dict["coupler_output_dir"]
    eki = JLD2.load_object(ClimaCalibrate.ekp_path(output_dir_root, iter))
    minibatch = EKP.get_current_minibatch(eki)
    config_dict["start_date"] = minibatch_to_start_date(minibatch)
    @info "Current minibatch: $minibatch"
    @info "Current start date: $(minibatch_to_start_date(minibatch))"
    spinup_days = 92
    nyears = length(minibatch)
    t_end_days = spinup_days + 365 * nyears
    config_dict["t_end"] = "$(t_end_days)days"

    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    config_dict["coupler_toml"] = [sampled_parameter_file]
    # Set member output directory
    member_output_dir = ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir
    @show config_dict
    sim = try
        setup_and_run(config_dict)
    catch e
        @info e
        println(catch_backtrace())
        rethrow(e)
    end

    @info "Completed member $member"
    return sim
end

function minibatch_to_start_date(batch)
    start_year = minimum(batch) + 2001
    @assert start_year >= 2002
    return "$(start_year)0901"
end
