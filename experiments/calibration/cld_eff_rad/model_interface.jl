import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))

function ClimaCalibrate.forward_model(iter, member)
    redirect_stderr(stdout) # This should live in ClimaCalibrate.set_worker_logger
    config_file = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "cld_eff_rad", "model_config.yml")
    config_dict = get_coupler_config_dict(config_file)

    output_dir_root = config_dict["coupler_output_dir"]

    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    config_dict["coupler_toml"] = [sampled_parameter_file]
    # Set member output directory
    member_output_dir = ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir
    sim = try
        setup_and_run(config_dict)
    catch e
        @error "Forward model error" exception = e
        bt = catch_backtrace()
        println("Stacktrace:")
        display(stacktrace(bt))
        nothing
    end
    @info "Completed member $member"
    return sim
end
