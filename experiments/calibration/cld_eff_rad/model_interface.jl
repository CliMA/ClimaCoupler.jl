import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
import LinearAlgebra

ENV["CLIMACOMMS_DEVICE"] = "CUDA"
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))

include("cre_leaderboard.jl")
include("plot_from_EKP.jl")
const SPINUP = 1

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
    sim = setup_and_run(config_dict)
    @info "Completed member $member"
    return sim
end

function ClimaCalibrate.analyze_iteration(ekp, g_ensemble, prior, output_dir, iteration)
    # Use iteration + 1 because iteration numbers are 0-indexed
    try
        # TODO: Make a plots directory if it doesn't exist and use that instead of output_dir
        plot_cre_leaderboard_from_iters(output_dir, SPINUP, iteration + 1)
        plot_constrained_params_and_errors(output_dir, ekp, prior)
    catch e
        @info e
    finally
        return nothing
    end
    return nothing
end
