import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
ENV["CLIMACOMMS_DEVICE"] = "CUDA"
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))

include("cre_leaderboard.jl")
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
        plot_cre_leaderboard_from_iters(output_dir, SPINUP, iteration + 1)
        plot_constrained_params_and_errors(output_dir, ekp, prior)
    catch
    finally
        return nothing
    end
    return nothing
end

"""
    plot_constrained_params_and_errors(output_dir, ekp, prior)

Save a figure in `output_dir` of the contrained parameters and the error over
the number of iterations.
"""
function plot_constrained_params_and_errors(output_dir, ekp, prior)
    dim_size = prod(length.(EKP.batch(prior)))
    # Add one more for the error plot
    fig = CairoMakie.Figure(size = ((dim_size + 1) * 500, 500));
    for i in 1:dim_size
        EKP.Visualize.plot_ϕ_over_iters(fig[1, i], ekp, prior, i)
    end
    EKP.Visualize.plot_error_over_iters(fig[1, dim_size], ekp)
    CairoMakie.save(joinpath(output_dir, "constrained_params_and_error.png"), fig)
end
