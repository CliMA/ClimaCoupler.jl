import ClimaCoupler
import ClimaCalibrate
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))

function ClimaCalibrate.forward_model(iter, member)
    config_file =
        joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "model_config.yml")
    config_dict = Input.get_coupler_config_dict(config_file)

    # Run for a shorter time if SHORT_RUN is set
    if SHORT_RUN
        config_dict["t_end"] = "480secs"
        config_dict["dt_rad"] = "480secs"
        map(diag -> diag["period"] = "480secs", config_dict["extra_atmos_diagnostics"])
    end

    output_dir_root = config_dict["coupler_output_dir"]
    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    config_dict["coupler_toml"] = [sampled_parameter_file]
    # Set member output directory
    member_output_dir =
        ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir
    return setup_and_run(config_dict)
end
