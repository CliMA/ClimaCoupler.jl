import ClimaCoupler
import ClimaCalibrate

include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))

function ClimaCalibrate.forward_model(iter, member)
    config_file = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip_config.yml")
    config_dict = get_coupler_config_dict(config_file)

    output_dir_root = config_dict["coupler_output_dir"]
    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    if haskey(config_dict, "toml")
        append!(config_dict, "toml", sampled_parameter_file)
    else
        config_dict["toml"] = sampled_parameter_file
    end

    # Set member output directory
    member_output_dir = ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir
    # TODO: disable default diagnostics
    return setup_and_run(config_dict)
end
