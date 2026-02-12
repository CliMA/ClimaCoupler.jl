ENV["CLIMACOMMS_DEVICE"] = "CUDA"
ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
import ClimaCoupler
import ClimaCalibrate
import CUDA
import EnsembleKalmanProcesses as EKP
using Dates: Date, Second
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))

# Note: api.jl and CALIBRATE_CONFIG are loaded/shared via @everywhere in run_calibration.jl
using Pkg
function ClimaCalibrate.forward_model(iter, member)
    Pkg.status()
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(CALIBRATE_CONFIG.config_file)
    output_dir_root = CALIBRATE_CONFIG.output_dir
    # Use the first (and typically only) sample date range for all iterations
    # EKP updates parameters across iterations, but uses the same observation period
    sample_date_range = first(CALIBRATE_CONFIG.sample_date_ranges)
    start_date = first(sample_date_range) - CALIBRATE_CONFIG.spinup
    start_date_str = replace(string(Date(start_date)), "-" => "")
    end_date = last(sample_date_range) + CALIBRATE_CONFIG.extend
    sim_length = Second(end_date - start_date)

    config_dict["start_date"] = start_date_str
    config_dict["t_end"] = "$(sim_length.value)secs"
    config_dict["checkpoint_dt"] = "900days"
    config_dict["dt"] = "120secs"
    config_dict["dt_cpl"] = "120secs"

    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    if haskey(config_dict, "coupler_toml")
        config_dict["coupler_toml"] =
            [config_dict["coupler_toml"]..., sampled_parameter_file]
    else
        config_dict["coupler_toml"] = [sampled_parameter_file]
    end
    # Set member output directory
    member_output_dir =
        ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir

    @info "Simulation dates" start_date end_date
    setup_and_run(config_dict)
    @info "Completed member $member"
    return nothing
end
