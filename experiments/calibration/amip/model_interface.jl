import ClimaCoupler
import ClimaCalibrate
using Pkg

# Include run_calibration.jl only if CALIBRATE_CONFIG is not already defined
# This allows other pipelines to include their own run_calibration.jl first
if !@isdefined(CALIBRATE_CONFIG)
    include(joinpath(@__DIR__, "run_calibration.jl"))
end

if !TEST_CALIBRATION
    ENV["CLIMACOMMS_DEVICE"] = "CUDA"
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
end

function ClimaCalibrate.forward_model(iter, member)
    Pkg.status()
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(CALIBRATE_CONFIG.config_file)
    output_dir_root = CALIBRATE_CONFIG.output_dir

    start_date =
        first(CALIBRATE_CONFIG.sample_date_ranges[iter + 1]) - CALIBRATE_CONFIG.spinup
    end_date = last(CALIBRATE_CONFIG.sample_date_ranges[iter + 1]) + CALIBRATE_CONFIG.extend
    CalibrationTools.update_timespan!(config_dict, start_date, end_date)

    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    CalibrationTools.add_parameter_filepath!(config_dict, sampled_parameter_file)

    # Set member output directory
    member_output_dir =
        ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir

    @info "Simulation dates" start_date end_date

    if !TEST_CALIBRATION
        ClimaCoupler.SimCoordinator.setup_and_run(config_dict)
    else
        @info "Emulating diagnostics for test calibration"
        CalibrationTools.setup_and_emulate_diagnostics(config_dict)
    end
    @info "Completed member $member"
    return nothing
end
