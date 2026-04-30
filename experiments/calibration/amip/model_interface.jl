import ClimaCoupler
import ClimaCoupler: CalibrationTools
import ClimaCalibrate
using Pkg

struct CouplerModelInterface <: ClimaCalibrate.AbstractModelInterface
    config::CalibrationTools.CalibrateConfig
end

function ClimaCalibrate.forward_model(interface::CouplerModelInterface, iter, member)
    Pkg.status()
    (; config) = interface
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(config.config_file)
    output_dir_root = config.output_dir

    (; sample_date_ranges, spinup, extend) = config
    start_date = first(sample_date_ranges[iter + 1]) - spinup
    end_date = last(sample_date_ranges[iter + 1]) + extend
    CalibrationTools.update_timespan!(config_dict, start_date, end_date)

    # Set member parameter file
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    CalibrationTools.add_parameter_filepath!(config_dict, sampled_parameter_file)

    # Set member output directory
    member_output_dir =
        ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir

    @info "Simulation dates" start_date end_date

    TEST_CALIBRATION = haskey(ENV, "TEST_CALIBRATION")
    if !TEST_CALIBRATION
        ClimaCoupler.SimCoordinator.setup_and_run(config_dict)
    else
        @info "Emulating diagnostics for test calibration"
        CalibrationTools.setup_and_emulate_diagnostics(config_dict)
    end
    @info "Completed member $member"
    return nothing
end

include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "observation_map.jl",
    ),
)

include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "post_analyze_iteration.jl",
    ),
)
