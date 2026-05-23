import ClimaCoupler
import ClimaCoupler: CalibrationTools
import ClimaCalibrate

include(joinpath(pkgdir(ClimaCoupler), "experiments", "AMIP", "code_loading.jl"))

import Pkg
import Statistics

"""
    CouplerModelInterface <: ClimaCalibrate.AbstractModelInterface

A model interface struct for running the AMIP calibration pipeline.

See the ClimaCalibrate.jl documentation for the methods that
`CouplerModelInterface` should implement.
"""
struct CouplerModelInterface <: ClimaCalibrate.AbstractModelInterface
    config::CalibrationTools.CalibrateConfig
end

"""
    ClimaCalibrate.forward_model(interface::CouplerModelInterface, iter, member)

Run a coupled model simulation.

This function may be called in parallel depending on the ClimaCalibrate backend
used.
"""
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
    # Ensure the simulation restarts automatically
    config_dict["detect_restart_file"] = true
    config_dict["output_dir_style"] = "activelink"
    config_dict["checkpoint_dt"] = "1days"

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

"""
    ClimaCalibrate.model_interface_filepath(::CouplerModelInterface)

Return a filepath to the definition of the `CouplerModelInterface` struct and
all its associated methods.

This is required to use the `ClimaCalibrate.HPCBackend`s.
"""
function ClimaCalibrate.model_interface_filepath(::CouplerModelInterface)
    return @__FILE__
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
