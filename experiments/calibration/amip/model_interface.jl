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
    # Use sample_date_ranges[iter], NOT [iter + 1]: at calibration loop-iteration
    # `iter` the EKP update compares this member's G against observation `iter`
    # (the minibatcher yields batch [iter], advancing after each update). Running
    # sample_date_ranges[iter + 1] here compares model output for one sample
    # against the observation for the previous sample — a harmless no-op when all
    # sample_date_ranges are identical (as in fixed-target runs), but it pairs the
    # wrong year with the wrong observation once the samples differ (e.g. a
    # multi-year minibatch). Indexing with `iter` keeps G and obs on the same sample.
    start_date = first(sample_date_ranges[iter]) - spinup
    end_date = last(sample_date_ranges[iter]) + extend
    CalibrationTools.update_timespan!(config_dict, start_date, end_date)

    # Set member parameter file. Sampled `<base>_E<index>` scalars calibrate
    # single elements of the vector parameter <base>: they are spliced into the
    # base vector from the run's own coupler_toml files before ClimaParams sees
    # them (any name/index mismatch errors, failing the member). Without `_E#`
    # entries the sampled file is used as is.
    sampled_parameter_file = ClimaCalibrate.parameter_path(output_dir_root, iter, member)
    sampled_parameter_file = CalibrationTools.write_spliced_parameter_file(
        sampled_parameter_file,
        CalibrationTools.parameter_dict(config_dict),
        joinpath(dirname(sampled_parameter_file), "parameters_spliced.toml"),
    )
    CalibrationTools.add_parameter_filepath!(config_dict, sampled_parameter_file)

    # Set member output directory
    member_output_dir =
        ClimaCalibrate.path_to_ensemble_member(output_dir_root, iter, member)
    config_dict["coupler_output_dir"] = member_output_dir
    config_dict["detect_restart_files"] = true
    config_dict["checkpoint_dt"] = "10days"

    ClimaCoupler.Input.update_t_start_for_restarts!(config_dict)

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
