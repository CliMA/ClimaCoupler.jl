# To test a calibration, set `export TEST_CALIBRATION=` in the terminal or
# `ENV["TEST_CALIBRATION"]=""` in the REPL. A test calibration differs from a
# full calibration in three ways: it runs for only a single iteration, uses only
# one prior, and emulates the production of the diagnostics. The diagnostics
# produced will have the correct metadata, but the values may be nonsensical.

# For the test calibration, we only care whether a single iteration can be
# completed. Since each iteration follows the same structure, successfully
# completing one iteration will catch most errors with the full calibration
# pipeline

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    # Coarse h_elem for the CI smoke test only
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config",
        "amip_configs",
        "amip_calibration_pipeline_test.yml",
    ),
    # Pressure-level fields need model output with pressure_coordinates: true in the yml.
    observables = [
        (name = "ta", reduce = zonal_average, model_error_fraction = 0.05),
        (name = "hur", reduce = zonal_average, model_error_fraction = 0.05),
    ],
    pressure_levels = 100.0 .* [200.0, 500.0, 850.0],  # Pa
    n_iterations = 1,
    # October 2010 has initial conditions in the wxquest_initial_conditions artifact. The
    # noise covariance is the spread over ten Octobers; identical dates would give none.
    sample_date_ranges = [(DateTime(2010, 10, 1), DateTime(2010, 10, 1))],
    covariance_date_ranges = [(DateTime(y, 10, 1), DateTime(y, 10, 1)) for y in 2001:2010],
    extend = Month(1),
    spinup = Day(0),
    # On Derecho, it is preferable to save the calibration output to scratch.
    output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration"),
    priors = EKP.combine_distributions([
        PD.constrained_gaussian("precipitation_timescale", 1200, 300, 300, 2400),
    ]),
    worker_walltime = 60,
)
