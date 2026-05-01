# To test a calibration, set `export TEST_CALIBRATION=` in the terminal or
# `ENV["TEST_CALIBRATION"]=""` in the REPL. A test calibration differs from a
# full calibration in three ways: it runs for only a single iteration, uses only
# one prior, and emulates the production of the diagnostics. The diagnostics
# produced will have the correct metadata, but the values may be nonsensical.

# For the test calibration, we only care whether a single iteration can be
# completed. Since each iteration follows the same structure, successfully
# completing one iteration will catch most errors with the full calibration
# pipeline

# Define which coupler config to use
config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Calibrate only on Jan 1 2010
sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:6]

# On Derecho, it is preferable to save the calibration output to the scratch
# directory (e.g. "/glade/derecho/scratch")
output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration")

# spinup and sample_date_ranges are chosen to match the only dates
# available in the wxquest_initial_conditions artifact
const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    # Note: Pressure-level variables require model output with
    # pressure_coordinates: true in config
    short_names = ["ta", "hur"],
    minibatch_size = 1,
    n_iterations = 1,
    sample_date_ranges = [
        (Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:6
    ],
    extend = Dates.Month(1),
    spinup = Dates.Day(0),
    output_dir,
    rng_seed = 42,
)

# Used in generate_observations.jl and observation_map.jl
# Units: Pa (not hPa)
const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]

# To disable normalization, update generate_observations.jl to not apply the
# normalization. You may want to do the same in the observation map as well.
const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

const CALIBRATION_PRIORS =
    [PD.constrained_gaussian("precipitation_timescale", 1200, 300, 300, 2400)]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
