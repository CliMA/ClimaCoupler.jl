# A worked zonal-mean example: liquid water path and the two cloud radiative
# effects. Every observable here is reduced to a zonal mean by the shared
# preprocessing.

# Define which coupler file to use
config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Calibrate on October 2010, the October with initial conditions in the
# wxquest_initial_conditions artifact. The noise covariance is the interannual
# spread of the observations, so it is estimated over ten Octobers; identical
# dates would give it no spread at all.
sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:6]
covariance_date_ranges =
    [(Dates.DateTime(year, 10, 1), Dates.DateTime(year, 10, 1)) for year in 2001:2010]

# On Derecho, it is preferable to save the calibration output to the scratch
# directory (e.g. "/glade/derecho/scratch")
output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_pressure_levels")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "swcre", "lwcre"],
    minibatch_size = 1,
    n_iterations = 5,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# No pressure-level variables here, so this selects nothing and the step is a no-op.
const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]

# To disable normalization, update generate_observations.jl to not apply the
# normalization. You may want to do the same in the observation map as well.
const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# How much of the model-data difference to treat as irreducible model error, as a
# variance, one value per observable in the order of `short_names`.
const NOISE_BETA = [0.10, 0.074, 0.074]

const CALIBRATION_PRIORS = [
    PD.constrained_gaussian("cloud_fraction_eps_rel", 0.05, 0.02, 0.001, 0.2),
    PD.constrained_gaussian("entr_coeff", 0.3, 0.15, 0.02, 1.0),
    PD.constrained_gaussian("cloud_liquid_rain_collision_efficiency", 0.8, 0.15, 0.1, 1.0),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1800, 300, 300, 3600),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
