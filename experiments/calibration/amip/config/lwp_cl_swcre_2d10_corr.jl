# 10 degree lwp+cl+swcre rerun with a spatially correlated noise floor.
#
# The diagonal-floor run (lwp_cl_swcre_2d10.jl) contracted to 0.01-0.03x
# prior spread by iteration 4, far tighter than the healthy zonal companion
# at the same point, with drift up to 2.6x spread per iteration. Same
# over-information mechanism as the lwp_pr 2-D pair: the diagonal floor
# counts every grid cell as independent. This config adds
# decorrelation_length = 800 km to both noise groups (see
# correlated_noise.jl) and changes nothing else.
#
# Predictions (falsifiable):
# 1. Final spread stays above 0.1x prior for all five parameters (the zonal
#    run holds 0.6-0.97x at iteration 1; the diagonal 2d10 run fell an order
#    of magnitude below it).
# 2. Posterior means agree with the zonal run within one final spread for
#    the parameters both constrain (q_liq first).
# 3. cl and swcre still fit to their floors (about 1 sigma), showing the
#    correlation removes fake information, not real signal.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir =
    joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_cl_swcre_2d10_corr")

const COARSEN_FACTOR = 4

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "cl", "swcre"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Identical floors to the diagonal run; the only change is the correlation.
const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [
    (
        short_names = ["lwp"],
        model_error_scale = 0.2,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["cl", "swcre"],
        model_error_scale = 0.5,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        5e-4, 3e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1800, 700, 300, 3600),
    PD.constrained_gaussian("snow_autoconversion_timescale", 1800, 700, 300, 3600),
    PD.constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    PD.constrained_gaussian("cloud_fraction_steepness_scale", 0.7, 0.2, 0.1, 2.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
