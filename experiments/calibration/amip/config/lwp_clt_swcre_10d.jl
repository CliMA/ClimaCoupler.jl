# Cloud calibration on a 10 degree grid: lwp + clt + swcre.
#
# Companion to lwp_clt_swcre_5d.jl: same observables and parameters, coarser
# grid, floors at the 10 degree base values (the 5 degree config carries
# 1.5x these, sized to its 4x constraint count). Running the pair tests the
# floor-scaling rule directly: if the accounting is right, both runs stay
# healthy and agree; if the 5 degree run collapses while this one holds,
# the scaling is too small.
#
# clt is CALIPSO/CloudSat total column cloud cover (type "any") against the
# model clt diagnostic (shortwave McICA cloud cover).
#
# Predictions (falsifiable):
# 1. No COLLAPSED flag; final spread above 0.05x prior for all parameters.
# 2. Residuals reach their floors (about 1 to 1.5 sigma).
# 3. q_liq in the healthy-run cluster (2 to 3e-4); posteriors agree with the
#    5 degree run within one final spread.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_clt_swcre_10d")

const COARSEN_FACTOR = 4

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "clt", "swcre"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [
    (
        short_names = ["lwp"],
        model_error_scale = 0.4,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["clt", "swcre"],
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
