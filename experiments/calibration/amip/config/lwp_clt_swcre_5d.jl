# Cloud calibration on a 5 degree grid: lwp + clt + swcre.
#
# Two changes from lwp_cl_swcre_2d10.jl:
# 1. clt (CALIPSO/CloudSat total column cloud cover, type "any") replaces the
#    level-resolved cl. The model side is the clt diagnostic (shortwave McICA
#    cloud cover). Both are detection or overlap definitions of "cloud
#    anywhere in the column"; the CALIPSO detection threshold makes the obs
#    an undercount, like cl.
# 2. COARSEN_FACTOR = 2 (5 degrees, 72x36) instead of 4. The coarse 10 degree
#    blocks merged coastal stratocumulus decks (Peru, Namibia, California)
#    with land-masked cells; 5 degrees resolves them.
#
# A 5 degree grid has 4x the constraints of 10 degrees, so the floors scale
# up by 1.5x from the 10 degree values the corr04 run validates
# (lwp_pr_2d10_corr04.jl). FINALIZE THE FLOORS AGAINST THE corr04 VERDICT
# BEFORE LAUNCHING. The correlated floor stays on: at 5 degree spacing
# (555 km) the L = 800 km kernel correlates neighbors at 0.5, unlike at
# 10 degrees where it measurably did nothing.
#
# Predictions (falsifiable):
# 1. No COLLAPSED flag; final spread above 0.05x prior for all parameters.
# 2. Residuals reach their floors (about 1 to 1.5 sigma) instead of
#    freezing above 2 sigma.
# 3. q_liq lands in the healthy-run cluster (2 to 3e-4); the cloud-fraction
#    parameters move clt as they moved cl (posterior eps_rel near 0.04 to
#    0.06).
# 4. The Sc deck blocks appear as resolved features in the bias maps rather
#    than merged coastal blocks.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_clt_swcre_5d")

const COARSEN_FACTOR = 2

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

# 1.5x the corr04-validated 10 degree floors (see header).
const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [
    (
        short_names = ["lwp"],
        model_error_scale = 0.6,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["clt", "swcre"],
        model_error_scale = 0.75,
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
