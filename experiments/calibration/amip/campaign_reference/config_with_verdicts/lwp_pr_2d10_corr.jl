# 10 degree lwp+pr rerun with a spatially correlated noise floor.
#
# The diagonal-floor run (lwp_pr_2d10.jl) collapsed: parameter spread reached
# 0.0x prior by iteration 5 while the residual was still 2.1-2.4 sigma, and
# the 15 degree companion showed the same signature more mildly. The diagonal
# floor treats all ~650 grid cells as independent constraints even though
# monthly means decorrelate over ~800 km. This config keeps everything else
# identical and adds decorrelation_length = 800 km to the floor (see
# correlated_noise.jl), so a cluster of cells inside one correlation patch
# counts as roughly one constraint.
#
# Predictions (falsifiable):
# 1. No COLLAPSED steering flag through 6 iterations; final spread stays
#    above 0.05x prior for both parameters.
# 2. q_liq lands near the zonal answer (2.8e-4) within half a prior sigma,
#    with rain_tau near 1450.
# 3. Residuals settle toward the floor (below ~1.5 sigma) instead of
#    freezing at 2 sigma with drift larger than the spread.
# If the run still collapses, the over-information is not (only) spatial
# correlation and hypothesis B needs revision.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_pr_2d10_corr")

const COARSEN_FACTOR = 4

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "pr"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Identical floor to the diagonal run; the only change is the correlation.
const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [(
    short_names = ["lwp", "pr"],
    model_error_scale = 0.2,
    decorrelation_length = OBS_DECORRELATION_LENGTH,
)]

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
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
