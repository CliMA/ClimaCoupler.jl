# 10 degree lwp+pr with a correlated noise floor AND a doubled floor (0.4).
#
# NOT YET LAUNCHED. Prepared 2026-07-29 after the L = 800 km correlation-only
# rerun (lwp_pr_2d10_corr.jl) tracked the diagonal run's collapse trajectory
# (0.02x prior spread at iteration 2 with a 2 to 3 sigma residual). Two
# reasons the correlation alone is not enough:
# 1. Parameter sensitivities are spatially smooth, and the exponential
#    kernel only inflates smooth directions by about 2 to 3x at this block
#    spacing (1100 km blocks, L = 800 km, neighbor correlation 0.25). The
#    contraction slows about 25%, not the ~10x needed.
# 2. The 2 sigma stuck residual means the model cannot fit 2-D lwp and pr
#    structure to 20% of the field mean. Zonal averaging cancels
#    compensating regional biases; a 2-D grid keeps them. The honest 2-D
#    structural floor is about twice the zonal one.
# This config keeps L = 800 km and raises model_error_scale 0.2 -> 0.4.
#
# Predictions (falsifiable):
# 1. No COLLAPSED flag; final spread above 0.05x prior.
# 2. Residuals reach the (new) floor, about 1 to 1.3 sigma.
# 3. q_liq within half a prior sigma of 2.8e-4.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_pr_2d10_corr04")

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

const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [(
    short_names = ["lwp", "pr"],
    model_error_scale = 0.4,
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
