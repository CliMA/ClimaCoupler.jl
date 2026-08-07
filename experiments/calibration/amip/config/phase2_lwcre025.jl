# PHASE-2 rerun with a tightened lwcre noise floor (0.5 -> 0.25).
#
# The first phase2 run finished with lwcre at 0.33 sigma from iteration 1
# through iteration 8 while every other observable sat at its floor
# (prediction 6 of phase2.jl said "lwcre stays near 1 sigma" and failed).
# A residual far below the floor means the assumed structural error is
# about twice the real one, so lwcre contributed almost no constraint.
# This config halves the lwcre floor and changes nothing else.
#
# Predictions (falsifiable):
# 1. lwcre lands near its new floor (0.6 to 1.1 sigma), no longer far below.
# 2. The extra lwcre information tightens the detrainment parameters:
#    final spreads for detr_massflux_vertdiv_coeff and detr_buoy_coeff
#    come in below the first run's 0.58 and 0.73 x prior.
# 3. lwp, pr, swcre stay at their floors (0.8 to 1.1 sigma) and q_liq
#    stays in the 2 to 3e-4 cluster; no COLLAPSED flag.
# If lwcre instead pins the fit and pushes other observables off their
# floors, 0.25 is too tight and the right floor is between 0.25 and 0.5.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:9]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_phase2_lwcre025")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "pr", "swcre", "lwcre"],
    minibatch_size = 1,
    n_iterations = 8,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Same floors as phase2.jl except lwcre, which gets its own group at 0.25.
const OBS_NOISE_GROUPS = [
    (short_names = ["lwp", "pr"], model_error_scale = 0.2),
    (short_names = ["swcre"], model_error_scale = 0.5),
    (short_names = ["lwcre"], model_error_scale = 0.25),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        3e-4, 1.5e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1450, 500, 300, 3600),
    PD.constrained_gaussian("detr_massflux_vertdiv_coeff", 0.35, 0.15, 0.0, 0.8),
    PD.constrained_gaussian("detr_buoy_coeff", 1.0, 0.6, 0.05, 5.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
