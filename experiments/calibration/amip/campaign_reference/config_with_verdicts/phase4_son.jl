# PHASE-4-SON: phase4's design on the SON season (Sep+Oct+Nov).
#
# STAGED; LAUNCH GATES: (1) the Aug-25 ICs finish downloading and
# preprocessing into son_ic (get_aug25_ics.sh), (2) phase4 OND validates
# the rotating 3-month minibatch machinery (drift below 1x spread/iter
# from iteration 3, no collapse).
#
# Each member runs Aug 25 to Dec 1 (~98 days: 7 day spinup, then
# Sep/Oct/Nov monthly means). Advantages over OND: CALIPSO has complete
# SON for ALL five years (OND lost 2009 to the missing December), so the
# covariance uses 5 realizations instead of 4; and members are ~20%
# cheaper (98 vs 122 days). Compared against phase4's OND posterior, the
# SON run replicates the seasonal-robustness test with a different month
# mix (Sep in, Dec out).
#
# Predictions (falsifiable):
# 1. No COLLAPSED flag; residuals reach their floors as in phase3.
# 2. The rotating target does not destabilize the update: parameter drift
#    stays below 1x spread/iter from iteration 3 on.
# 3. q_liq lands in the 2.2-2.5e-4 cluster AND within one final spread of
#    phase4 OND's answer. Agreement across two independent season mixes
#    is the strongest seasonal-robustness evidence this campaign can
#    produce.
# 4. September (the month nearest its initial condition, 7 day spinup)
#    shows no elevated residual relative to Oct/Nov in the final
#    iteration; a hot September means the spinup is too short.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_3mo.yml",
)

# SON of one year per iteration; CALIPSO covers SON for all five years.
const SAMPLE_YEARS = [2010, 2009, 2008, 2007, 2006, 2010]

sample_date_ranges =
    [(Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1)) for y in SAMPLE_YEARS]

const COVARIANCE_DATE_RANGES =
    [(Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1)) for y in 2006:2010]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_phase4_son")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "pr", "swcre", "lwcre", "cl"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Floors carried over from phase3 (validated per observable).
const OBS_NOISE_GROUPS = [
    (short_names = ["lwp", "pr"], model_error_scale = 0.2),
    (short_names = ["swcre"], model_error_scale = 0.5),
    (short_names = ["lwcre"], model_error_scale = 0.25),
    (short_names = ["cl"], model_error_scale = 0.5),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# Same priors as phase3: the time extension should re-derive the answer,
# not inherit it.
const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        3e-4, 1.5e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1450, 500, 300, 3600),
    PD.constrained_gaussian("detr_massflux_vertdiv_coeff", 0.35, 0.15, 0.0, 0.8),
    PD.constrained_gaussian("detr_buoy_coeff", 1.0, 0.6, 0.05, 5.0),
    PD.constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    PD.constrained_gaussian("cloud_fraction_steepness_scale", 0.7, 0.2, 0.1, 2.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
