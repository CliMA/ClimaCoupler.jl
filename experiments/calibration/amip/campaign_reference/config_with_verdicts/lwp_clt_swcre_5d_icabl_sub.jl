# IC-mode ablation, SUB arm. Pair partner: lwp_clt_swcre_5d_icabl_amip.jl.
# Both arms are graded on IDENTICAL observations (Sep+Oct 2010 monthly
# means, covariance Sep-Oct 2006-2010, same 8 params and checked priors);
# they differ only in how members reach September:
#   SUB:  ERA5 initialization 2010-08-25 (son_ic), 7-day spinup.
#   AMIP: default IC from 2010-08-01, 1-month spinup.
# The posterior difference is therefore attributable to IC mode (atmosphere
# + land state), with season and averaging window held fixed - the
# confound in the relax vs relax_amip comparison.
#
# Predictions (falsifiable):
# 1. The arms reproduce the relax vs relax_amip split at matched data:
#    SUB eps_rel near 0.05, AMIP eps_rel near 0.09-0.10. If instead both
#    land together, the split was driven by season/window, not ICs.
# 2. q_liq: SUB near 1.9e-4, AMIP near 1.35e-4 (same direction as the
#    unmatched pair).
# 3. Residuals reach the same floors in both arms (the data is identical,
#    so differing fit quality would indicate reachability, not ICs).
#
# Verdict (run 2026-08-02, stopped converged after iteration 5; spreads
# ~0.05x prior, posterior unchanged from iteration 4):
# Posterior: q_liq 1.34e-4, rain_tau 1550, snow_tau 1895, cond_evap 117,
# subl_dep 333, margin 1.118, eps_rel 0.0655, steepness 0.522.
# Triangulation vs relax (ERA5, Oct only) and relax_amip (default IC,
# Aug-Oct):
# 1. q_liq: 1.34e-4 here vs 1.35e-4 in relax_amip, with matched tightened
#    priors. The apparent IC split (1.93e-4 vs 1.35e-4) was the broken
#    wide prior in relax, not initialization. q_liq is IC-insensitive.
#    (Prediction 2 falsified, for a good reason.)
# 2. eps_rel: 0.049 (ERA5, Oct) -> 0.0655 (ERA5, Sep+Oct) -> 0.096 (free,
#    Aug-Oct). Roughly one third of the subseasonal-to-AMIP shift is
#    window/season, two thirds is IC mode. Prediction 1 partially holds:
#    the split is real but not purely IC.
# 3. margin: 1.066 -> 1.118 -> 1.25, same one-third/two-thirds shape.
# 4. cond_evap is protocol-invariant (106/117/106) - the most
#    transferable parameter in the set.
# 5. rain_tau (1167 -> 1550 -> 1733) is contaminated: relax used the
#    sigma=700 prior, the other two sigma=400, so its window share is
#    inflated. Not attributable cleanly.
# 6. Prediction 3 PASS: identical floors reached (0.6-0.7 sigma).
# Implication for the ladder: eps_rel and the margin are confirmed
# environment-compensating (wide priors when transferring); q_liq and
# cond_evap transfer. Arm B (ERA5 Jul 1 long lead) would split the IC
# share into attractor drift vs cold start; decide after relax_amip
# finishes.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_3mo.yml",
)

sample_date_ranges =
    [(Dates.DateTime(2010, 9, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(
    pkgdir(ClimaCoupler),
    "amip_calibration_icabl_sub",
)

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

# lwp/swcre keep the 1.5x-scaled 10 degree floors; clt drops to 0.35
# (1.4x the 10d companion's 0.25) to make it informative.
const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [
    (
        short_names = ["lwp"],
        model_error_scale = 0.6,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["clt"],
        model_error_scale = 0.35,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["swcre"],
        model_error_scale = 0.75,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# checked_constrained_gaussian: EKP's constrained_gaussian silently ignores
# targets with magnitudes below ~1e-2 (see prior_tools.jl). Every earlier
# q_liq prior was in fact a unit normal over its bounds (mean = bounds
# center 5e-4), which is why iteration-1 q_liq always started there.
include(joinpath(@__DIR__, "..", "prior_tools.jl"))

const CALIBRATION_PRIORS = [
    checked_constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        # Tightened around the reproduced 2D answer (1.6-1.9e-4 across four
        # runs). The zonal cluster (2.2-2.55e-4) stays within ~1 sigma.
        1.8e-4, 0.6e-4, 0.0, 1e-3,
    ),
    checked_constrained_gaussian("rain_autoconversion_timescale", 1800, 400, 300, 3600),
    checked_constrained_gaussian("snow_autoconversion_timescale", 1800, 400, 300, 3600),
    checked_constrained_gaussian("condensation_evaporation_timescale", 150, 100, 20, 3600),
    checked_constrained_gaussian("sublimation_deposition_timescale", 500, 400, 20, 3600),
    checked_constrained_gaussian("cloud_fraction_floor_release_abs_margin", 1, 0.5, 0, 2),
    checked_constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    checked_constrained_gaussian("cloud_fraction_steepness_scale", 0.7, 0.2, 0.1, 2.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
