# PHASE-5 combined calibration: phase3's 6 parameters x 6 observables
# (phase3's five + clt), zonal, on phase4's 3-month multi-year minibatch.
#
# Launch gates, all passed 2026-07-30:
#   1. phase3 finished healthy (02:10 MDT; reference posterior in
#      identifiability_map.md).
#   2. The clt-floor rerun validated clt 0.25 at 10 degrees (clt 0.75
#      sigma, steering, cloud-fraction spreads tightened); adopted here
#      for the zonal clt floor.
#   3. The Aug-25 ICs are preprocessed into son_ic.
# Phase4 (SON without clt) is SKIPPED by decision (2026-07-30): it stays
# staged as the on-demand ablation if this run's q_liq lands outside the
# 2.2-2.5e-4 cluster and the season-vs-clt attribution matters. The
# rotating-minibatch machinery is self-validated here by the iteration-3
# kill rule (drift below 1x spread/iter, no collapse).
#
# Predictions (falsifiable):
# 1. No COLLAPSED flag; residuals reach their floors; the rotating target
#    keeps drift below 1x spread/iter from iteration 3.
# 2. q_liq lands in the 2.2-2.55e-4 zonal cluster. If it lands near the
#    2D answer (1.6-1.9e-4) instead, run the staged phase4 ablation.
# 3. eps_rel arbitration: phase3 gave 0.0335, the clt reruns 0.058. If
#    clt information dominates, expect 0.04-0.06; if the detrainment
#    trade dominates, expect ~0.033. Either way it must stay separable
#    from the detrainment pair.
# 4. clt sits near its floor (0.7-1.3 sigma) without pushing cl off its
#    own floor.
# 5. rain_tau finally tightens: final spread below phase3's 0.74x prior
#    (3 monthly means x rotating years is 3x phase3's data).
# 6. December... does not apply (SON); September carries the 7-day
#    spinup: its per-month residual must not exceed Oct/Nov's by more
#    than 2x (else the spinup is too short and phase6 moves to Aug-18
#    ICs).
#
# Verdict (judged 2026-07-31 on iterations 1-5; iteration 6 sims rerunning
# after a walltime kill and can only refine, not update, these numbers):
# 1. PARTIAL. No collapse; all six observables at or below floor by iter 5
#    (lwp 0.92, pr 1.05, cl 0.86, swcre 0.67, lwcre 0.68, clt 0.33 sigma).
#    Drift exceeded 1x spread/iter at iters 3 AND 4 (rain_tau 1.12,
#    eps_rel 1.21, steepness 1.13) and settled only at iter 5 (max 0.8).
#    The year rotation kicks the target harder than predicted; parameters
#    oscillate mildly rather than run away.
# 2. PASS. q_liq 2.36e-4, stable (2.63 -> 2.41 -> 2.36), inside the
#    zonal cluster. The phase4 SON ablation stays parked.
# 3. Leans zonal. eps_rel 0.0392, separable from the detrainment pair
#    throughout. The grid tension narrowed but persists: this run 0.039
#    vs the 5d relax run 0.0488 (was 0.0335 vs 0.058).
# 4. HALF. cl held its floor (0.86 sigma) but clt fell to 0.33 sigma,
#    below the 0.7-1.3 band. Second independent confirmation the zonal
#    clt floor 0.25 is ~2x generous; phase6 drops it to ~0.12-0.15.
# 5. PASS. rain_tau spread 0.33x prior (< 0.74x), mean 1429.
# 6. PASS, inverted. September is the BEST month, so the 7-day spinup is
#    adequate. Residuals instead grow with lead time, sharpest in pr:
#    Sep 0.72 -> Oct 0.85 -> Nov 1.43 sigma (same gradient at iters 3-5,
#    also visible in lwcre 0.55 -> 0.85). November pr holds learnable
#    signal; candidate phase6 lever alongside the clt floor.
# Posterior after iter 5: q_liq 2.36e-4, rain_tau 1429, dmfvd 0.563,
# detr_buoy 1.36, eps_rel 0.0392, steepness 0.616.
# Final posterior (run completed 2026-07-31, 6 iterations): q_liq 2.33e-4,
# rain_tau 1385, dmfvd 0.561, detr_buoy 1.65 (back at phase3's 1.61),
# eps_rel 0.0365 (further toward the zonal 0.0335; verdict 3 firms up),
# steepness 0.589. Iter-6 residuals: lwp 1.01, pr 1.25, cl 0.84 at floor;
# swcre 0.65, lwcre 0.67, clt 0.32 below. All verdicts stand.
#
# Prior caveat (found 2026-07-31): the q_liq prior below never took
# effect. EKP's constrained_gaussian silently returns a unit normal over
# the bounds for small-magnitude targets (see prior_tools.jl), so this
# run's effective q_liq prior had mean 7.5e-4 (bounds center) and was
# wide, not (3e-4, 1.5e-4). The data pulled q_liq DOWN to 2.36e-4 against
# that higher prior mean, which strengthens, not weakens, prediction 2's
# conclusion.
#
# Sampling is phase4 SON's: ~98-day members (Aug 25 IC, 7 day spinup,
# Sep+Oct+Nov monthly means), target year rotating over all five years
# (CALIPSO SON is complete 2006-2010, so the covariance keeps 5
# realizations). clt (CALIPSO column cover vs model McICA clt)
# complements cl: cl constrains the vertical profile at 3 altitudes, clt
# the total-column overlap.
#
# Predictions (falsifiable): filled in at launch, after the gates. Expect
# phase4's predictions to carry over, plus: clt sits at its floor without
# pushing cl off its own, and the cloud-fraction posteriors tighten below
# phase3's final spreads.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_3mo.yml",
)

const SAMPLE_YEARS = [2010, 2009, 2008, 2007, 2006, 2010]

sample_date_ranges =
    [(Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1)) for y in SAMPLE_YEARS]

const COVARIANCE_DATE_RANGES =
    [(Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1)) for y in 2006:2010]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_phase5_clt")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "pr", "swcre", "lwcre", "cl", "clt"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Floors from the validated runs: lwp/pr 0.2 (phase2), swcre 0.5 (irreducible
# brightness bias), lwcre 0.25 (phase2_lwcre025), cl 0.5 (cls zonal),
# clt 0.25 (clt-floor rerun, 10 degree arm).
const OBS_NOISE_GROUPS = [
    (short_names = ["lwp", "pr"], model_error_scale = 0.2),
    (short_names = ["swcre"], model_error_scale = 0.5),
    (short_names = ["lwcre"], model_error_scale = 0.25),
    (short_names = ["cl"], model_error_scale = 0.5),
    (short_names = ["clt"], model_error_scale = 0.25),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# Priors reuse the donor runs' (not their posteriors): the merged run should
# re-derive the answer from the data, not inherit confidence.
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
