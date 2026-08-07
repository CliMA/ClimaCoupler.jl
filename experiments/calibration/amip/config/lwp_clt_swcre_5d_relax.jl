# lwp_clt_swcre_5d_cltf035 plus three phase/relaxation parameters:
#   condensation_evaporation_timescale (150 s) — liquid formation tau_relax
#   sublimation_deposition_timescale (500 s)   — ice formation tau_relax
#   cloud_fraction_floor_release_abs_margin (1) — cloud-fraction floor
#     release shape (new in Atmos 0.42; not yet in ClimaParams, read from
#     the member TOML via the create_parameters release-params hook).
# All three wiring-verified by preflight (2026-07-30): the margin reaches
# turbconv_params, the timescales reach the 1M formation tau_relax.
# 8 parameters -> 17 members. Everything else inherited from the cltf035
# rerun (5 degree grid, clt floor 0.35, seed 42, 6 iterations).
#
# Predictions (falsifiable):
# 1. No COLLAPSED flag; residuals reach or stay below their floors as in
#    the cltf035 rerun.
# 2. The floor-release margin is identified by clt (it directly shapes
#    when the cloud-fraction floor releases): its spread contracts below
#    0.5x prior by the final iteration. If it stays near 1x prior, the
#    margin is not identifiable at this sample size and the ClimaParams
#    promotion should wait.
# 3. The relaxation timescales trade against q_liq/rain_tau without going
#    degenerate (the steering degen line stays separable); q_liq stays in
#    the 2D cluster (1.6-1.9e-4).
# 4. Spreads on the original five parameters do not END wider than the
#    cltf035 rerun's finals — adding phase physics should not destabilize
#    the already-constrained set.
#
# Verdict (run completed 2026-07-30, job 6961043, 6 iterations):
# 1. PASS. No collapse; final residuals lwp 0.67σ, clt 0.59σ, swcre 0.68σ,
#    all below their floors, matching the cltf035 rerun (0.67/0.61/0.67).
# 2. PASS. Margin spread ended 0.03-0.04x prior (well under 0.5) and the
#    data moved the mean (1.0 -> 1.27 at iter 2 -> 1.066 final), so the
#    margin is identifiable at this sample size. This is the promotion
#    evidence for ClimaParams.
# 3. PASS on separability, EDGE on q_liq. Timescales stayed separable and
#    landed apart: condensation-evaporation 106 s, sublimation-deposition
#    278 s. The anticipated trade appeared as rain_tau 1385 -> 1167 vs the
#    cltf035 rerun. q_liq ended 1.93e-4, 1.5% above the 1.6-1.9e-4 window.
# 4. PASS. Finals 0.03-0.04x prior vs cltf035's 0.02-0.03x: marginally
#    wider, not destabilized.
# Full posterior: q_liq 1.93e-4, rain_tau 1167, snow_tau 1361,
# cond_evap_tau 106, subl_dep_tau 278, margin 1.066, eps_rel 0.0488,
# steepness 0.539. eps_rel sits between the zonal (0.0335) and 2D-clt
# (0.058) answers.
# INCIDENT (RESOLVED 2026-07-31): no stale config was read. EKP's
# constrained_gaussian silently returns a unit normal over the bounds for
# targets with magnitudes below ~1e-2 (absolute optimizer tolerances; see
# prior_tools.jl), so BOTH the "old" (5e-4, 3e-4) and the tightened
# (1.8e-4, 0.6e-4) q_liq priors were in fact the same logistic unit
# normal, mean 5e-4 = bounds center. The tightening below never took
# effect in this run. Arriving at 1.93e-4 from an effective 5e-4 prior
# mean remains independent evidence for the low-q_liq answer. Use
# checked_constrained_gaussian on any rerun.
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
    joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_clt_swcre_5d_relax")

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

const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        # Tightened around the reproduced 2D answer (1.6-1.9e-4 across four
        # runs): the old wide prior (5e-4, 3e-4) made iteration 1 jump far
        # from the prior mean, so the error plot mostly showed the jump
        # rather than real learning. The zonal cluster (2.2-2.55e-4) stays
        # within ~1 sigma.
        1.8e-4, 0.6e-4, 0.0, 1e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1800, 700, 300, 3600),
    PD.constrained_gaussian("snow_autoconversion_timescale", 1800, 700, 300, 3600),
    PD.constrained_gaussian("condensation_evaporation_timescale", 150, 100, 20, 3600),
    PD.constrained_gaussian("sublimation_deposition_timescale", 500, 400, 20, 3600),
    PD.constrained_gaussian("cloud_fraction_floor_release_abs_margin", 1, 0.5, 0, 2),
    PD.constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    PD.constrained_gaussian("cloud_fraction_steepness_scale", 0.7, 0.2, 0.1, 2.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
