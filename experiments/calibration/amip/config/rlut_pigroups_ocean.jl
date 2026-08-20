# Pi-groups entrainment calibration against OCEAN-ONLY outgoing longwave
# radiation - the follow-up to rlut_pigroups (2026-08-19/20), which returned
# a clean null (posterior = prior, loss +10%, contraction 5.6x). Post-mortem
# findings this config acts on:
#
#   1. OCEAN ONLY. The rlut error is a convection-placement dipole whose
#      worst lobes sit over land (Amazon, Maritime Continent, +30..40 W/m^2).
#      ClimaAnalysis.apply_landmask (threshold 0.5) is applied as the LAST
#      preprocessing step before the sample collection (OCEAN_ONLY = true in
#      generate_observations.jl), dropping n from 2592 to ~1800. The sim side
#      needs no matching code: GEnsembleBuilder fills G with
#      ClimaAnalysis.flatten(sim_var, obs_metadata), so the observation's
#      drop_mask restricts the simulation automatically (ensemble_builder.jl
#      line ~257, verified). A mismatch fails at the observation-map stage,
#      which is cheaply retryable (member output already on disk).
#   2. model_error_scale 0.05 -> 0.02 (floor ~4.5 W/m^2, ~1x interannual).
#      The 0.05 floor equalled the measured TOTAL error and declared the
#      residual irreducible (whitened misfit 0.34 sigma -> nothing to chase).
#      The scale sweep (diff_model_error_scale/) puts the classically
#      well-specified point at 0.02 (whitened RMS ~0.93 sigma), which also
#      matches the irreducible-only error estimate from the 17-member
#      decomposition. EKP will now genuinely feel the ocean misfit; whether
#      the parameters can act on it is the question this run asks.
#   3. fixed_cloud_liquid_terminal_velocity reverted to the ClimaParams
#      default 0.01 (entry removed from the pigroups TOML; the zero-fclv
#      variant is preserved in git 238b0b10 and the archived run's dumps).
#   4. n_iterations = 4 (down from 6): cheap first look; extend by raising
#      the count if the trajectory warrants it.
#
# Go/no-go: run go_no_go.jl with GNG_MIN_SPREAD=2 (the 5 W/m^2 default was
# tuned to the old 11 W/m^2 floor).
#
# Inherited design (see rlut_pigroups.jl for full rationale):
#
# Loss: rlut (CERES EBAF Ed4.2.1, `toa_lw_all_mon`) alone, land AND ocean,
# 5 degree (144x72 -> 72x36), n = 2592 with ZERO missing points - rlut needs
# no coverage mask, unlike MAC lwp which drops 54.6% over land. An ocean-only
# variant is wanted later and has no mechanism yet: coverage masks are derived
# from the data's own NaNs, so it needs a land mask built from scratch.
#
# Target: October 2010, the SAME month every iteration. There is no
# minibatching over years here, deliberately. Grading identical weather every
# iteration makes the residual trajectory directly readable - the micro_edmf
# run cycled five Septembers and cost real interpretive effort disentangling
# weather from parameter response. The price is a posterior tuned to one
# October's weather, with no multi-year averaging to blunt it.
#
# MODEL: ClimaAtmos 0.42.6 from #main (tree ae62c083), pinned in
# experiments/AMIP/Manifest-v1.11.toml. Checked against the 0.42.4 the
# micro_edmf run used: src/prognostic_equations/edmfx_entr_detr.jl is
# BYTE-IDENTICAL between the two, so the closure below is exactly the code
# that was read when this config was designed. Commit the manifest with this
# config and do not update packages while the run or its analysis is live.
#
# PARAMETERS (3 -> 7 TransformUnscented members). Under
# `edmfx_entr_model: "PiGroups"` (set in the coupler YAML, which overrides
# climaatmos_progedmf_1m.yml's "Generalized"):
#
#   entr_vel_scale = limiter * max(0, c1|Pi_1| + ... + c5|Pi_5| + c6) / (z - z_sfc)
#
# with c = entr_param_vec[1:6]; elements 7-12 are never read. We sample
# elements 2, 3 and 6 through the `<base>_E<index>` convention added in
# commit dee61355: CalibrationTools splices each sampled scalar into the base
# vector from the run's coupler_toml BEFORE ClimaParams sees it, and errors on
# any name or index mismatch. Elements 1, 4 and 5 stay at the base value of 0.
#
#   Pi_2 = TKE / (w_j - w_0)^2, clamped to [-1, 1]  -> c2 weights turbulence
#   Pi_3 = sqrt(area fraction)                      -> c3 weights plume size
#   c6   = the constant term
#
# PRIOR CENTRE IS THE SHIPPED MODEL. With c2 = c3 = 0 and c6 = 0.3, pi_sum is
# identically 0.3, and entr_coeff in the production TOML is also 0.3 - so the
# centre reproduces the current InvZEntrainment closure exactly. That is a
# strong safety property the micro_edmf run did not have (its centre was an
# unstable configuration that NaN'd all 17 members).
#
# PRIOR WIDTH IS AN OPEN CONCERN, recorded rather than silently narrowed.
# c2, c3 ~ N(0, 0.5) on [-5, 5] as specified. Because both Pi_2 and Pi_3 are
# non-negative after abs() and bounded (Pi_2 <= 1, Pi_3 = sqrt(a) <= 0.84 at
# EDMF_max_area 0.7), a member with strongly negative c2 and c3 clamps pi_sum
# to 0 over much of the column - NO entrainment at all - while the positive
# tail reaches pi_sum ~ 9.5, about 32x the production 0.3. The micro_edmf run
# died at roughly 3x its tuned entrainment. At +-1 sigma the excursion is mild
# (pi_sum in roughly 0.3 +- 0.9), and TransformUnscented only ever evaluates
# sigma points, not tail draws - but the 7 sigma points must be audited and
# the centre smoke-tested before any ensemble launch. See PREDICTIONS below.
#
# SPINUP is 7 days, forced to start 2010-09-24. son_ic holds only Aug 25 and
# Sep 1 for 2006-2010, so this date must be GENERATED (get_initial_conditions.py
# then preprocessing.jl; the download stage needs CDS credentials including
# reanalysis-era5-complete access). start_date in the YAML is set to 20100924
# so a missing IC fails at preflight, not on a worker.
#
# OBSERVATION MODEL, all three settings measured rather than inherited:
#
#   model_error_scale = 0.05  (floor 11.2 W/m^2)
#     Measured from the finished micro_edmf run's own rlut diagnostic against
#     CERES: residual RMS 11.6-13.5 W/m^2 across members and iterations, so
#     the floor matches the mismatch a good parameter set actually leaves.
#     Reusing an lwcre-style 0.25 would assert 57 W/m^2 on a 225 W/m^2 field.
#     Caveat carried forward: D = (scale * mean field)^2 has the WRONG SHAPE -
#     |error| correlates only +0.30 with the mean field, and the floor comes
#     out ~2x too small in the ITCZ (error 22 vs floor 12) and 3-7x too big
#     poleward. At 0.05 that leaves tropical leverage ~1.8 and extratropical
#     ~0.4, which is the right emphasis for an ENTRAINMENT calibration even
#     though the shape is wrong for the general case.
#
#   COVARIANCE_DATE_RANGES = all 26 CERES Octobers
#     5 years gave SVD rank 4, and rank-4 captured "100%" of the variance only
#     because 4 was all the rank available; against 26 years those same four
#     directions account for 48.6%. Interannual std/pt rises 3.97 -> 4.61
#     W/m^2. Costs nothing (the SVD is of a 2592x26 matrix) and gives EKP real
#     noise in 25 directions instead of 4. The window carries a +3.35 W/m^2
#     forced trend which will land in a leading EOF; that inflates the
#     covariance, which is the conservative direction.
#
#   decorrelation_length = 800 km, identity_fraction at its 0.05 default
#     Tested, not assumed. The measured autocorrelation of the irreducible
#     residual field has a large nugget (C ~ 0.3 at the closest resolved
#     separation, so ~70% of the variance is uncorrelated at 5 degrees) and a
#     long tail on the remainder (fit L ~ 2300 km). Matching both gives
#     n_eff = 46 (Bayley-Hammersley) / 386 (participation ratio) against 75 /
#     211 for the configured kernel - the two standard measures disagree on
#     the SIGN of the correction, so 800 km sits inside the bracket the data
#     supports. Kept as-is. If it is ever revisited, identity_fraction is the
#     better-supported target than L.
#
# WHY THE GRID IS COARSENED. Not a statistical argument any more: n_eff is
# 75-211 at 5 degrees and 82-267 at 2.5, so quadrupling the points buys 10-25%
# more information for 16x the dense-covariance memory (54 MB -> 860 MB per
# observation). The correlated floor, not the grid, sets the information
# content. Coarsening is a compute decision and should be described as one.
#
# PREDICTIONS (record a verdict after the run - the micro_edmf config never
# got one, and that omission is why its lessons had to be re-derived):
# 1. Smoke test: a single centre member integrates 10 simulated days with no
#    NaN. The centre is the production closure, so failure here means
#    something is wrong with PiGroups wiring, not with the priors.
# 2. Go/no-go at iteration 1: rlut leverage ratio (residual / ensemble
#    response spread) <= ~3, AND reachable spread >= ~5 W/m^2. The micro_edmf
#    run's rlut reachable spread was 3.20 W/m^2 median against an 11.2 W/m^2
#    floor - S/N ~ 0.3 - which is the most likely explanation for its flat
#    convergence. If the Pi-groups parameters cannot beat that, six iterations
#    will not help and the run should be stopped.
# 3. Identifiability: at least one of c2, c3, c6 displaces >= 0.5 prior sigma
#    by the final iteration. All three flat means Pi-groups entrainment is
#    invisible to monthly 5-degree OLR and a different observable is needed.
# 4. No collapse: ensemble contraction stays under ~10x. The micro_edmf run
#    contracted 71x, which is what an over-informative covariance looks like.
#    Rank 25 instead of rank 4 should help; if it still collapses, the
#    decorrelation length is the next suspect.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_pigroups.yml",
)

# One fixed target, repeated. model_interface.jl indexes
# sample_date_ranges[iter] for iter in 1:n_iterations, so this vector must be
# at least n_iterations long; n_iterations + 1 matches the previous configs.
# Every entry is identical, so the minibatcher hands EKP the same observation
# every iteration and the residual trajectory is a same-weather comparison.
const _TARGET = (Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1))
sample_date_ranges = fill(_TARGET, 7)

# Every CERES October. SVDplusD requires each sample date to be one of these;
# 2010 is index 11.
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2000:2025
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_rlut_pigroups_ocean")

const COARSEN_FACTOR = 2

# Ocean-only loss: consumed by generate_observations.jl after coarsening.
const OCEAN_ONLY = true

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["rlut"],
    minibatch_size = 1,
    n_iterations = 4,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [(
    short_names = ["rlut"],
    model_error_scale = 0.02,
    decorrelation_length = OBS_DECORRELATION_LENGTH,
)]

# Unused by a 2-D-only loss, but preprocessing.jl calls
# select_pressure_levels / select_altitude_levels unconditionally.
const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

include(joinpath(@__DIR__, "..", "prior_tools.jl"))

# checked_constrained_gaussian, not the plain EKP constructor: EKP's silently
# returns a unit normal for small-magnitude targets, and c2/c3 are centred at
# exactly 0. Names follow the `<base>_E<index>` convention; run_calibration.jl
# calls CalibrationTools.check_element_priors on them before submitting
# anything, so a typo or an out-of-range index fails at launch rather than on
# a worker.
const CALIBRATION_PRIORS = [
    checked_constrained_gaussian("entr_param_vec_E2", 0, 0.5, -5, 5),
    checked_constrained_gaussian("entr_param_vec_E3", 0, 0.5, -5, 5),
    checked_constrained_gaussian("entr_param_vec_E6", 0.3, 0.1, 0, 1),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
