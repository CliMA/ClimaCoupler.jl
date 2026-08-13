# Post-rebase reproducibility check of the floor-release calibration.
#
# Rerun of lwp_clt_swcre_release.jl at HALF the iterations (3 of 6), after
# rebasing this branch onto main (2026-08-12): registered ClimaCalibrate
# 0.4.0 replaces the ne/worker-packing 0.3.2 branch, generate_observations
# was ported to the SampleBuilder API, and the environment moved to main's
# Manifest. Everything scientific is byte-identical to the original config:
# same coupler YAML, dates, grid, floors, priors, and seed. Only
# n_iterations (3) and output_dir differ.
#
# Purpose: confirm the migrated pipeline reproduces the original run's
# trajectory before trusting it for new experiments.
#
# Environment deltas vs the original run (2026-08-06/07, tmux driver):
#   ClimaCalibrate 0.3.2 ne/worker-packing -> 0.4.0 registered
#   ClimaAtmos     0.42.3                  -> 0.42.4
#   ClimaAnalysis  0.5.22                  -> 0.5.23
#
# Predictions (falsifiable):
# 1. The observation vector is identical to the original run's (checked
#    against /glade/derecho/scratch/nefrathe/amip_calibration_release_out
#    before launch; the pre/post-rebase A/B already matched to ~1e-13).
# 2. Iteration-1 member parameters are IDENTICAL to the original run's
#    (same priors + seed 42 + TransformUnscented sigma points are
#    deterministic; check iteration_001/member_*/parameters.toml).
# 3. The 3-iteration trajectory tracks the original within the weather
#    floor: same displacement signs, residuals within ~0.1 sigma of the
#    original per iteration. Larger drift implicates the ClimaAtmos
#    0.42.3 -> 0.42.4 bump (the campaign once measured 27% ocean-lwp
#    drift from an Atmos upgrade at identical parameters; record, do not
#    average over it).
#
# VERDICT (run 2026-08-12 18:13 -> 22:15 MDT, tmux driver, PBS workers
# 7104267/70/71/72, 3 iterations, driver exited 0; run dir
# /glade/derecho/scratch/kphan/amip_calibration_release_repro_out,
# comparison via compare_trajectory.jl in that dir):
# 1. PASS. Observation vector identical to the original run's (all
#    samples, covariance blocks, masks; diffs <= ~1e-13 relative).
# 2. PASS. Iteration-1 parameters byte-identical, 13/13 members.
# 3. FAIL as stated, with the cause identified and NOT the machinery.
#    Loss ran 6-18% below the original (0.0410/0.0335/0.0298 vs
#    0.0437/0.0408/0.0360); five of six parameters tracked within
#    0.03-0.3 sigma with matching signs, but abs_margin FLIPPED
#    (repro ~0.0 sigma vs original +0.60): the eps_rel-decoupling story
#    of the original run does not hold under the new physics.
#    ISOLATION: at byte-identical iteration-1 parameters, G changed
#    23% relative RMS overall (lwp 41%, clt 23%, swcre 23%), i.e.
#    0.66-0.91x the ensemble's own parameter-driven spread. Cause: the
#    rebase onto main brought cloud_ice_formation: TemperatureDependent
#    into climaatmos_progedmf_1m.yml (commit 674ccf36, 2026-08-08) plus
#    ClimaAtmos 0.42.3 -> 0.42.4. The forward model is a DIFFERENT
#    MODEL; this is the guide's environment-drift pitfall measured
#    again (previously: 27% ocean lwp from an Atmos upgrade).
# CONCLUSION: the migrated calibration machinery reproduces exactly
# (deterministic chain bit-identical end to end); the campaign's
# scientific posteriors do NOT transfer to the new physics unchanged.
# Treat campaign posteriors as priors to re-verify on this model, per
# the ladder's transfer rules, and record the physics flag with any
# future comparison to campaign-era runs.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_release.yml",
)

const _SEP_YEARS = [2006, 2010, 2008, 2007, 2009]
sample_date_ranges = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 9, 1)) for
    y in _SEP_YEARS[1 .+ ((0:6) .% length(_SEP_YEARS))]
]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 9, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_release_repro")

const COARSEN_FACTOR = 2

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "clt", "swcre"],
    minibatch_size = 1,
    # Half of the original 6: enough to compare the early trajectory
    # (iterations 1-3 grade Sep 2006, 2010, 2008, as in the original).
    n_iterations = 3,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# clt 0.35 -> 0.25: with leverage established (sep_multi ratio 2.2), clt's
# remaining problem is priority against swcre, and this is the veto.
const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [
    (
        short_names = ["lwp"],
        model_error_scale = 0.6,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["clt"],
        model_error_scale = 0.25,
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

include(joinpath(@__DIR__, "..", "prior_tools.jl"))

const CALIBRATION_PRIORS = [
    checked_constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    # Release shape (see _compute_cloud_fraction): margin scales the
    # equilibrium-width term of the release width, abs_margin the
    # eps_rel*q_sat term (freeing it breaks the eps_rel coupling),
    # sharpness the transition exponent. Defaults are all 1.
    checked_constrained_gaussian("cloud_fraction_floor_release_margin", 1, 0.5, 0, 3),
    checked_constrained_gaussian("cloud_fraction_floor_release_abs_margin", 1, 0.5, 0, 2),
    checked_constrained_gaussian("cloud_fraction_floor_release_sharpness", 1, 0.5, 0.2, 4),
    # residual caps CF below 1 inside fully-released saturated decks, so
    # the prior keeps its mass at small values and the bound stays well
    # under where it visibly caps stratocumulus. sigma 0.04 not 0.05:
    # EKP rejects mu - sigma at the lower bound.
    checked_constrained_gaussian("cloud_fraction_floor_residual", 0.05, 0.04, 0, 0.3),
    # The optical knob: N0 multiplying the aerosol-modulated droplet number
    # behind the Liu-Hallett effective radius. Prior from row3_optics.
    checked_constrained_gaussian(
        "prescribed_cloud_droplet_number_concentration",
        1e8, 8e7, 1e6, 1e9,
    ),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
