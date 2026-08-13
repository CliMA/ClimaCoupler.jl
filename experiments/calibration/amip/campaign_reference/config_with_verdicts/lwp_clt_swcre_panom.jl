# Pattern-anomaly clt calibration: the zsw run with ONE change, clt
# graded on its ZONAL ANOMALY (the field minus its zonal mean) instead
# of the full field.
#
# Why. The EDMF run closed the scalar campaign: clt sits at its assumed
# floor and no parameter family reaches it under the full-field metric.
# Metric scoping (2026-08-08, existing iteration-1 ensembles vs CALIPSO
# Septembers, zero GPU) showed why and found the exception: most of the
# full-field clt misfit is zonal-mean amplitude the runs already fit;
# the PATTERN component alone has leverage ratio ~0.9-1.0 for the
# closure/optics family (full field ~1.3-2.3 in the same framework;
# Sc-deck box means are unreachable at ratio 1.5-3 and are rejected).
# The anomaly deletes the amplitude term that is already converged and
# spends clt's weight entirely on the component parameters can reach.
#
# Design:
# - ANOM_SHORT_NAMES = ["clt"]: zonal anomaly on the native grid, then
#   COARSEN_FACTOR 2, correlated 800 km floor, scale 0.25 unchanged.
#   The per-point floor (scale x seasonal mean) floors the anomaly
#   relative to local pattern strength; QuantileRegularization guards
#   the zero-crossing points.
# - clt zonal-mean amplitude is NOT in the loss. Guard: prediction 2
#   watches the clt global bias in the iteration bias plots.
# - lwp full-field 0.6 and swcre zonal 0.3 unchanged from zsw.
# - Same six closure/optics priors and rng_seed 42 as release/zsw: in
#   config space this is one change against zsw. NOTE the model is
#   ClimaAtmos 0.42.3 (post-rebase); sigma units are NOT comparable to
#   the 0.42.2 runs, so all predictions are within-run (iterations 1
#   and 6 both grade Sep 2006, giving same-weather endpoints).
# - Workers run in the preempt queue (0.2x charge) via
#   CALIBRATION_WORKER_QUEUE; checkpoints are off (900 days), so a
#   preempted member retries as a fresh start in a new output
#   increment. Preemption costs compute, not correctness.
#
# Predictions (falsifiable):
# 1. THE POINT. The clt pattern residual improves: its final value ends
#    below its iteration-1 value by more than the weather floor
#    (~0.13 sigma). This is the direct test of the scoping claim that
#    the pattern is reachable; if it fails, the reduced-space ratio ~1
#    was quadrature optimism and clt is closed at this rung for good.
# 2. The optimizer does not game the anomaly by wrecking amplitude:
#    clt global-mean bias at iteration 6 within 0.03 (clt fraction) of
#    iteration 1.
# 3. Zonal swcre stays within the weather floor of its start (the
#    cloud-removal trade must not reappear through the swcre block).
# 4. At least one parameter displaces >= 0.5 prior sigma: pattern
#    information updates the posterior, not just the residual.
#
# VERDICT (2026-08-08, 6 iterations, driver exited 0, preempt workers
# uneventful; iterations 1 and 6 both grade Sep 2006).
# clt pattern residual FLAT: 1.41 -> 1.40 sigma. lwp 0.68 -> 0.68,
# zonal swcre 0.90 -> 0.79. Ensemble COLLAPSED (0.0x prior spread).
# Frozen posterior (NOT adopted): eps_rel 0.078 (+2.0 sigma, the
# campaign's largest excursion), sharpness +0.45, abs_margin +0.41,
# Nd -0.4, residual -0.33, margin +0.2.
# P1 FAIL, decisively: with clt's entire weight on the pattern, the
#    residual did not move. The scoping ratio ~1 was quadrature
#    optimism: parameters move the pattern strongly but not TOWARD the
#    observed pattern (see CALIBRATION_LESSONS B9).
# P2 NOT JUDGED numerically; moot at zero spread, and eps_rel +2.0
#    implies the amplitude was spent chasing the pattern anyway.
# P3 PASS trivially (swcre at floor).
# P4 PASS on its metric (eps_rel +2.0) and meaningless: displacement at
#    zero spread is over-information, not information (LESSONS B10: the
#    anomaly's near-zero per-point floors over-inform the block).
# CONCLUSION: the last calibration route to clt at this rung is closed.
# Every family against the full field, and closure/optics against the
# isolated pattern: the clt pattern error is directionally orthogonal
# to the entire parameter space. Model development or the SCM rung.
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

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_panom")

const COARSEN_FACTOR = 2
const ZONAL_SHORT_NAMES = ["swcre"]
const ANOM_SHORT_NAMES = ["clt"]

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
    (short_names = ["swcre"], model_error_scale = 0.3),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

include(joinpath(@__DIR__, "..", "prior_tools.jl"))

const CALIBRATION_PRIORS = [
    checked_constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    checked_constrained_gaussian("cloud_fraction_floor_release_margin", 1, 0.5, 0, 3),
    checked_constrained_gaussian("cloud_fraction_floor_release_abs_margin", 1, 0.5, 0, 2),
    checked_constrained_gaussian("cloud_fraction_floor_release_sharpness", 1, 0.5, 0.2, 4),
    checked_constrained_gaussian("cloud_fraction_floor_residual", 0.05, 0.04, 0, 0.3),
    checked_constrained_gaussian(
        "prescribed_cloud_droplet_number_concentration",
        1e8, 8e7, 1e6, 1e9,
    ),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
