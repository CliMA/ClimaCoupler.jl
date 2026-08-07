# Zonal-swcre calibration: the release run with ONE change, swcre graded
# on its zonal mean instead of the full field.
#
# Why. The release run (lwp_clt_swcre_release.jl) ended with the same
# trade as every predecessor: swcre improved 20 percent, clt degraded 22
# percent. Its posterior leverage showed why. After Nd fixed the
# amplitude bias, the remaining swcre error is spatial pattern (ratio 2.9,
# a few percent attainable at best), and the only scalar-parameter route
# to reflection pattern is cloud amount, so the optimizer kept paying clt
# for unfixable swcre structure. Grading swcre zonally keeps the
# physically fixable amplitude constraint (the Nd pathway) and deletes
# the pattern term that finances cloud removal. lwp and clt keep the
# full 5 degree field. Restore full-field swcre when pattern-capable
# parameters (EDMF entrainment/detrainment) enter the set.
#
# Design, applying the so_jan lessons:
# - ZONAL_SHORT_NAMES = ["swcre"]: zonal mean on the native grid, 72
#   latitude values per sample. lwp and clt keep COARSEN_FACTOR 2.
# - swcre gets a DIAGONAL floor (no decorrelation_length): the 800 km
#   correlated floor double-counts under zonal reduction (so_jan null).
# - swcre model_error_scale 0.3, not 0.75: the reduction cuts swcre from
#   2592 constraints to 72, about 8x less pull after the measured 1.7-2.9x
#   per-point S/N gain, and (0.75/0.3)^2 = 6.25x restores the balance so
#   the amplitude bias still steers.
# - 72 values per variable, so the default regularization quantile 0.05
#   (needs 20) is fine.
#
# Everything else is byte-identical to the release run ON PURPOSE:
# same six priors, same rng_seed, same five Septembers, same fixed TOML
# (calibration_release_fixed.toml via the same yml). Same priors + seed
# give identical iteration-1 sigma points, so iteration 1 REUSES the
# release run's member outputs (copied member dirs + completed
# checkpoint markers; the WorkerBackend skips completed members and the
# driver recomputes G with this config's observation map). One of six
# GPU iterations saved, and the comparison to the release run is
# controlled: any outcome difference is the swcre representation.
#
# Predictions (falsifiable):
# 1. THE POINT. clt does not degrade: its final residual stays at or
#    below its iteration-1 value (0.76 sigma at the 0.25 scale). If clt
#    still degrades with the swcre pattern term gone, the cloud-removal
#    pressure comes from lwp or from clt's own noise floor, not from
#    swcre, and the pattern/amplitude explanation is wrong.
# 2. The zonal swcre amplitude still improves: final zonal swcre
#    residual at or below the release run's (its zonal-mean equivalent),
#    driven by Nd moving at least as far as the release's -0.94 sigma.
# 3. eps_rel stays below the release run's +1.17 sigma displacement,
#    since the swcre pattern chase was the remaining pressure on it.
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

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_zsw")

const COARSEN_FACTOR = 2
const ZONAL_SHORT_NAMES = ["swcre"]

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
    # Zonal: diagonal floor, rebalanced scale (see header).
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
