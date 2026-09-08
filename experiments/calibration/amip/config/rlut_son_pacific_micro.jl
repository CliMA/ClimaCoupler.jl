# rlut-ONLY ablation of rlut_rsut_son_pacific_micro: identical priors,
# simulation and domain; the loss drops the rsut block.
#
# WHY. The joint run (2026-08-27, 4/4 iterations) ended at a 1.7% improvement
# with reachability 5.8% and near-equal smooth response norms across all seven
# parameters. Because THIS config keeps the joint run's priors, seed and
# simulation unchanged, its UTKI sigma points are identical, so:
#   - iteration 1 is FREE: copy (or symlink) the joint run's
#     iteration_001/member_* into this run dir before launching; workers skip
#     on the `completed` checkpoints and the driver re-flattens the same
#     NetCDFs under the rlut-only metadata (mechanism proven on
#     rlut_son_pacific). All 15 members completed there - nothing to doctor.
#   - any difference between this run's parameter trajectory and the joint
#     run's isolates EXACTLY what the rsut block contributed to the update
#     (same iteration-1 information, different loss).
#
# NOISE MODEL: single rlut group, scale 0.02 (sweep-validated; floor 5.3
# W/m^2, whitened center residual ~0.9 sigma on this box).
#
# PREDICTIONS (grade after the run):
# 1. The parameter trajectory stays close to the joint run's (the joint run's
#    per-parameter response norms were nearly equal across blocks, so the
#    rsut block mostly repeated the rlut block's information).
# 2. If it DIVERGES, the divergence is in the microphysics/fall-speed
#    parameters (the SW-side levers), not E3/E6.
# 3. Contraction gentle; 15/15 members survive (identical corners to a run
#    that had zero deaths in four iterations).
#
# WALLTIME/LAUNCH: iteration 1 costs minutes; launch 1 with
# CALIBRATION_N_ITERATIONS=3 (recycle + iters 2-3 ~ 8.5 h, inside the 12 h
# worker walltime; the pacific run measured exactly this), launch 2 with 4.
# NO smoke test needed: the prior center is byte-identical to the joint
# run's, which was smoke-tested and then survived four full iterations.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_pigroups_son.yml",
)

const _TARGET = (Dates.DateTime(2010, 9, 1), Dates.DateTime(2010, 11, 1))
sample_date_ranges = fill(_TARGET, 7)

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1)) for y in 2000:2025
]

output_dir =
    joinpath(pkgdir(ClimaCoupler), "amip_calibration_rlut_son_pacific_micro")

const COARSEN_FACTOR = 2
const OCEAN_ONLY = true
const SEASONAL_MEAN = true
const REGION_LAT = (-30.0, 30.0)
const REGION_LON = [(120.0, 180.0), (-180.0, -90.0)]

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["rlut"],
    minibatch_size = 1,
    n_iterations = parse(Int, get(ENV, "CALIBRATION_N_ITERATIONS", "4")),
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Month(1),
    output_dir,
    rng_seed = 42,
)

const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [(
    short_names = ["rlut"],
    model_error_scale = 0.02,
    decorrelation_length = OBS_DECORRELATION_LENGTH,
)]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

include(joinpath(@__DIR__, "..", "prior_tools.jl"))

# IDENTICAL to rlut_rsut_son_pacific_micro.jl - this is what makes the sigma
# points, and therefore the recycled iteration 1, byte-identical. Do not edit
# one file without the other (or the reuse silently breaks).
const CALIBRATION_PRIORS = [
    checked_constrained_gaussian("entr_param_vec_E3", 0, 0.3, -5, 5),
    checked_constrained_gaussian("entr_param_vec_E6", 0.3, 0.1, 0, 1),
    checked_constrained_gaussian("rain_autoconversion_timescale", 1200, 300, 100, 7200),
    checked_constrained_gaussian(
        "cloud_ice_specific_humidity_autoconversion_threshold",
        5e-6,
        1e-6,
        1e-6,
        1e-5,
    ),
    checked_constrained_gaussian("fixed_snow_terminal_velocity", 1, 0.3, 0.1, 3),
    checked_constrained_gaussian("fixed_cloud_ice_terminal_velocity", 0.01, 0.003, 0, 0.2),
    checked_constrained_gaussian("detr_massflux_vertdiv_coeff", 0.3, 0.1, 0.1, 5.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
