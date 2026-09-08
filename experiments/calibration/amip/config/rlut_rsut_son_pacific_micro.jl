# Joint rlut + rsut calibration over the TROPICAL PACIFIC box with
# microphysics + detrainment parameters - AMIP mode, SON 2010 seasonal mean,
# ocean points. Successor to rlut_son_pacific.
#
# WHY THIS RUN. The entrainment-only campaign is closed: rlut_son_pacific
# (2026-08-26) confirmed the pre-registered null - the three Pi-group
# entrainment parameters see ~10% of the box residual in the whitened metric,
# reproduce the global run's answer (c6 -> 0.26 under both losses), and cannot
# touch the placement dipole. This run changes the LEVERS, not the domain:
#   - 1M microphysics parameters change cloud optical properties WHERE cloud
#     exists (condensate-weighted fingerprints, not broad tint);
#   - detr_massflux_vertdiv_coeff reshapes where mass flux detrains
#     vertically (anvil height/amount);
#   - rsut joins the loss so the update cannot fix OLR by thinning anvils at
#     albedo's expense - the LW/SW compensation we previously only caught
#     post-hoc is now penalized inside the update.
#
# PARAMETERS (7 -> UTKI 2p+1 = 15 members):
#   entr_param_vec_E3, _E6 as before; E2 is DROPPED and pinned at the toml
#   default 0 - it was degenerate with E3 (learned anticorrelation valley)
#   and its +sqrt(3)sigma corner is the only proven NaN crash (3 runs).
#   Five plain (non-spliced) parameters, all verified present in the
#   ClimaParams registry and on this config's active code path
#   (microphysics_model "1M" + prognostic EDMF):
#     rain_autoconversion_timescale              default 1000   prior N(1200, 300)  on [100, 7200]
#     cloud_ice_..._autoconversion_threshold     default 1e-6   prior N(5e-6, 1e-6) on [1e-6, 1e-5]
#     fixed_snow_terminal_velocity               default 1.0    prior N(1, 0.3)     on [0.1, 3]
#     fixed_cloud_ice_terminal_velocity          default 0.01   prior N(0.01, 0.003) on [0, 0.2]
#     detr_massflux_vertdiv_coeff                default 1      prior N(0.3, 0.1)   on [0.1, 5]
#   NOTE two prior means sit far from the model defaults (ice autoconversion
#   threshold 5x default; detr coeff 0.3 vs default 1, whole prior below the
#   default). The iteration-1 center is therefore a DIFFERENT model from every
#   previous run's center - deliberate, per the run design.
#   All priors use checked_constrained_gaussian: the 5e-6 and 0.01 targets are
#   in the regime where EKP's constrained_gaussian silently returns a unit
#   normal (the July 2026 q_liq incident).
#
# NOISE MODEL. Two independent SVDplusD blocks (rlut, rsut), each with the
# 800 km correlated floor; cross-covariance between the blocks (ENSO moves
# both) is NOT modeled - known simplification, same as the lwp+swcre runs.
# Scales are PER-GROUP because floor = scale x field-mean does not transfer
# across variables (rlut box mean ~263 W/m^2, rsut ~100):
#   rlut: 0.02, measured optimal by the 2026-08-25 covariance sweep
#         (floor/interannual 0.64, whitened iteration-1 residual 0.91 sigma).
#   rsut: set by measure_rsut_scale.jl before prep (targets floor/interannual
#         0.6-1.0 and ~1 sigma whitened center-member residual; the previous
#         runs' members already output monthly rsut, so the misfit anchor is
#         free). Reusing 0.02 here would put the floor near ~2 W/m^2 against
#         comparable interannual spread - the lwcre failure mode.
#   A wrong scale is recoverable WITHOUT re-running members: sigma points
#   depend only on the priors, so regenerate observations with the corrected
#   scale and recycle iteration_001 via the checkpoint-copy procedure proven
#   on rlut_son_pacific. Pre-registered decision rule: whitened iteration-1
#   residual outside [0.5, 2] sigma per block -> stop, rescale, recycle.
#
# PRE-REGISTERED PREDICTIONS (grade after iteration 1 / after the run):
# 1. Whitened reachable variance of the joint box residual across the 14
#    sigma-point response directions exceeds the entrainment-only 10.4%.
#    If it does not, no parameter set of this family reaches the dipole.
# 2. The rsut response is dominated by the three fall-speed/autoconversion
#    parameters; the entrainment pair contributes mostly through rlut.
# 3. Contraction stays gentle (< ~2x per axis per iteration) at the measured
#    scales; no lwcre-style collapse.
# 4. No LW/SW seesaw: the center-member rlut and rsut residuals do not move
#    in opposite directions across iterations.
# 5. All 15 members survive iteration 1 (the proven crash corner is gone;
#    moderate confidence - twelve corners are untested).
#
# WALLTIME/LAUNCH. Simulation length unchanged (Aug-Nov 2010) -> ~4.3 h per
# iteration with all 15 workers up (4 nodes at 4 workers/node). Same
# clean-boundary policy as before: 2 iterations per 12 h launch via
# CALIBRATION_N_ITERATIONS (2, then 4). No recycled first iteration - the
# priors changed, so previous sigma-point outputs do not apply.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_pigroups_son.yml",
)

# One fixed target, repeated (see rlut_son_pacific.jl).
const _TARGET = (Dates.DateTime(2010, 9, 1), Dates.DateTime(2010, 11, 1))
sample_date_ranges = fill(_TARGET, 7)

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1)) for y in 2000:2025
]

output_dir =
    joinpath(pkgdir(ClimaCoupler), "amip_calibration_rlut_rsut_son_pacific_micro")

const COARSEN_FACTOR = 2

# Ocean-only + tropical Pacific box, applied to BOTH variables (the
# preprocessing maps over vars).
const OCEAN_ONLY = true
const SEASONAL_MEAN = true
const REGION_LAT = (-30.0, 30.0)
const REGION_LON = [(120.0, 180.0), (-180.0, -90.0)]

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["rlut", "rsut"],
    minibatch_size = 1,
    n_iterations = parse(Int, get(ENV, "CALIBRATION_N_ITERATIONS", "4")),
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Month(1),
    output_dir,
    rng_seed = 42,
)

const OBS_DECORRELATION_LENGTH = 8.0e5

# Per-variable noise groups: separate SVDplusD blocks and separate floors.
# rlut scale is the sweep-validated 0.02. The rsut scale 0.10 was MEASURED
# 2026-08-26 (measure_rsut_scale.jl): box mean 86.6 W/m^2, interannual
# sigma ~8.4, model misfit rms 17.1-17.5 (two member anchors) -> at 0.10 the
# floor/interannual ratio is 0.82 and the whitened center residual 0.90
# sigma/DOF, matching the healthiest known operating point (rlut at 0.91).
# 0.02 here would give ratio 0.21 and 4.2 sigma - the lwcre failure mode.
const RSUT_MODEL_ERROR_SCALE =
    parse(Float64, get(ENV, "CALIBRATION_RSUT_SCALE", "0.10"))
const OBS_NOISE_GROUPS = [
    (
        short_names = ["rlut"],
        model_error_scale = 0.02,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["rsut"],
        model_error_scale = RSUT_MODEL_ERROR_SCALE,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
]

# Unused by a 2-D-only loss, but preprocessing.jl calls
# select_pressure_levels / select_altitude_levels unconditionally.
const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

include(joinpath(@__DIR__, "..", "prior_tools.jl"))

# checked_constrained_gaussian for every prior (not just the _E# ones): the
# 5e-6 and 0.01 targets are exactly the small-magnitude regime where EKP's
# plain constructor silently ignores the requested moments.
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
