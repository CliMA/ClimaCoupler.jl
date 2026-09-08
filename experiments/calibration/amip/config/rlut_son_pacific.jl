# Pi-groups entrainment calibration against rlut over the TROPICAL PACIFIC
# only - AMIP mode, SON 2010 seasonal mean, ocean points.
#
# WHY A REGION. Every previous run graded a near-global field in which the
# reachable fraction of the residual was single-digit: the misfit is a
# convection-PLACEMENT dipole, and entrainment coefficients rescale plume
# dilution rather than relocating convection. Restricting the loss to the box
# where BOTH lobes of that dipole live - the Maritime Continent / warm pool
# through the East Pacific ITCZ - removes extratropical and Atlantic/Indian
# points that contribute noise and irreducible bias but no reachable signal.
# The parameters act on tropical deep convection; the loss should be measured
# there.
#
#   REGION_LAT = (-30, 30)
#   REGION_LON = [(120, 180), (-180, -90)]   two boxes: the -180..180 grid
#                                            splits the Pacific at the dateline
#
# Applied by apply_region_mask in preprocessing.jl after the ocean mask and
# before the seasonal mean. It NaNs everything outside the box; the
# simulation side needs no matching code, since the observation metadata's
# drop_mask carries the mask into GEnsembleBuilder's flatten (same mechanism
# the ocean mask uses, verified).
#
# EXPECT n TO DROP HARD: roughly 12 of 36 latitudes and 30 of 72 longitudes
# survive, so ~360 grid points before the ocean mask and perhaps ~300 after,
# versus 1696 for rlut_son_ocean. Consequences to check at prep time rather
# than assume:
#   - The rank-25 interannual covariance is now a much larger fraction of a
#     ~300-dimensional space, so the floor governs proportionally less and
#     conditioning improves.
#   - model_error_scale is INHERITED AT 0.02 but was calibrated on the
#     near-global field. The tropics have both a larger model error (this is
#     where the dipole is) and a different interannual spread, so the
#     floor/interannual ratio will move. Read it off the prep log and adjust
#     if it lands outside ~0.8-2.2, the band used for the previous runs.
#   - The 800 km correlated floor is unchanged, but with a smaller domain the
#     effective number of independent constraints falls too; watch the
#     iteration-1 contraction.
#
# STILL OCEAN-ONLY (inherited). Within this box that removes the Maritime
# Continent islands and the Central/South American coast. If the intent is to
# grade tropical Pacific convection wherever it occurs, land included, drop
# OCEAN_ONLY - it is one line.
#
# Everything else is inherited from rlut_son_ocean: AMIP mode, SON 2010
# seasonal mean (1-month spinup from 2010-08-01, Aug discarded by the season
# window), covariance over SON 2000-2025, c2/c3 priors at sigma 0.3, scale
# 0.02, 800 km kernel, 4 iterations, seed 42.
#
# WALLTIME: unchanged at ~4.3 h per iteration (the simulation is identical;
# only the loss changes), so 2 iterations per launch via
# CALIBRATION_N_ITERATIONS, or the relay.
#
# PREDICTIONS (grade after the run):
# 1. Reachable fraction of the residual is substantially higher than the
#    ~4-7% measured on the near-global fields - that is the entire premise.
#    If it is not, the problem is the parameters, not the domain.
# 2. All 7 members survive iteration 1 (unchanged priors and simulation).
# 3. c6 moves DOWN as in rlut_son_ocean (that run had rlut biased +2.65 W/m^2
#    and wanted less entrainment). A tropics-only loss should want this more
#    strongly, since the box is where the positive bias is concentrated.
# 4. lwcre continues to improve alongside rlut (no compensating error), as it
#    did in rlut_son_ocean. A tropics-only loss makes the cloud response the
#    dominant term, so a reversal here would be informative.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_pigroups_son.yml",
)

# One fixed target, repeated. model_interface.jl indexes
# sample_date_ranges[iter] for iter in 1:n_iterations, so this vector must be
# at least n_iterations long; n_iterations + 1 matches the previous configs.
# Every entry is identical, so the minibatcher hands EKP the same observation
# every iteration and the residual trajectory is a same-weather comparison.
const _TARGET = (Dates.DateTime(2010, 9, 1), Dates.DateTime(2010, 11, 1))
sample_date_ranges = fill(_TARGET, 7)

# Every CERES October. SVDplusD requires each sample date to be one of these;
# 2010 is index 11.
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1)) for y in 2000:2025
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_rlut_son_pacific")

const COARSEN_FACTOR = 2

# Ocean-only loss: consumed by generate_observations.jl after coarsening.
const OCEAN_ONLY = true

# Grade seasonal time means over the sample/covariance windows (see header).
const SEASONAL_MEAN = true

# Tropical Pacific box. Two longitude ranges because the -180..180 grid the
# regridder produces splits the Pacific at the dateline; apply_region_mask
# takes their union. Everything outside is NaN'd.
const REGION_LAT = (-30.0, 30.0)
const REGION_LON = [(120.0, 180.0), (-180.0, -90.0)]

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["rlut"],
    minibatch_size = 1,
    # Overridable per launch so the run can be split at CLEAN ITERATION
    # BOUNDARIES instead of letting worker walltime kill members mid-flight:
    # a rerun member restarts from its 10-day checkpoint
    # (detect_restart_files true) but ClimaDiagnostics' monthly accumulation
    # state is NOT checkpointed, so its monthly means - the G entries - would
    # be silently computed from partial months. Launch 1 with
    # CALIBRATION_N_ITERATIONS=2 (exits cleanly ~9 h, inside the 12 h worker
    # walltime); launch 2 with 4 (resume skips iterations 1-2). Never point
    # the relay at this config with a target beyond what one walltime holds.
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
    checked_constrained_gaussian("entr_param_vec_E2", 0, 0.3, -5, 5),
    checked_constrained_gaussian("entr_param_vec_E3", 0, 0.3, -5, 5),
    checked_constrained_gaussian("entr_param_vec_E6", 0.3, 0.1, 0, 1),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
