# Pi-groups entrainment calibration against OCEAN-ONLY SON-SEASON rlut in
# AMIP MODE - the successor run applying every lesson of the campaign:
#
#   1. AMIP mode (evolving observed SST/sea ice) replaces subseasonal
#      persistence forcing, making a multi-month graded window legitimate.
#      No WeatherModel/era5_initial_condition_dir - AMIP mode ignores them.
#   2. Target: the SON 2010 SEASONAL MEAN (time average of the Sep, Oct and
#      Nov monthly means; the Aug spinup month is discarded by the window).
#      SEASONAL_MEAN = true applies ClimaAnalysis.average_season_across_time
#      on BOTH sides (no local averaging code): each season becomes one
#      time-mean slice, the time dimension is KEPT with SON stamped at its
#      first date (Sep 1), and date selection picks the SON slices ->
#      n = 1696 ocean constraints. Verified on CERES: the 2010 SON slice is
#      bit-equal to the hand mean of the Sep+Oct+Nov monthly maps.
#      Averaging 3 months cuts weather noise ~sqrt(3).
#   3. Spinup 1 month (Aug 1 2010 start); atmosphere equilibrates from the
#      generic AMIP IC well within it. Land does NOT (see YAML note) -
#      accepted risk under the ocean-only loss.
#   4. c2/c3 prior sigma 0.5 -> 0.3: sigma points at +-0.52, inside the
#      measured stability envelope (+-0.87 NaN'd members 2/3 in BOTH ocean
#      runs, and the FIRST update is where everything gets decided - it
#      must see all 7 members).
#   5. Everything else inherited from rlut_pigroups_ocean: ocean mask,
#      scale 0.02 (floor ~4.5 W/m^2; note the SON-mean interannual spread
#      is smaller than single-month, so the effective floor/interannual
#      ratio rises toward ~1.5-2 - conservative side), 800 km kernel,
#      4 iterations (T=0.4), seed 42.
#
# WALLTIME: 122-day members ~ 4-4.5 h -> 4 iterations ~ 17 h > the 12 h
# worker walltime. LAUNCH THROUGH THE RELAY (2 iterations per launch):
#   tmux new-session -d -s cal_rlut_son \
#     'bash experiments/calibration/amip/calibration_relay.sh <run_dir> 4 \
#        bash <run_dir>/driver_rlut_son_ocean.sh'
#
# MODEL: ClimaAtmos 0.42.7 #main (tree 1f031c77, depot slug p0ScH) - the
# working-tree manifest pin at staging time; edmfx_entr_detr.jl and the
# entr_param_vec mapping verified BYTE-IDENTICAL to the 0.42.4/0.42.6
# trees all previous runs used.
#
# PREDICTIONS (grade after):
# 1. All 7 members survive iteration 1 (sigma points inside the envelope).
# 2. Whitened iteration-1 residual ~1 sigma (scale 0.02 vs SON-mean-ish
#    spread); contraction < 30x total. If >> that, the floor story repeats.
# 3. The genuine open question after the lwcre falsification: with honest
#    conditioning and all members alive, does c6 still move up? A c6 rise
#    that also IMPROVES ocean rlut bias without degrading lwcre (watch the
#    free rlutcs/rlut diagnostics) would be the first defensible parameter
#    result; c6 flat reproduces the null cleanly.
# 4. Weather-noise reduction: iteration-1 physical member spread per point
#    smaller than the single-October runs' (~4-5 W/m^2).
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

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_rlut_son_ocean")

const COARSEN_FACTOR = 2

# Ocean-only loss: consumed by generate_observations.jl after coarsening.
const OCEAN_ONLY = true

# Grade seasonal time means over the sample/covariance windows (see header).
const SEASONAL_MEAN = true

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["rlut"],
    minibatch_size = 1,
    n_iterations = 4,
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
