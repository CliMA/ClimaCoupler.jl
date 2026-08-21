# Pi-groups entrainment calibration against OCEAN-ONLY LONGWAVE CLOUD
# RADIATIVE EFFECT (lwcre = rlutcs - rlut, CERES-derived) - the controlled
# falsification companion to rlut_pigroups_ocean. EVERYTHING except the
# observable is deliberately identical to that run (same priors, seed, mask,
# scale 0.02, 4 iterations, YAML, TOML), so any difference in the posterior
# is attributable to the loss alone.
#
# WHY (measured 2026-08-20 from the rlut_pigroups_ocean output, zero GPU):
# the rlut run's c6 rise 0.299 -> 0.341 removed the rlut bias (-0.99 -> -0.07
# W/m^2) while degrading lwcre MONOTONICALLY (bias -4.70 -> -5.46, RMS 9.26
# -> 9.78): a compensating-error fit. Delta rlutcs ~ -6.5 W/m^2 means the
# OLR deficit is CLEAR-SKY (humidity/temperature), and the calibration paid
# for it by thinning already-deficient high cloud. lwcre isolates the cloud
# part and should therefore pull the opposite way.
#
# PREDICTIONS (grade after the run):
# 1. Members 2 and 3 (c2/c3 = +0.87 sigma corners) NaN in iteration 1 at the
#    same simulated days - iteration 1's simulations are identical to the
#    rlut_pigroups_ocean run's (the loss does not touch the forward model).
# 2. Over-confident noise model: the scale-0.02 floor is 0.52 W/m^2 =
#    0.11x the 4.0 W/m^2 interannual spread (the honest lwcre scale is ~0.2;
#    the campaign used 0.25). Whitened iteration-1 residual ~2-2.5 sigma;
#    contraction well beyond the rlut run's 24.8x.
# 3. THE HEADLINE: c6 is pulled DOWN, below 0.30 (rlut run: UP to 0.34).
#    Conflicting single-observable posteriors = neither is cloud physics;
#    the productive successor is a JOINT rlut+lwcre (+rlutcs) loss.
# 4. Apparent success: the -5.5 W/m^2 mean offset is 31% of lwcre's MSE and
#    reachable, so lwcre bias shrinks dramatically and RMS improves - while
#    rlut degrades in mirror image (watch it free in the saved diagnostics).
#
# lwcre facts (ocean-masked, Octobers 2000-2025): mean 25.0 W/m^2 (5..72),
# interannual 4.0 W/m^2 (19% of local mean; rlut: 2%), model RMS 9.8 = 39%
# of mean. Sim side computes lwcre = rlutcs - rlut in observation_map.jl
# from diagnostics every member already writes.
#
# Inherited design: see rlut_pigroups_ocean.jl and rlut_pigroups.jl.
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

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwcre_pigroups_ocean")

const COARSEN_FACTOR = 2

# Ocean-only loss: consumed by generate_observations.jl after coarsening.
const OCEAN_ONLY = true

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwcre"],
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
    short_names = ["lwcre"],
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
