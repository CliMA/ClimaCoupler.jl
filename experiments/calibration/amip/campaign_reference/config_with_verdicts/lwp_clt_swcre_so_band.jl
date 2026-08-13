# Southern Ocean band calibration: 45S to 65S only, lwp + clt + swcre,
# October 2010, subseasonal with a 7-day spinup.
#
# Why a band. The global vterm run showed these parameters barely move clt
# (residual 3.7x the ensemble spread, and 6.5x after removing the weather
# floor), and that lwp and swcre gains came mostly from the tropics and
# subtropics. The Southern Ocean is where fall speeds, cloud cover and
# shortwave cloud radiative effect should couple, so a band restricted to
# it removes regions that dilute the signal. Everything else follows the
# vterm run so the two can be compared.
#
# Parameters: 5 free, 11 members. The autoconversion threshold is FIXED at
# 2e-4 in toml/calibration_so_fixed.toml, and the cloud-fraction shape and
# phase timescales stay at the relax posterior in
# toml/calibration_relax_posterior.toml.
#
# Two deliberate departures from the global runs:
#
# 1. COARSEN_FACTOR = 1, the native 2.5 degree grid. The band holds only
#    about 11 percent of the globe, so coarsening to 5 degrees would leave
#    roughly 290 clt points. Coarsening raises the valid FRACTION (a block
#    survives if any sub-cell does) but always lowers the point COUNT: in
#    60-55S, lwp has about 180 valid cells at 2.5 degrees against 54 at 5
#    degrees. Count is what matters against the internal-variability floor.
#
# 2. OBS_DECORRELATION_LENGTH = 4.0e5, down from 8.0e5. The band is only
#    about 2200 km wide, so an 800 km correlated floor would treat a
#    band-wide cloud bias as expected error. That bias is the whole point
#    of this run. The clt diagnostic showed the correlated floor absorbing
#    exactly this kind of smooth error: clt's residual was 0.62 sigma on
#    the covariance diagonal but 0.33 sigma whitened.
#
# Data note. MAC lwp is passive microwave over ocean and October is near
# the Antarctic sea-ice maximum, so lwp coverage in October 2010 runs 96
# percent at 50-45S, 94 at 55-50S, 63 at 60-55S, 14 at 65-60S and 0 south
# of 65S. clt (CALIPSO) and swcre (CERES) cover the band fully. lwp
# therefore constrains only the northern three quarters of the band. This
# is a limit of the instrument, not of the processing, and no resolution
# choice recovers it.
#
# Predictions (falsifiable):
# 1. THE TEST. clt's leverage ratio (residual / ensemble spread at
#    iteration 1, from leverage.jl) falls below the global run's 3.7. If it
#    does not, these parameters do not move Southern Ocean cloud cover
#    either, and the cloud-fraction shape parameters are the only lever.
# 2. clt improves rather than degrades. It got worse in every global run
#    (0.620 to 0.654). Any monotone fall here means the global degradation
#    was a competition between regions, not a property of clt.
# 3. v_snow stays high. It moved +1.45 prior sigma globally, and if that
#    came from mixed-phase Southern Ocean cloud it should be at least as
#    strong in the band. If it collapses to the prior mean, the global
#    signal came from elsewhere.
# 4. The weather floor is worse here. The band has roughly a tenth of the
#    global points, so internal variability averages down less. Expect the
#    iteration-6 ensemble spread to sit at a higher fraction of iteration
#    1 than the global run's 0.8 to 1.6. If the signal-to-noise falls below
#    the global values, the band is too small for 38-day members and the
#    fix is longer members, not different parameters.
#
# Verdict (run 2026-08-04/05, job 7014203, STOPPED after 3 of 6 iterations
# once the leverage table settled; the remaining 3 would have cost about 35
# GPU-hours to confirm a flat line):
# Leverage ratios are stable and all WORSE than the global vterm run:
#   iteration     1     2     3   global vterm
#   lwp         3.7   3.9   3.8            1.5
#   clt         5.7   6.0   5.8            3.7
#   swcre       3.6   3.7   3.6            2.6
# Residuals moved about 2 percent in total: lwp 0.537 -> 0.526, clt 0.493 ->
# 0.482, swcre 0.325 -> 0.314. Loss 0.0049 -> 0.0046, which back-converts to
# a whitened residual of 0.157 sigma, so the band observation sits even
# further inside its assumed noise than clt did globally (0.35).
# Displacements after 3 iterations, in prior sigma: v_snow -0.24, v_ice
# -0.17, subl_dep +0.05, v_rain +0.04, rain_tau +0.03.
# 1. FALSIFIED. clt's ratio ROSE from the global 3.7 to 5.7, so narrowing to
#    the Southern Ocean made leverage worse, not better, on all three
#    observables. Fixing q_liq removed the largest single lever (lwp's
#    spread fell hardest, 0.456 -> 0.145), and the region is independently
#    less responsive to fall speeds than the tropics and subtropics.
# 2. PASS in sign, negligible in size. clt fell monotonically (0.493,
#    0.488, 0.482) instead of degrading as it did globally, but 0.011 sigma
#    at a ratio of 5.8 is not distinguishable from noise.
# 3. FALSIFIED, and this is the useful part. v_snow moved -0.24 sigma here
#    against +1.45 sigma globally. The global snow signal is NOT a Southern
#    Ocean signal; the band weakly pulls the other way. Look for it where
#    lwp has real leverage, that is, the subtropics and mid-latitudes.
# 4. UNRESOLVED. The G spread was flat across the three iterations (lwp
#    0.145, 0.133, 0.138), which is consistent with a weather floor, but
#    three iterations is too few for the spread-contraction test that
#    measured the floor in the global run.
# KEY METHODOLOGICAL POINT. The residual/spread ratio is invariant to the
# noise model: both terms are divided by sigma, so scaling
# model_error_scale cancels. Lowering the floors would have moved the
# residuals into a comfortable range while leaving leverage exactly as bad.
# Reweighting cannot buy leverage. Do not narrow the region again to chase
# it either; that is what this run tested.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_so.yml",
)

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_so_band")

# Southern Ocean only. Both the observation and the simulation read this.
const LAT_WINDOW = (-65, -45)

# Native 2.5 degree grid: the band is too small to give away 4x the points.
const COARSEN_FACTOR = 1

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

const OBS_DECORRELATION_LENGTH = 4.0e5
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

include(joinpath(@__DIR__, "..", "prior_tools.jl"))

const CALIBRATION_PRIORS = [
    checked_constrained_gaussian("rain_autoconversion_timescale", 1800, 400, 300, 3600),
    checked_constrained_gaussian("sublimation_deposition_timescale", 500, 400, 20, 3600),
    checked_constrained_gaussian("fixed_cloud_ice_terminal_velocity", 0.1, 0.05, 0, 1),
    checked_constrained_gaussian("fixed_rain_terminal_velocity", 5, 1, 0, 10),
    checked_constrained_gaussian("fixed_snow_terminal_velocity", 1, 0.5, 0, 5),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
