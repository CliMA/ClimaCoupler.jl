# LWP-only, SON 2010 (3 months: Sep+Oct+Nov) calibration config.
#
# Refinement of lwp_only.jl to raise signal-to-noise. The single-October run
# converged the parameters consistently but with a noisy (weather-dominated)
# error trajectory. Using three months of monthly means concatenated into one
# observation triples the (semi-independent) data, shrinking the parameter
# posterior ~1/sqrt(3) and smoothing the per-iteration error. The recipe windows
# [start, end] and flattens all monthly slices in it (it does NOT average), so a
# Sep1..Nov1 range yields the Sep/Oct/Nov monthly means stacked.
#
# Select via ENV["CALIBRATION_CONFIG"] = this file.

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Sample = SON 2010. The range spans month-start labels Sep1, Oct1, Nov1 so the
# time window captures all three monthly means (monthly slices are start-labeled,
# as confirmed by the single-month config selecting October via (Oct1, Oct1)).
# Length must be >= n_iterations + 1.
sample_date_ranges = [
    (Dates.DateTime(2010, 9, 1), Dates.DateTime(2010, 11, 1)) for _ in 1:7
]

# Covariance from the same SON season across 2006-2010 (interannual spread of the
# stacked 3-month vector). Must contain the sample range.
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_son")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp"],
    minibatch_size = 1,
    # Each member now runs ~Aug25 -> Dec1 (~98 days, ~2.5x the 1-month cost), so
    # fewer iterations fit the 12h walltime. Params converged by ~iter 6 in the
    # 1-month run and the richer data should converge at least as fast.
    n_iterations = 6,
    sample_date_ranges,
    # extend = Month(1): end = last(range)=Nov1 + 1mo = Dec1, so November is a
    # complete monthly mean. spinup pushes the start back to ~Aug25.
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        5e-4, 3e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1800, 700, 300, 3600),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
