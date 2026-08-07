# LWP-only, OND 2010 (3 months: Oct+Nov+Dec) calibration config.
#
# 3-month refinement of lwp_only.jl to raise signal-to-noise (3x the monthly data
# concatenated -> ~1/sqrt(3) tighter parameter posterior and a smoother error
# trajectory). OND rather than SON because ERA5 initial conditions only exist on
# sparse dates: the run start (sample_start - spinup) must match an available IC,
# and 2010 has no Aug/early-Sep IC. OND starts from the 2010-09-24 IC with the
# usual 7-day spinup; SON would have forced a Jul-1 start (~153-day runs).
#
# The recipe windows [start, end] and flattens all monthly slices in it (it does
# NOT average), so an Oct1..Dec1 range yields the Oct/Nov/Dec monthly means
# stacked. Select via ENV["CALIBRATION_CONFIG"] = this file.

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Sample = OND 2010. Range spans month-start labels Oct1, Nov1, Dec1 (monthly
# slices are start-labeled). Length must be >= n_iterations + 1.
sample_date_ranges = [
    (Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 12, 1)) for _ in 1:7
]

# Covariance from the same OND season across 2006-2010. Must contain the sample.
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 12, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_ond")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp"],
    minibatch_size = 1,
    # Each member runs ~Sep24 -> Jan1 (~99 days, ~2.5x the 1-month cost). Params
    # converged by ~iter 6 in the 1-month run; richer data converges at least as
    # fast, so 6 iterations (walltime permitting ~4-5) is enough.
    n_iterations = 6,
    sample_date_ranges,
    # extend = Month(1): end = Dec1 + 1mo = Jan1, so December is a complete
    # monthly mean. spinup = 7d puts the start at 2010-09-24 (a valid ERA5 IC).
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
