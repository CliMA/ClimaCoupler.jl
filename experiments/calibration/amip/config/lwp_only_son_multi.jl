# LWP-only, SON multi-year minibatch (2006-2010) calibration config.
#
# The original target season (SON), done as a multi-year minibatch to beat down
# single-realization weather noise. Each ensemble member runs one year's SON
# (Sep1 -> Dec1) from that year's downloaded ERA5 IC (Sep-1, 1-deg, in son_ic),
# with the persistence forcing extended to cover the window. Distinct-year
# samples with minibatch_size = 1 IS EKP minibatching: each iteration targets a
# different year, so parameters converge to the multi-year SON climatology
# instead of overfitting one year's weather. No forward_model code change needed.
#
# Note: ~92-day members (~3.2h) mean only ~3 iterations fit the 12h walltime, so
# a single pass samples ~3 of the 5 years. Full 5-year-per-update averaging would
# need minibatch_size>1 (multi-date forward_model + ~25 workers).
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration_son.yml")

# Distinct-year SON samples, cycled to fill n_iterations+1. Ordered to span the
# range early (so the first few iterations see diverse years).
const _SON_YEARS = [2006, 2010, 2008, 2007, 2009]
sample_date_ranges = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1))
    for y in _SON_YEARS[1 .+ ((0:9) .% length(_SON_YEARS))]
]

# Covariance from the interannual spread of SON 2006-2010.
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_son_multi")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp"],
    minibatch_size = 1,
    n_iterations = 9,
    sample_date_ranges,
    # end = Nov1 + 1mo = Dec1 (November complete); spinup = 0 so start = Sep1
    # (each year's Sep-1 ERA5 IC). Forcing covers Aug1..Feb1 (extended).
    extend = Dates.Month(1),
    spinup = Dates.Day(0),
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
