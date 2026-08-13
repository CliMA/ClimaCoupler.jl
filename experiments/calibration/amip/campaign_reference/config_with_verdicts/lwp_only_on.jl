# LWP-only, ON 2010 (2 months: Oct+Nov) calibration config.
#
# 2-month refinement of lwp_only.jl to raise signal-to-noise (~sqrt(2) ≈ 1.4x)
# without any new data. A 3-month (OND/SON) run is impossible here: the
# subseasonal SST/sea-ice forcing (wxquest artifact) only spans IC-1mo..IC+2mo,
# and the Sep-24 IC's forcing ends 2010-11-24 (which killed the OND attempt).
# The Oct-1 IC's forcing covers Sep 1 -> Dec 1, which is exactly enough for an
# Oct+Nov run.
#
# Because the run must start on an available ERA5 IC date, we start from the
# Oct-1 IC with spinup = 0 (start = Oct 1). The first ~2-3 days of October
# therefore carry cloud spin-up transients — a ~4% contamination of the 61-day
# Oct+Nov mean, acceptable for LWP. Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Sample = ON 2010. Range spans month-start labels Oct1, Nov1 so the time window
# captures the Oct and Nov monthly means. Length must be >= n_iterations + 1.
sample_date_ranges = [
    (Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 11, 1)) for _ in 1:7
]

# Covariance from the same ON season across 2006-2010. Must contain the sample.
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 11, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_on")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp"],
    minibatch_size = 1,
    # Each member runs Oct1 -> Dec1 (~61 days, ~1.6x the 1-month cost). ~2h/iter
    # warm, so ~5-6 iterations fit the 12h walltime.
    n_iterations = 6,
    sample_date_ranges,
    # extend = Month(1): end = last(range)=Nov1 + 1mo = Dec1, so November is a
    # complete monthly mean and the run stays within the Oct-1 IC forcing window.
    extend = Dates.Month(1),
    # spinup = 0: start = Oct1 - 0 = Oct1, matching the Oct-1 ERA5 IC date (the
    # forward model requires start to be an available IC date).
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
