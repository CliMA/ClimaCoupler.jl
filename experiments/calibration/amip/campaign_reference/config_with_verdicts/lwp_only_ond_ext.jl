# LWP-only, OND 2010 (Oct+Nov+Dec) with EXTENDED forcing.
#
# True 3-month run enabled by extending the wxquest SST/SIC/albedo forcing to
# cover through Jan 1 2011 (downloaded from ERA5 via CDS and reprocessed with
# WeatherQuest's preprocess_sst/sic/albedo — see wxquest_ic_extended/). Uses a
# custom coupler yml (amip_calibration_ond.yml) that sets era5_initial_condition_dir
# to that extended directory.
#
# Starts from the Oct-1 ERA5 IC (spinup = 0, start = Oct 1); the run goes
# Oct 1 -> Jan 1 (~99 days), fully within the extended Sep1..Jan1 forcing window.

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration_ond.yml")

# Sample = OND 2010 (Oct/Nov/Dec monthly means).
sample_date_ranges = [
    (Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 12, 1)) for _ in 1:7
]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 12, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_ond_ext")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp"],
    minibatch_size = 1,
    # ~99-day members (~2.5x baseline), so ~3-4 iterations fit the 12h walltime.
    n_iterations = 6,
    sample_date_ranges,
    # end = Dec1 + 1mo = Jan1 (December complete); spinup = 0 so start = Oct1 (the
    # Oct-1 ERA5 IC date). Forcing covers Sep1..Jan1 (extended).
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
