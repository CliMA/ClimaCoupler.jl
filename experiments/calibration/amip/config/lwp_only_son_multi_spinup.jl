# LWP-only, SON multi-year minibatch (2006-2010) WITH 7-day spin-up.
#
# Identical to lwp_only_son_multi.jl EXCEPT spinup = 7 days, to isolate the
# spin-up confound: the spinup=0 SON multi-year run converged q_liq to ~3.3e-4,
# well below the single-year (7-day spin-up) ~7e-4. With no spin-up, the SON mean
# includes the transient where the model's prognostic cloud liquid equilibrates
# from an ERA5 IC that has none, biasing LWP. Here each run starts Aug 25
# (= Sep 1 - 7 days) from a dedicated Aug-25 ERA5 IC (son_ic_spinup), spins up 7
# days, then the Sep-Nov sample is taken. Same seasons/years/resolution as the
# spinup=0 run, so the two differ only by spin-up.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration_son_spinup.yml")

const _SON_YEARS = [2006, 2010, 2008, 2007, 2009]
sample_date_ranges = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1))
    for y in _SON_YEARS[1 .+ ((0:9) .% length(_SON_YEARS))]
]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 11, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_son_multi_spinup")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp"],
    minibatch_size = 1,
    n_iterations = 9,
    sample_date_ranges,
    # spinup = 7 days: start = Sep1 - 7d = Aug25 (each year's Aug-25 ERA5 IC);
    # end = Nov1 + 1mo = Dec1. Run Aug25 -> Dec1 (~98 days incl spin-up).
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
