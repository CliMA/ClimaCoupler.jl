# 2-D preprocessing validation at 10 degrees: lwp+pr on a block-averaged
# longitude-latitude grid instead of the zonal mean.
#
# Same parameters, platform, and noise floors as the validated zonal lwp+pr
# calibration. The only change is the spatial reduction: COARSEN_FACTOR = 4
# block-averages the 2.5 degree comparison grid to 10 degrees (36 x 18), which
# is near the decorrelation scale of monthly means and resolves the subtropical
# stratocumulus decks as 2 to 4 pixels.
#
# The known answer is q_liq near 2.8e-4 with rain_tau near 1450 and both
# observables at their noise floors. Reproducing it without collapse validates
# 2-D coarse observations for later phases. The companion config at 15 degrees
# (lwp_pr_2d15.jl) is the fallback if 10 degrees over-informs.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_pr_2d10")

const COARSEN_FACTOR = 4

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "pr"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
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
