# Rerun of lwp_clt_swcre_5d.jl with clt in its own noise group at floor
# 0.35 (was 0.75, shared with swcre). The first run showed clt 2x
# under-floored: its residual sat at 0.33 sigma from iteration 1, it
# steered nothing, and the update traded it away (clt RMSE degraded
# 0.201 -> 0.218, bias flipped +0.078 -> -0.117) while lwp and swcre
# improved. Companion: lwp_clt_swcre_10d_cltf025.jl. Nothing else changes
# (same seed, dates, priors, 6 iterations).
#
# Predictions (falsifiable):
# 1. clt lands near its floor (0.7 to 1.3 sigma), no longer far below.
# 2. clt is not traded away: final clt RMSE does not exceed the initial
#    0.201.
# 3. The cloud-fraction parameters tighten: eps_rel and steepness final
#    spreads come in below the first run's 0.05 x prior values.
# 4. q_liq arbitration: if q_liq stays near 1.7-1.8e-4 with clt now
#    pulling, the low value is real for the 2D configuration; if it moves
#    toward 2.2-2.5e-4, the first run's low value was a clt-neglect
#    artifact.
# 5. No COLLAPSED flag. The extra clt information will contract spreads
#    at least as fast as the first run (which ended 0.03-0.05x prior), so
#    final spreads here carry no uncertainty meaning; only the means count.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir =
    joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_clt_swcre_5d_cltf035")

const COARSEN_FACTOR = 2

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

# lwp/swcre keep the 1.5x-scaled 10 degree floors; clt drops to 0.35
# (1.4x the 10d companion's 0.25) to make it informative.
const OBS_DECORRELATION_LENGTH = 8.0e5
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

const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        5e-4, 3e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1800, 700, 300, 3600),
    PD.constrained_gaussian("snow_autoconversion_timescale", 1800, 700, 300, 3600),
    PD.constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    PD.constrained_gaussian("cloud_fraction_steepness_scale", 0.7, 0.2, 0.1, 2.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
