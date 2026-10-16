# Stage 1 of the EarthCARE cloud calibration: frozen water path against EarthCARE, with
# CERES cloud radiative effect and GPCP precipitation, over four 1M snow and rain
# parameters. See docs/src/earthcare_calibration.md.
#
# The EarthCARE record is Dec 2025 to May 2026. Its months are relabelled to 2010 in the
# target file so each target pairs with a 2010 model run; the comparison is
# climatological. The noise covariance comes from ERA5 cloud ice plus snow over 2001 to
# 2010, one covariance per calendar month.

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Targets alternate April and January 2010 (EarthCARE April and January 2026); iteration
# `i` runs and is scored against entry `i`.
sample_date_ranges =
    [(Dates.DateTime(2010, m, 1), Dates.DateTime(2010, m, 1)) for m in repeat([4, 1], 3)]
covariance_date_ranges = [
    (Dates.DateTime(year, m, 1), Dates.DateTime(year, m, 1)) for m in (1, 4) for
    year in 2001:2010
]
const COVARIANCE_BY_MONTH = true

# Member output is tens of GB, so it goes to local scratch rather than the NFS home.
output_dir =
    get(ENV, "CALIBRATION_OUTPUT_DIR", "/scratch/nefrathe/amip_calibration_earthcare_fwp")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["fwp", "swcre", "lwcre", "pr"],
    minibatch_size = 1,
    n_iterations = 2,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# Model-error variance per observable, in the order of `short_names`. The frozen water
# path carries the radar retrieval uncertainty (about 40 percent), CRE the values of the
# zonal_cloud config, precipitation the GPCP uncertainty.
const NOISE_BETA = [0.15, 0.074, 0.074, 0.10]

# Target and covariance files, built by ~/earthcare_eval/build_targets.py until they
# become ClimaArtifacts.
const EARTHCARE_TARGETS_NC = get(
    ENV,
    "EARTHCARE_TARGETS_NC",
    "/home/nefrathe/earthcare_eval/targets/earthcare_targets.nc",
)
const ERA5_FWP_NC =
    get(ENV, "ERA5_FWP_NC", "/home/nefrathe/earthcare_eval/targets/era5_fwp_2001_2010.nc")
const EXTRA_DATA_LOADERS = [
    "fwp" => CalibrationTools.LocalNetCDFDataLoader(
        EARTHCARE_TARGETS_NC,
        ["fwp" => "fwp", "cl1km" => "cl1km", "cl3km" => "cl3km", "cl8km" => "cl8km"],
    ),
]
const COVARIANCE_DATA_LOADERS =
    ["fwp" => CalibrationTools.LocalNetCDFDataLoader(ERA5_FWP_NC, ["fwp" => "fwp"])]

# 1M microphysics parameters that set how much frozen water is held as snow and how fast
# it falls out, plus the warm-rain timescale that shares the precipitation path.
# Snow sediments with the Chen2022 fall speed by default (`snow_terminal_velocity`), whose
# coefficients are set by the snow apparent density; `fixed_snow_terminal_velocity` and
# `snow_aspect_ratio` have no effect in this configuration.
const CALIBRATION_PRIORS = [
    PD.constrained_gaussian("snow_autoconversion_timescale", 1360, 500, 200, 6000),
    PD.constrained_gaussian("cloud_ice_snow_collision_efficiency", 0.1, 0.07, 0.005, 1.0),
    PD.constrained_gaussian("snow_apparent_density", 100, 40, 30, 400),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1200, 300, 300, 3600),
]
const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
