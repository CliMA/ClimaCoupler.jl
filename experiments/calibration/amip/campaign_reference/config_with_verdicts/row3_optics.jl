# ROW-3 sweep support config: phase3's observables PLUS clt, zonal,
# October 2010. Used by run_parameter_sweep.jl (Row 3, cloud optics) for
# observation generation, wiring preflight, and sweep scoring. Not a
# calibration config; the sweep specification lives in
# run_parameter_sweep.jl.
#
# Priors list the six phase3 parameters (pinned at the phase3 posterior in
# the sweep) plus the two Row-3 optics parameters, so preflight checks the
# optics wiring.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:9]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_row3_optics")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "pr", "swcre", "lwcre", "cl", "clt"],
    minibatch_size = 1,
    n_iterations = 8,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Floors from the validated runs: lwp/pr 0.2 (phase2), swcre 0.5 (irreducible
# brightness bias), lwcre 0.25 (phase2_lwcre025), cl 0.5 (cls zonal).
const OBS_NOISE_GROUPS = [
    (short_names = ["lwp", "pr"], model_error_scale = 0.2),
    (short_names = ["swcre"], model_error_scale = 0.5),
    (short_names = ["lwcre"], model_error_scale = 0.25),
    (short_names = ["cl"], model_error_scale = 0.5),
    (short_names = ["clt"], model_error_scale = 0.3),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# Priors reuse the donor runs' (not their posteriors): the merged run should
# re-derive the answer from the data, not inherit confidence.
const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        3e-4, 1.5e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1450, 500, 300, 3600),
    PD.constrained_gaussian("detr_massflux_vertdiv_coeff", 0.35, 0.15, 0.0, 0.8),
    PD.constrained_gaussian("detr_buoy_coeff", 1.0, 0.6, 0.05, 5.0),
    PD.constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    PD.constrained_gaussian("cloud_fraction_steepness_scale", 0.7, 0.2, 0.1, 2.0),
    PD.constrained_gaussian(
        "prescribed_cloud_droplet_number_concentration",
        1e8, 8e7, 1e6, 1e9,
    ),
    PD.constrained_gaussian("ice_cloud_effective_radius", 25e-6, 1e-5, 5e-6, 1e-4),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
