# Define which coupler file to use
config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Calibrate only on Jan 1 2010
sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:6]

# On Derecho, it is preferable to save the calibration output to the scratch
# directory (e.g. "/glade/derecho/scratch")
output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    # Note: Pressure-level variables require model output with
    # pressure_coordinates: true in config
    short_names = ["ta", "hur"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Used in generate_observations.jl and observation_map.jl
# Units: Pa (not hPa)
const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]

# To disable normalization, update generate_observations.jl to not apply the
# normalization. You may want to do the same in the observation map as well.
const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

const CALIBRATION_PRIORS = [
    PD.constrained_gaussian("precipitation_timescale", 1200, 300, 300, 2400),
    PD.constrained_gaussian("Tq_correlation_coefficient", 0.4, 0.4, -1.0, 1.0),
    PD.constrained_gaussian("mixing_length_eddy_viscosity_coefficient", 0.2, 0.1, 0, 1.0),
    PD.constrained_gaussian("mixing_length_diss_coeff", 0.22, 0.15, 0.0, 10.0),
    PD.constrained_gaussian("mixing_length_tke_surf_flux_coeff", 8.0, 4.0, 0, 100.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
