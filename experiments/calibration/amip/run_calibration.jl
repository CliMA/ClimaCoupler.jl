using Dates
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "observation_map.jl",
    ),
)

model_interface = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "amip",
    "model_interface.jl",
)

# CALIBRATION CONFIGURATION

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "subseasonal_configs",
    "amip_monthly_calibration.yml",
)

# Calibrate only on Jan 1 2010
sample_date_ranges = [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:6]

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    # Note: Pressure-level variables require model output with
    # pressure_coordinates: true in config
    short_names = ["ta", "hur", "swcre", "lwcre"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir = "/home/ext_nefrathe_caltech_edu/clima/ClimaCoupler.jl/cc/wxquest_v4_final/test_amip_output",
    rng_seed = 42,
)

# Used in generate_observations.jl and observation_map.jl
# Units: Pa (not hPa)
const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]

const CALIBRATION_PRIORS = [
    # Atmospheric parameters
    # PD.constrained_gaussian("entr_inv_tau", 0.002, 0.0015, 0.0, 0.01),

    # PiGroup linear regression coefficients
    # PD.ParameterDistribution(
    #     PD.VectorOfParameterized([PD.Normal(0.0, 5.0), PD.Normal(0.0, 5.0), PD.Normal(0.4, 0.2)]),
    #     repeat([PD.no_constraint()], 3),
    #     "entr_param_vec",
    # ),

    # PD.constrained_gaussian("detr_buoy_coeff", 0.12, 0.06, 0.0, 1.0),
    # PD.constrained_gaussian("detr_vertdiv_coeff", 0.6, 0.25, 0.0, 5.0),
    # PD.constrained_gaussian("EDMF_surface_area", 0.1, 0.03, 0, 1),

    # PD.constrained_gaussian("precipitation_timescale", 600, 300, 100, 2000),
    # PD.constrained_gaussian("precipitation_timescale", 1200, 300, 300, 2400),
    PD.constrained_gaussian("diagnostic_covariance_coeff", 2.1, 0.5, 0.0, 10.0),
    # PD.constrained_gaussian("Tq_correlation_coefficient", 0.0, 0.5, -1.0, 1.0),


    # PD.constrained_gaussian("mixing_length_eddy_viscosity_coefficient", 0.2, 0.1, 0, 1.0),
    # PD.constrained_gaussian("mixing_length_diss_coeff", 0.22, 0.15, 0.0, 10.0),
    # PD.constrained_gaussian("mixing_length_tke_surf_flux_coeff", 8.0, 4.0, 0, 100.0),

    # Land parameters
    # PD.constrained_gaussian("pmodel_cstar", 0.30, 0.15, 0.0, 1.0),
    # PD.constrained_gaussian("leaf_Cd", 0.01, 0.006, 0.0, 0.1),

    # gravity wave parameters
    # PD.constrained_gaussian("nogw_Bt_0", 0.0043, 0.003, 0.001, 0.01),
    # PD.constrained_gaussian("ogw_mountain_height_width_exponent", 0.4, 0.3, 0.0, 1.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)

if abspath(PROGRAM_FILE) == @__FILE__

    observation_vector_filepath = joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "observation_vec.jld2",
    )
    isfile(observation_vector_filepath) || error(
        "Filepath to observation vector is not valid. Update the filepath or generate the observations using generate_observations.jl",
    )
    observation_vector = JLD2.load_object(observation_vector_filepath)

    (; sample_date_ranges, minibatch_size) = CALIBRATE_CONFIG
    obs_series = EKP.ObservationSeries(
        Dict(
            "observations" => observation_vector,
            "names" => [
                string(Dates.year(start_date)) for
                (start_date, stop_date) in sample_date_ranges
            ],
            "minibatcher" => ClimaCalibrate.minibatcher_over_samples(
                length(observation_vector),
                minibatch_size,
            ),
        ),
    )

    rng_seed = CALIBRATE_CONFIG.rng_seed
    rng = Random.MersenneTwister(rng_seed)

    ekp = EKP.EnsembleKalmanProcess(
        obs_series,
        EKP.TransformUnscented(PRIORS, impose_prior = true);
        verbose = true,
        rng,
        scheduler = EKP.DataMisfitController(terminate_at = 1000000),
    )

    if ClimaCalibrate.get_backend() == ClimaCalibrate.DerechoBackend
        backend = ClimaCalibrate.DerechoBackend(
            model_interface,
            verbose = true,
            hpc_kwargs = Dict(
                # Options include "premium", "regular", "economy", "preempt"
                :job_priority => "regular", # {}
                # 720 minutes is 12 hours
                :time => 720,
                :ntasks => 1,
                :cpus_per_task => 12,
                :gpus_per_task => 1,
            ),
        )
    elseif ClimaCalibrate.get_backend() == ClimaCalibrate.GCPBackend
        backend = ClimaCalibrate.GCPBackend(;
            model_interface,
            verbose = true,
            hpc_kwargs = Dict(
                :ntasks => 1,
                :gpus_per_task => 1,
                :cpus_per_task => 12,
                :time => 720,
                :partition => "a3",
            ),
        )
    else
        error("Unsupported backend: $(ClimaCalibrate.get_backend())")
    end

    (; n_iterations, output_dir) = CALIBRATE_CONFIG
    eki = ClimaCalibrate.calibrate(backend, ekp, n_iterations, PRIORS, output_dir;)
end
