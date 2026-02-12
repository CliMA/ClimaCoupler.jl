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
        "experiments/calibration/subseasonal_weekly/observation_map.jl",
    ),
)
model_interface = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "subseasonal_weekly",
    "model_interface.jl",
)


# CALIBRATION CONFIGURATION

const BASE_DATE_RANGE = (DateTime(2010, 1, 1), DateTime(2010, 1, 31))
const N_ITERATIONS = 6


# --- 1-day test run ---
# const BASE_DATE_RANGE = (DateTime(2010, 1, 1), DateTime(2010, 1, 1))
# const N_ITERATIONS = 3
# extend = Dates.Day(1)
# spinup = Dates.Day(0)

# Repeat the date range for each iteration until we're using multiple months
sample_date_ranges = fill(BASE_DATE_RANGE, N_ITERATIONS)

# Directory containing ERA5 weekly observation files (not used for CERES-only runs)
const ERA5_OBS_DIR = "/glade/campaign/univ/ucit0011/cchristo/wxquest_data/daily_weekly_stats/weekly"

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config/subseasonal_configs/wxquest_diagedmf_weekly_calibration.yml",
    ),
    # short_names = ["rsut", "rlut"],
    # Note: Pressure-level variables require model output with pressure_coordinates: true
    short_names = [
        "ta_850hPa",
        "ta_500hPa",
        "ta_200hPa",
        "hur_850hPa",
        "hur_500hPa",
        "hur_200hPa",
    ],
    minibatch_size = 1,
    n_iterations = N_ITERATIONS,
    sample_date_ranges,
    extend = Dates.Day(1),
    spinup = Dates.Day(7),
    output_dir = "/glade/derecho/scratch/cchristo/calibration/exp35",
    rng_seed = 42,
)

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "calibration_setup.jl"))
    prior = CALIBRATION_PRIOR

    observation_vector = JLD2.load_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal_weekly/obs_vec.jld2"),
    )

    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    minibatch_size = CALIBRATE_CONFIG.minibatch_size
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
        EKP.TransformUnscented(prior, impose_prior = true);
        verbose = true,
        rng,
        scheduler = EKP.DataMisfitController(terminate_at = 1000000),
    )

    # backend = ClimaCalibrate.ClimaGPUBackend(;
    #     hpc_kwargs = ClimaCalibrate.kwargs(gpus = 1, time = 60 * 12),
    #     model_interface,
    #     verbose = true,
    # )
    if ClimaCalibrate.get_backend() == ClimaCalibrate.DerechoBackend
        backend = ClimaCalibrate.DerechoBackend(;
            model_interface,
            verbose = true,
            hpc_kwargs = Dict(
                :job_priority => "regular", # {"premium", "regular", "economy", "preempt"}
                :time => 720,           # 12 hours in minutes
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
                :gpus_per_task => 1,
                :cpus_per_task => 12,
                :time => 720,
                :partition => "a3",
            ),
        )
    else
        error("Unsupported backend: $(ClimaCalibrate.get_backend())")
    end

    eki = ClimaCalibrate.calibrate(
        backend,
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir;
    )
end
