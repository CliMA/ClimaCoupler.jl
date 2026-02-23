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

# Sample date range for calibration (repeated for each iteration)
# const BASE_DATE_RANGE = (DateTime(2010, 1, 1), DateTime(2010, 1, 7))
const BASE_DATE_RANGE = (DateTime(2010, 1, 1), DateTime(2010, 1, 1))
const N_ITERATIONS = 3

# Repeat the date range for each iteration so we can reuse subseasonal's forward_model
# which indexes by sample_date_ranges[iter + 1]
sample_date_ranges = fill(BASE_DATE_RANGE, N_ITERATIONS)

# Directory containing ERA5 weekly observation files
const ERA5_OBS_DIR = "/glade/campaign/univ/ucit0011/cchristo/wxquest_data/daily_weekly_stats/weekly"

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config/subseasonal_configs/wxquest_diagedmf_weekly_calibration.yml",
    ),
    # short_names = ["tas", "mslp", "rsut", "rlut"],
    short_names = ["tas", "rlut"],
    minibatch_size = 1,
    n_iterations = N_ITERATIONS,
    sample_date_ranges,
    extend = Dates.Day(1),
    spinup = Dates.Day(0),
    output_dir = "/glade/derecho/scratch/cchristo/calibration/exp30",
    rng_seed = 42,
)

if abspath(PROGRAM_FILE) == @__FILE__
    # Load priors from shared config (single source of truth)
    include(joinpath(@__DIR__, "calibration_setup.jl"))
    prior = CALIBRATION_PRIOR

    observation_vector = JLD2.load_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal_weekly/obs_vec.jld2"),
    )

    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    minibatch_size = CALIBRATE_CONFIG.minibatch_size
    # Structure observations into an ObservationSeries
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

    backend = ClimaCalibrate.DerechoBackend(
        model_interface = model_interface,
        verbose = true,
        hpc_kwargs = Dict(
            :job_priority => "regular", # {"premium", "regular", "economy", "preempt"}
            :time => 720,           # 12 hours in minutes
            :ntasks => 1,
            :cpus_per_task => 12,
            :gpus_per_task => 1,
        ),
    )

    # backend = ClimaCalibrate.WorkerBackend()

    eki = ClimaCalibrate.calibrate(
        backend,
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir;
    )
end
