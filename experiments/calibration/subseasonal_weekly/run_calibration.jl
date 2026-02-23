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

# ==========================================================================
# CALIBRATION CONFIGURATION
# ==========================================================================
# Monthly calibration with spinup (using Jan 1 IC file):
#   - BASE_DATE_RANGE = (Jan 8, Jan 8): calibration period START (after spinup)
#   - spinup = 7 days: model starts 7 days BEFORE sample_date_range = Jan 1 (IC date)
#   - extend = 21 days: model runs 21 days AFTER sample_date_range = Jan 29
#   
# Timeline: Model runs Jan 1 -> Jan 29 (28 days total)
#           Calibration uses Jan 8 -> Jan 29 (21 days = 3 weeks)
#           Compared against January CERES monthly mean
#
const BASE_DATE_RANGE = (DateTime(2010, 1, 8), DateTime(2010, 1, 8))
const N_ITERATIONS = 6

# 1-day test run ---
# const BASE_DATE_RANGE = (DateTime(2010, 1, 1), DateTime(2010, 1, 1))
# const N_ITERATIONS = 3
# extend = Dates.Day(1)
# spinup = Dates.Day(0)

# Repeat the date range for each iteration so we can reuse subseasonal's forward_model
# which indexes by sample_date_ranges[iter + 1]
sample_date_ranges = fill(BASE_DATE_RANGE, N_ITERATIONS)

# Directory containing ERA5 weekly observation files (not used for CERES-only runs)
const ERA5_OBS_DIR = "/glade/campaign/univ/ucit0011/cchristo/wxquest_data/daily_weekly_stats/weekly"

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config/subseasonal_configs/wxquest_diagedmf_weekly_calibration.yml",
    ),
    short_names = ["rsut", "rlut"],
    minibatch_size = 1,
    n_iterations = N_ITERATIONS,
    sample_date_ranges,
    # Monthly run: 7-day spinup + 21-day calibration = 28 days total
    # Model starts at (Jan 8 - 7 days) = Jan 1, ends at (Jan 8 + 21 days) = Jan 29
    extend = Dates.Day(21),
    spinup = Dates.Day(7),
    # 1-day test run ---
    # extend = Dates.Day(1),
    # spinup = Dates.Day(0),
    output_dir = "/glade/derecho/scratch/cchristo/calibration/exp34",
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
