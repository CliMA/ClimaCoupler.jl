using Dates
using Distributed
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaComms
import ClimaCoupler
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2


include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "api.jl"))
include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments/calibration/subseasonal/observation_map.jl",
    ),
)
model_interface = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "subseasonal",
    "model_interface.jl",
)

# Weekly sample date ranges (start_date, end_date) for calibration
# Each range corresponds to a 7-day period matching the ERA5 weekly files
sample_date_ranges = [
    # (DateTime(2023, 1, 15), DateTime(2023, 1, 21)),
    # test a day 
    (DateTime(2023, 1, 15), DateTime(2023, 1, 15)),
    # Add more weekly ranges as needed:
    # (DateTime(2023, 1, 22), DateTime(2023, 1, 28)),
]

# Directory containing ERA5 weekly observation files
const ERA5_OBS_DIR = "/glade/campaign/univ/ucit0011/cchristo/wxquest_data/daily_weekly_stats/weekly"

const CALIBRATE_CONFIG = CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config/subseasonal_configs/wxquest_diagedmf.yml",
    ),
    short_names = ["tas"],  # Start with tas only
    # short_names = ["tas", "mslp", "pr"],  # Uncomment to add more variables
    minibatch_size = 1,
    n_iterations = 8,
    sample_date_ranges,
    extend = Dates.Day(1),  # Add 1 day so simulation covers full 7-day diagnostic period
    spinup = Dates.Day(0),
    output_dir = "output/subseasonal/exp4",
    obs_dir = ERA5_OBS_DIR,
    rng_seed = 42,
)

if abspath(PROGRAM_FILE) == @__FILE__
    # Priors for calibration parameters
    priors = [
        # Inverse entrainment timescale: mean=0.002, std=0.001, bounds=[0.0, 0.01]
        PD.constrained_gaussian("entr_inv_tau", 0.002, 0.0015, 0.0, 0.01),
        # Precipitation timescale: mean=600, std=300, bounds=[100, 1000]
        # PD.constrained_gaussian("precipitation_timescale", 600, 300, 100, 1000),
    ]
    prior = EKP.combine_distributions(priors)
    ensemble_size = 2 * length(priors) + 1  # TransformUnscented ensemble size

    # Add PBS workers with GPU resources on Derecho
    if nworkers() == 1
        @info "Adding PBS workers with GPU resources..."
        addprocs(
            ClimaCalibrate.PBSManager(ensemble_size);
            q = "main",
            A = "UCIT0011",
            l_select = "1:ncpus=12:ngpus=1",
            l_walltime = "06:00:00",
        )
    end

    # Load api.jl on all workers first (defines CalibrateConfig type)
    api_file = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "api.jl")
    @everywhere include($api_file)
    
    # Share CALIBRATE_CONFIG with workers and load model interface
    @everywhere const CALIBRATE_CONFIG = $CALIBRATE_CONFIG
    @everywhere include($model_interface)

    observation_vector = JLD2.load_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/obs_vec.jld2"),
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
        scheduler = EKP.DataMisfitController(terminate_at = 1000),
    )

    # Use WorkerBackend with PBS workers
    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.WorkerBackend(),
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir;
    )
end
