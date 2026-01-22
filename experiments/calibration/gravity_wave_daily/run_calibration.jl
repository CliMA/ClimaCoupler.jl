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
        "experiments/calibration/gravity_wave_daily/observation_map.jl",
    ),
)
model_interface = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "gravity_wave_daily",
    "model_interface.jl",
)

# Daily sample date ranges (start_date, end_date) for calibration
# Each range corresponds to a single day for daily calibration testing
sample_date_ranges = [
    (DateTime(2023, 1, 15), DateTime(2023, 1, 15)),
    (DateTime(2023, 1, 16), DateTime(2023, 1, 16)),
]

# Directory containing ERA5 daily observation files
const ERA5_OBS_DIR = "/glade/campaign/univ/ucit0011/cchristo/wxquest_data/daily_weekly_stats/"

const CALIBRATE_CONFIG = CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config/amip_configs/amip_land_daily.yml",
    ),
    short_names = ["tas", "mslp"],  # 2m temperature + mean sea level pressure
    minibatch_size = 1,
    n_iterations = 2,  # Reduced for testing
    sample_date_ranges,
    extend = Dates.Day(1),
    spinup = Dates.Day(0),
    output_dir = "output/gravity_wave_daily",
    obs_dir = ERA5_OBS_DIR,
    rng_seed = 42,
)

if abspath(PROGRAM_FILE) == @__FILE__
    # Priors for calibration parameters
    priors = [
        # Gravity wave momentum flux parameter: mean=0.0043, std=0.002, bounds=[0.001, 0.01]
        PD.constrained_gaussian("nogw_Bt_0", 0.0043, 0.002, 0.001, 0.01),
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
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/gravity_wave_daily/obs_vec.jld2"),
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
