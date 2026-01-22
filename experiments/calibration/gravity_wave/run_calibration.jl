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
        "experiments/calibration/gravity_wave/observation_map.jl",
    ),
)
model_interface = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "gravity_wave",
    "model_interface.jl",
)

years = 2010:2012
sample_date_ranges = [(DateTime(yr, 2, 1), DateTime(yr, 2, 1)) for yr in years]
const CALIBRATE_CONFIG = CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config/amip_configs/amip_land.yml",
    ),
    short_names = ["ta", "ua", "va"],
    minibatch_size = 1,
    n_iterations = 5,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Month(1),
    output_dir = "output/gravity_wave",
    rng_seed = 42,
)

if abspath(PROGRAM_FILE) == @__FILE__
    priors = [PD.constrained_gaussian("nogw_Bt_0", 0.0043, 0.002, 0.001, 0.01)]
    prior = EKP.combine_distributions(priors)

    observation_vector = JLD2.load_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/gravity_wave/obs_vec.jld2"),
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
        scheduler = EKP.DataMisfitController(terminate_at = 100),
    )

    # Calculate ensemble size for TransformUnscented
    ensemble_size = 2 * length(priors) + 1

    # Add local workers on the same node if not already present
    if nworkers() == 1
        @info "Adding $ensemble_size local workers with GPU assignments..."

        # Add local workers
        addprocs(ensemble_size; exeflags="--project=experiments/ClimaEarth")

        # Assign each worker to a specific GPU
        @everywhere begin
            worker_id = myid() - 1  # Worker IDs start at 2, so subtract 1 to get 0-indexed
            gpu_id = worker_id % 4  # Cycle through 4 GPUs (0, 1, 2, 3)
            ENV["CUDA_VISIBLE_DEVICES"] = string(gpu_id)
            @info "Worker $(myid()) assigned to GPU $gpu_id"
        end
    end

    # Load model interface on all workers
    api_file = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "api.jl")
    @everywhere include($api_file)
    @everywhere const CALIBRATE_CONFIG = $CALIBRATE_CONFIG
    @everywhere using Dates
    @everywhere include($model_interface)

    # Use WorkerBackend for single-node calibration
    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.WorkerBackend(),
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir;
    )
end
