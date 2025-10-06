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
include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/seasonal/observation_map.jl"))
model_interface = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "seasonal", "model_interface.jl")
include(model_interface)
years = 2010
sample_date_ranges = repeat([(DateTime(yr, 3, 1), DateTime(yr, 5, 1)) for yr in years], 5)

const CALIBRATE_CONFIG = CalibrateConfig(;
config_file = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/seasonal/amip_config.yml"),
    short_names = ["sw_cre", ],
    minibatch_size = 1,
    n_iterations = 5,
    sample_date_ranges,
    extend = Dates.Week(7),
    spinup = Dates.Day(14),
    # TODO: Use this in the model_interface
    output_dir = "/glade/derecho/scratch/zhaoyi/tmp/amip_calibration",
    rng_seed = 42,
)


if abspath(PROGRAM_FILE) == @__FILE__

    priors = [
        EKP.constrained_gaussian("prescribed_cloud_droplet_number_concentration", 3e8, 5e7, 1e7, 1e9)
        #EKP.constrained_gaussian("ice_cloud_effective_radius", 4e-5, 5e-6, 5e-6, 9e-5)
    ]
    prior = EKP.combine_distributions(priors)

    observation_vector =
        JLD2.load_object(joinpath(pkgdir(ClimaCoupler),"experiments/calibration/weatherquest_obs_vec.jld2"))

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
        scheduler = EKP.DataMisfitController(terminate_at = 100),
    )

    # Using ETKI
    # initial_ensemble = EKP.construct_initial_ensemble(rng, prior, 10)
    # ekp = EKP.EnsembleKalmanProcess(
    #     initial_ensemble,
    #     obs_series,
    #     EKP.TransformInversion(prior);
    #     verbose = true,
    #     rng,
    #     scheduler = EKP.DataMisfitController(terminate_at = 100),
    # )

    hpc_kwargs = ClimaCalibrate.kwargs(time = 60*12,
        ntasks = 2,
        gpus_per_task = 1,
        cpus_per_task = 4,
        l_job_priority = "premium",
        q = "main")
    exeflags = "--threads=4"
    eki = ClimaCalibrate.calibrate(ClimaCalibrate.DerechoBackend,
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir; 
        model_interface, hpc_kwargs, exeflags
    )

end
