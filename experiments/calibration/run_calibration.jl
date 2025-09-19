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
include(joinpath(pkgdir(ClimaCoupler), "experiments/calibration/coarse_amip/observation_map.jl"))
model_interface = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "coarse_amip", "model_interface.jl")

years = 2018:2024
sample_date_ranges = [(DateTime(yr, 9, 29), DateTime(yr, 10, 6)) for yr in years]
# sample_date_ranges = [(DateTime(yr, 9, 29), DateTime(yr, 9, 29)) for yr in years]

# sample_date_ranges = [("$yr-12-01", "$yr-12-01") for yr in years]

const CALIBRATE_CONFIG = CalibrateConfig(;
config_file = joinpath(pkgdir(ClimaCoupler), "config/subseasonal_configs/wxquest_diagedmf.yml"),
    short_names = ["tas", "mslp"],
    minibatch_size = 1,
    n_iterations = 7,
    sample_date_ranges,
    extend = Dates.Week(1),
    spinup = Dates.Day(28),
    output_dir = "/glade/derecho/scratch/nefrathe/tmp/output_quick",
    rng_seed = 42,
)


if abspath(PROGRAM_FILE) == @__FILE__

    priors = [
        # PD.constrained_gaussian("albedo_coefficient", 1.0, 0.4, 0, 3),
        PD.constrained_gaussian("mixing_length_diss_coeff", 0.22, 0.07, 0, 1),
        # PD.constrained_gaussian("precipitation_timescale", 919.3827604731249, 150.0, 0, Inf),
        PD.constrained_gaussian("EDMF_surface_area", 0.10928882001604676, 0.03, 0, Inf),
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
    ClimaCalibrate.initialize(ekp, prior, CALIBRATE_CONFIG.output_dir)

    hpc_kwargs = ClimaCalibrate.kwargs(time = 60*5,
        ntasks = 4,
        gpus_per_task = 1,
        cpus_per_task = 4,
        l_job_priority = "premium",
        q = "main")
    exeflags = "--threads=4"
    ekp = ClimaCalibrate.load_latest_ekp(CALIBRATE_CONFIG.output_dir)
    eki = ClimaCalibrate.calibrate(ClimaCalibrate.DerechoBackend,
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir; 
        model_interface, hpc_kwargs, exeflags
    )

end