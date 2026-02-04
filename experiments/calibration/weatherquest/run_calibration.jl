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
        "experiments/calibration/coarse_amip/observation_map.jl",
    ),
)
model_interface = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "coarse_amip",
    "model_interface.jl",
)

years = 2018:2024
sample_date_ranges = [(DateTime(yr, 10, 1), DateTime(yr, 10, 1)) for yr in years]
const CALIBRATE_CONFIG = CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config/subseasonal_configs/wxquest_diagedmf.yml",
    ),
    short_names = ["tas", "tas - ta", "hfls", "hfss", "rsns", "rlns"],
    minibatch_size = 1,
    n_iterations = 7,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Month(1),
    output_dir = "/glade/derecho/scratch/nefrathe/tmp/output_oceanmask",
    rng_seed = 42,
)

if abspath(PROGRAM_FILE) == @__FILE__

    priors = [
        PD.constrained_gaussian("bucket_soil_heat_capacity", 5e6, 2.5e6, 1e6, 10e6),
        PD.constrained_gaussian("bucket_capacity_fraction", 2, 1.5, 0, 10),
        PD.constrained_gaussian("mixing_length_diss_coeff", 0.25, 0.1, 0, 1),
        PD.constrained_gaussian("EDMF_surface_area", 0.14, 0.05, 0, 1),
        PD.constrained_gaussian(
            "mixing_length_eddy_viscosity_coefficient",
            0.22,
            0.1,
            0,
            1,
        ),
        PD.constrained_gaussian("mixing_length_tke_surf_flux_coeff", 7, 3, 0, 20),
    ]
    prior = EKP.combine_distributions(priors)

    observation_vector = JLD2.load_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/weatherquest_obs_vec.jld2"),
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
        scheduler = EKP.DataMisfitController(terminate_at = 100),
    )

    hpc_kwargs = ClimaCalibrate.kwargs(
        time = 60 * 6,
        ntasks = 2,
        gpus_per_task = 1,
        cpus_per_task = 4,
        q = "main",
    )
    # backend = ClimaCalibrate.DerechoBackend(;model_interface, hpc_kwargs, verbose=true)
    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.DerechoBackend,
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir;
        model_interface,
        hpc_kwargs,
        exeflags = "--threads=4",
    )
end
