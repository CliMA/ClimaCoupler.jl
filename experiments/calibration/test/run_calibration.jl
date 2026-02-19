using Dates
using Distributed
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaComms
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

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

years = 2018:2021
sample_date_ranges = [(DateTime(yr, 9, 1), DateTime(yr, 9, 1)) for yr in years]
const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config/subseasonal_configs/wxquest_diagedmf.yml",
    ),
    short_names = ["hfls", "hfss", "rsus", "rlus"],
    minibatch_size = 1,
    n_iterations = 3,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Month(0),
    output_dir = "output/subseasonal",
    rng_seed = 42,
)

if abspath(PROGRAM_FILE) == @__FILE__
    priors = [PD.constrained_gaussian("precipitation_timescale", 600, 300, 100, 1000)]
    prior = EKP.combine_distributions(priors)

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
        scheduler = EKP.DataMisfitController(terminate_at = 100),
    )

    backend = ClimaCalibrate.WorkerBackend()

    eki = ClimaCalibrate.calibrate(
        backend,
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir;
    )
end
