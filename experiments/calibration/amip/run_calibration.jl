using Dates
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

model_interface_filepath = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "amip",
    "model_interface.jl",
)
include(model_interface_filepath)

# Choose which calibration config to use
config_dir = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip", "config")
default_config_path = joinpath(config_dir, "pressure_levels.jl")

test_calibration_config_path = joinpath(config_dir, "pipeline_test.jl")
const TEST_CALIBRATION = haskey(ENV, "TEST_CALIBRATION")

# Allow selecting an alternate config via ENV (e.g. the LWP-only calibration),
# so multiple calibration setups can share this driver. Falls back to the test
# config or the default lwp+cl config.
config_path = if haskey(ENV, "CALIBRATION_CONFIG")
    ENV["CALIBRATION_CONFIG"]
elseif TEST_CALIBRATION
    test_calibration_config_path
else
    default_config_path
end

@info "Using calibration configuration in: $config_path"
include(config_path)

(; output_dir) = CALIBRATE_CONFIG
isdir(output_dir) || mkdir(output_dir)

if abspath(PROGRAM_FILE) == @__FILE__
    # The observation vector lives in the config's output_dir so that different
    # calibration setups (e.g. lwp+cl vs LWP-only) keep independent observations.
    observation_vector_filepath = joinpath(output_dir, "observation_vec.jld2")
    isfile(observation_vector_filepath) || error(
        "Filepath to observation vector is not valid. Update the filepath or generate the observations using generate_observations.jl",
    )
    observation_vector = JLD2.load_object(observation_vector_filepath)

    (; sample_date_ranges, minibatch_size) = CALIBRATE_CONFIG
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

    (; rng_seed) = CALIBRATE_CONFIG
    rng = Random.MersenneTwister(rng_seed)

    # Use a gentle fixed timestep instead of DataMisfitController. DMC has no
    # step-size knob and, with an over-informative observation, collapsed the
    # ensemble covariance ~3 orders of magnitude in a single step — pinning the
    # parameters near the prior mean before any weak signal could be exploited.
    # A small fixed Δt spreads the update over many iterations: accumulated
    # pseudo-time T reaches ~1 (the full posterior) over 1/Δt steps, so set
    # n_iterations ≈ 1/Δt in the config (Δt = 0.1 → n_iterations = 10).
    ekp = EKP.EnsembleKalmanProcess(
        obs_series,
        EKP.TransformUnscented(PRIORS, impose_prior = true);
        verbose = true,
        rng,
        scheduler = EKP.DefaultScheduler(0.1),
    )

    coupler_model_interface = CouplerModelInterface(CALIBRATE_CONFIG)
    (; n_iterations, output_dir) = CALIBRATE_CONFIG

    # Async workers (ClimaCalibrate ne/async branch): persistent GPU workers run
    # ensemble members concurrently and are reused across iterations, instead of
    # submitting one PBS job per member per iteration. See experiments/AMIP
    # Project.toml [sources] which pins ClimaCalibrate to the ne/async branch.
    walltime = TEST_CALIBRATION ? 60 : 720
    ClimaCalibrate.add_workers(
        EKP.get_N_ens(ekp);
        device = :gpu,
        time = walltime,
        env = ["CLIMACOMMS_CONTEXT" => "SINGLETON", "CLIMACOMMS_DEVICE" => "CUDA"],
    )
    # Load the model interface on every worker, including any that join later.
    ClimaCalibrate.@worker_setup include($model_interface_filepath)

    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.WorkerBackend(),
        ekp,
        coupler_model_interface,
        n_iterations,
        PRIORS,
        output_dir,
    )
end
