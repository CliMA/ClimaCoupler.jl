using Dates
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

using Distributed

model_interface_filepath = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "amip",
    "model_interface.jl",
)
include(model_interface_filepath)

# Choose which calibration config to use
config_path = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip", "config")
calibration_config_path = joinpath(config_path, "pressure_levels.jl")

test_calibration_config_path = joinpath(config_path, "pipeline_test.jl")
const TEST_CALIBRATION = haskey(ENV, "TEST_CALIBRATION")

config_path = TEST_CALIBRATION ? test_calibration_config_path : calibration_config_path

@info "Using calibration configuration in: $config_path"
include(config_path)

(; output_dir) = CALIBRATE_CONFIG
isdir(output_dir) || mkdir(output_dir)

if abspath(PROGRAM_FILE) == @__FILE__
    observation_vector_filepath = joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "observation_vec.jld2",
    )
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

    ekp = EKP.EnsembleKalmanProcess(
        obs_series,
        EKP.TransformUnscented(PRIORS, impose_prior = true);
        verbose = true,
        rng,
        scheduler = EKP.DataMisfitController(terminate_at = 1000000),
    )

    coupler_model_interface = CouplerModelInterface(CALIBRATE_CONFIG)

    if false
        addprocs(ClimaCalibrate.SlurmManager())
        @everywhere include($model_interface_filepath)
        (; n_iterations, output_dir) = CALIBRATE_CONFIG
        ClimaCalibrate.calibrate(
            ClimaCalibrate.WorkerBackend(),
            ekp,
            coupler_model_interface,
            n_iterations,
            PRIORS,
            output_dir,
        )
        return nothing
    elseif ClimaCalibrate.get_backend() == ClimaCalibrate.DerechoBackend
        backend = ClimaCalibrate.DerechoBackend(;
            directives = [
                # Options include "premium", "regular", "economy", "preempt"
                :job_priority => "regular",
                # 720 minutes is 12 hours
                :time => 720,
                :ntasks => 1,
                :cpus_per_task => 12,
                :gpus_per_task => 1,
            ],
            modules = ["climacommon"],
            env_vars = ["CLIMACOMMS_CONTEXT" => "SINGLETON", "CLIMACOMMS_DEVICE" => "CUDA"],
        )
    elseif ClimaCalibrate.get_backend() == ClimaCalibrate.GCPBackend
        backend = ClimaCalibrate.GCPBackend(;
            directives = [
                :ntasks => 1,
                :gpus_per_task => 1,
                :cpus_per_task => 12,
                :time => 720,
                :partition => "a3mega",
            ],
            modules = ["climacommon"],
            env_vars = ["CLIMACOMMS_CONTEXT" => "SINGLETON", "CLIMACOMMS_DEVICE" => "CUDA"],
        )
    elseif ClimaCalibrate.get_backend() == ClimaCalibrate.ClimaGPUBackend
        backend = ClimaCalibrate.ClimaGPUBackend(;
            directives = [:ntasks => 1, :cpus_per_task => 6, :time => 60],
            modules = ["climacommon"],
            env_vars = ["CLIMACOMMS_CONTEXT" => "SINGLETON", "CLIMACOMMS_DEVICE" => "CPU"],
        )
    else
        error("Unsupported backend: $(ClimaCalibrate.get_backend())")
    end

    (; n_iterations, output_dir) = CALIBRATE_CONFIG
    eki = ClimaCalibrate.calibrate(
        backend,
        ekp,
        coupler_model_interface,
        n_iterations,
        PRIORS,
        output_dir,
        model_interface_filepath;
        experiment_dir = ClimaCalibrate.project_dir(),
    )
end
