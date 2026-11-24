using Dates
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

model_interface_filepath = joinpath(@__DIR__, "model_interface.jl")
include(model_interface_filepath)
include(joinpath(@__DIR__, "step_size_rule.jl"))

# Choose which calibration config to use
config_dir = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip", "config")
default_config_path = joinpath(config_dir, "year2010_rsut_rlut_pr.jl")

test_calibration_config_path = joinpath(config_dir, "pipeline_test.jl")
const TEST_CALIBRATION = haskey(ENV, "TEST_CALIBRATION")

config_path =
    TEST_CALIBRATION ? test_calibration_config_path :
    get(ENV, "CALIBRATION_CONFIG", default_config_path)

@info "Using calibration configuration in: $config_path"
include(config_path)

(; output_dir) = CALIBRATE_CONFIG
isdir(output_dir) || mkdir(output_dir)

if abspath(PROGRAM_FILE) == @__FILE__
    observation_vector_filepath = observation_vector_path(output_dir)
    isfile(observation_vector_filepath) || error(
        "No observations at $observation_vector_filepath. Generate them for this config with generate_observations.jl",
    )
    observation_vector = JLD2.load_object(observation_vector_filepath)

    (; sample_date_ranges, rng_seed, priors, observables) = CALIBRATE_CONFIG
    # One target per iteration, in order.
    obs_series = ClimaCalibrate.observation_series_from_samples(
        observation_vector,
        1,
        [string(Dates.year(start_date)) for (start_date, stop_date) in sample_date_ranges],
    )

    rng = Random.MersenneTwister(rng_seed)

    ekp = EKP.EnsembleKalmanProcess(
        obs_series,
        EKP.TransformUnscented(priors, impose_prior = true);
        verbose = true,
        rng,
        scheduler = all(o -> isnothing(o.effective_entries), observables) ?
                    EKP.DataMisfitController(terminate_at = 1000000) :
                    EffectiveEntriesController(
                        effective_entries(first(observation_vector), observables),
                    ),
    )

    coupler_model_interface =
        CouplerModelInterface(CALIBRATE_CONFIG; emulate_diagnostics = TEST_CALIBRATION)

    (; n_iterations) = CALIBRATE_CONFIG

    # The smoke test runs the members serially in this process so that CI does
    # not submit Slurm jobs from inside its own allocation. A full calibration
    # puts one ensemble member on each GPU.
    backend = if TEST_CALIBRATION
        ClimaCalibrate.JuliaBackend()
    else
        # Pre-fetch and interpolate any missing ERA5 initial conditions on the driver
        # before parallel workers start, so workers do not race on the same cache file.
        config_dict =
            ClimaCoupler.Input.get_coupler_config_dict(CALIBRATE_CONFIG.config_file)
        if get(config_dict, "initial_condition", nothing) == "WeatherModel"
            target_levels = collect(0.0:300.0:Float64(config_dict["z_max"]))
            era5_dir = get(config_dict, "era5_initial_condition_dir", nothing)
            for (start_date, _) in unique(sample_date_ranges)
                ClimaAtmos.weather_model_data_path(
                    start_date - CALIBRATE_CONFIG.spinup,
                    target_levels,
                    era5_dir,
                )
            end
        end
        n_ens = EKP.get_N_ens(ekp)
        ClimaCalibrate.add_workers(
            n_ens;
            device = :gpu,
            time = CALIBRATE_CONFIG.worker_walltime,
            max_relaunches = max(n_iterations - 1, 1) * n_ens,
            env = ["CLIMACOMMS_CONTEXT" => "SINGLETON", "CLIMACOMMS_DEVICE" => "CUDA"],
            mem = "32GB",
        )
        ClimaCalibrate.@worker_setup include($model_interface_filepath)
        ClimaCalibrate.WorkerBackend()
    end

    eki = ClimaCalibrate.calibrate(
        backend,
        ekp,
        coupler_model_interface,
        n_iterations,
        priors,
        output_dir,
    )
end
