using Dates
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

# Workers on Vast.ai meet the driver in a bucket, see the backend below
const ON_VAST = haskey(ENV, "CLIMA_BUCKET_ENDPOINT")
if ON_VAST
    using AWSS3
end

model_interface_filepath = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "amip",
    "model_interface.jl",
)
# The file picks the config, by TEST_CALIBRATION, and evaluates to the interface
coupler_model_interface = include(model_interface_filepath)

(; output_dir) = CALIBRATE_CONFIG
isdir(output_dir) || mkdir(output_dir)

if abspath(PROGRAM_FILE) == @__FILE__
    observation_vector_filepath = observation_vector_path(output_dir)
    isfile(observation_vector_filepath) || error(
        "No observations at $observation_vector_filepath. Generate them for this config with generate_observations.jl",
    )
    observation_vector = JLD2.load_object(observation_vector_filepath)

    (; sample_date_ranges) = CALIBRATE_CONFIG
    # One target per iteration, in order.
    obs_series = ClimaCalibrate.observation_series_from_samples(
        observation_vector,
        1,
        [string(Dates.year(start_date)) for (start_date, stop_date) in sample_date_ranges],
    )

    (; rng_seed, priors) = CALIBRATE_CONFIG
    rng = Random.MersenneTwister(rng_seed)

    ekp = EKP.EnsembleKalmanProcess(
        obs_series,
        EKP.TransformUnscented(priors, impose_prior = true);
        verbose = true,
        rng,
        scheduler = EKP.DataMisfitController(terminate_at = 1000000),
    )

    (; n_iterations) = CALIBRATE_CONFIG

    # On Vast.ai each member runs on a rented instance that pulls it from the
    # bucket. The smoke test runs the members serially in this process so that
    # CI does not submit Slurm jobs from inside its own allocation. A full
    # calibration puts one ensemble member on each GPU.
    backend = if ON_VAST
        ic_files = ClimaCalibrate.SharedFiles.BucketFiles(
            ENV["CLIMA_BUCKET_ENDPOINT"],
            ENV["CLIMA_BUCKET"],
            "era5_initial_conditions",
        )
        for start_date in unique(first.(sample_date_ranges) .- CALIBRATE_CONFIG.spinup)
            names = ClimaCoupler.Input.ERA5.output_filenames(start_date)
            all(name -> ispath(ic_files, name), names) && continue
            ic_dir = ClimaCoupler.Input.ERA5.fetch_initial_conditions(start_date)
            for name in names
                ispath(ic_files, name) ||
                    ClimaCalibrate.SharedFiles.upload(ic_files, name, joinpath(ic_dir, name))
            end
        end

        ClimaCalibrate.SharedFilesBackend(;
            lease = 900.0,
            files = ClimaCalibrate.SharedFiles.BucketFiles(
                ENV["CLIMA_BUCKET_ENDPOINT"],
                ENV["CLIMA_BUCKET"],
            ),
            outputs = [joinpath(get_job_id(CALIBRATE_CONFIG), "output_active")],
            workers = ClimaCalibrate.Backend.Vast.Workers(
                parse(Int, get(ENV, "VAST_WORKERS", string(EKP.get_N_ens(ekp))));
                gpu = get(ENV, "VAST_GPU", "RTX_3090"),
                disk = 60,
                image = "nefrathenrici/clima-base:pr-5-bf60034",
                env = Dict(
                    "CLIMACOMMS_CONTEXT" => "SINGLETON",
                    "CLIMACOMMS_DEVICE" => TEST_CALIBRATION ? "CPU" : "CUDA",
                    (
                        TEST_CALIBRATION ? ("TEST_CALIBRATION" => "1",) :
                        ("CLIMA_ARTIFACTS" => "historical_sst_sic",)
                    )...,
                ),
            ),
        )
    elseif TEST_CALIBRATION
        ClimaCalibrate.JuliaBackend()
    else
        ClimaCalibrate.add_workers(
            EKP.get_N_ens(ekp);
            device = :gpu,
            time = CALIBRATE_CONFIG.worker_walltime,
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
