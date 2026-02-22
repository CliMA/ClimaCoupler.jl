using Dates
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

# Override JLD2's default_iotype to use IOStream instead of MmapIO
# This avoids Bus errors from memory-mapped files on Lustre filesystem
JLD2.default_iotype() = IOStream

include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments/calibration/subseasonal_weekly/observation_map.jl",
    ),
)
model_interface = joinpath(
    pkgdir(ClimaCoupler),
    "experiments",
    "calibration",
    "subseasonal_weekly",
    "model_interface.jl",
)

# Sample date range for calibration (repeated for each iteration)
const BASE_DATE_RANGE = (DateTime(2010, 1, 1), DateTime(2010, 1, 7))
const N_ITERATIONS = 6

# Repeat the date range for each iteration so we can reuse subseasonal's forward_model
# which indexes by sample_date_ranges[iter + 1]
sample_date_ranges = fill(BASE_DATE_RANGE, N_ITERATIONS)

# Directory containing ERA5 weekly observation files
const ERA5_OBS_DIR = "/glade/campaign/univ/ucit0011/cchristo/wxquest_data/daily_weekly_stats/weekly"

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config/subseasonal_configs/wxquest_diagedmf_monthly_calibration.yml",
    ),
    short_names = ["tas", "mslp", "rsut", "rlut"],
    minibatch_size = 1,
    n_iterations = N_ITERATIONS,
    sample_date_ranges,
    extend = Dates.Day(1),
    spinup = Dates.Day(0),
    output_dir = "/glade/derecho/scratch/cchristo/calibration/exp28",
    rng_seed = 42,
)

if abspath(PROGRAM_FILE) == @__FILE__
    # Load priors from shared config (single source of truth)
    include(joinpath(@__DIR__, "calibration_priors.jl"))
    prior = CALIBRATION_PRIOR

    observation_vector = JLD2.load_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal_weekly/obs_vec.jld2"),
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
        scheduler = EKP.DataMisfitController(terminate_at = 1000000),
    )

    # # ==========================================================================
    # # OPTION 2: TransformInversion (robust to failures, flexible ensemble size)
    # # - Loads pre-computed inputs from ekp_inputs.jld2 (run precompute_ekp_inputs.jl first!)
    # # - To use: comment out OPTION 1 above, uncomment below
    # # ==========================================================================
    # ekp_inputs_path = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal_weekly/ekp_inputs.jld2")
    # if !isfile(ekp_inputs_path)
    #     error("ekp_inputs.jld2 not found! Run precompute_ekp_inputs.jl first (use: qsub precompute.pbs)")
    # end
    # @info "Loading pre-computed EKP inputs from $ekp_inputs_path"
    # ekp_inputs = JLD2.load(ekp_inputs_path)
    # y = ekp_inputs["y"]
    # noise_scalar = ekp_inputs["noise_scalar"]
    # initial_ensemble = ekp_inputs["initial_ensemble"]
    # @info "Loaded: y=$(length(y)) points, ensemble=$(size(initial_ensemble)), noise=$noise_scalar"
    
    # # Use UniformScaling for noise covariance (efficient - no huge matrix!)
    # Γ = noise_scalar * I
    
    # ekp = EKP.EnsembleKalmanProcess(
    #     initial_ensemble,
    #     y,
    #     Γ,
    #     EKP.TransformInversion();
    #     verbose = true,
    #     rng,
    #     scheduler = EKP.DataMisfitController(terminate_at = 10000),
    # )

    # ==========================================================================
    # OPTION 3: Basic Inversion (simple, flexible)
    # - Same constructor as OPTION 2, just different process
    # - To use: uncomment obs/y/Γ/initial_ensemble from OPTION 2, then uncomment below
    # ==========================================================================
    # ekp = EKP.EnsembleKalmanProcess(
    #     initial_ensemble,
    #     y,
    #     Γ,
    #     EKP.Inversion();
    #     verbose = true,
    #     rng,
    #     scheduler = EKP.DataMisfitController(terminate_at = 1000),
    # )
    # backend = ClimaCalibrate.ClimaGPUBackend(;
    #     hpc_kwargs = ClimaCalibrate.kwargs(gpus = 1, time = 60 * 12),
    #     model_interface,
    #     verbose = true,
    # )
    backend = ClimaCalibrate.WorkerBackend()

    eki = ClimaCalibrate.calibrate(
        backend,
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir;
    )
end
