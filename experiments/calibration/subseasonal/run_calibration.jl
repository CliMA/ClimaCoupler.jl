using Dates
using Distributed
using LinearAlgebra: I
import Random
import ClimaCalibrate
import ClimaAnalysis
import ClimaComms
import ClimaCoupler
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

# Override JLD2's default_iotype to use IOStream instead of MmapIO
# This avoids Bus errors from memory-mapped files on Lustre filesystem
JLD2.default_iotype() = IOStream


include(joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "api.jl"))
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

# Sample date ranges (start_date, end_date) for calibration
# Uncomment ONE of the options below:
sample_date_ranges = [
    # === OPTION 1: Single day (fast testing, ~30 min/iteration) ===
    # (DateTime(2013, 3, 18), DateTime(2013, 3, 18)),
    
    # === OPTION 2: 7-day period (production, ~3 hours/iteration) ===
    (DateTime(2023, 1, 15), DateTime(2023, 1, 21)),
    # (DateTime(2013, 3, 18), DateTime(2013, 3, 24)),
    
    # Add more ranges as needed:
    # (DateTime(2013, 3, 25), DateTime(2013, 3, 31)),
]

# Directory containing ERA5 weekly observation files
const ERA5_OBS_DIR = "/glade/campaign/univ/ucit0011/cchristo/wxquest_data/daily_weekly_stats/weekly"

const CALIBRATE_CONFIG = CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config/subseasonal_configs/wxquest_diagedmf.yml",
    ),
    # short_names = ["tas"],  # Start with tas only
    short_names = ["tas", "mslp"],# "pr"], # ADD TOA
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Day(1),  # Add 1 day so simulation covers full 7-day diagnostic period
    spinup = Dates.Day(0),
    # Use scratch filesystem - more reliable for JLD2/HDF5 on Lustre
    output_dir = "/glade/derecho/scratch/cchristo/calibration/exp27",  # Full gridpoint calibration
    obs_dir = ERA5_OBS_DIR,
    rng_seed = 42,
)

if abspath(PROGRAM_FILE) == @__FILE__
    # Load priors from shared config (single source of truth)
    include(joinpath(@__DIR__, "calibration_priors.jl"))
    prior = CALIBRATION_PRIOR
    # Ensemble size:
    # - TransformUnscented/Unscented: automatically 2*n_params + 1 (determined by process)
    # - TransformInversion/Inversion: uses CALIBRATION_ENSEMBLE_SIZE from calibration_priors.jl
    ensemble_size = CALIBRATION_ENSEMBLE_SIZE

    # Add PBS workers with GPU resources on Derecho
    if nworkers() == 1
        @info "Adding PBS workers with GPU resources..."
        addprocs(
            ClimaCalibrate.PBSManager(ensemble_size);
            q = "main",
            A = "UCIT0011",
            l_select = "1:ncpus=12:ngpus=1",
            l_walltime = "12:00:00",
        )
    end

    # Load api.jl on all workers first (defines CalibrateConfig type)
    api_file = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "api.jl")
    @everywhere include($api_file)
    
    # Share CALIBRATE_CONFIG with workers and load model interface
    @everywhere const CALIBRATE_CONFIG = $CALIBRATE_CONFIG
    @everywhere include($model_interface)

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

    # ==========================================================================
    # OPTION 1: TransformUnscented (efficient, quick convergence) [KNOWN WORKING]
    # - Uses ObservationSeries constructor
    # - Ensemble size automatically 2*n_params + 1
    # ==========================================================================
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
    # ekp_inputs_path = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/ekp_inputs.jld2")
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

    # Use WorkerBackend with PBS workers
    eki = ClimaCalibrate.calibrate(
        ClimaCalibrate.WorkerBackend(),
        ekp,
        CALIBRATE_CONFIG.n_iterations,
        prior,
        CALIBRATE_CONFIG.output_dir;
    )
end
