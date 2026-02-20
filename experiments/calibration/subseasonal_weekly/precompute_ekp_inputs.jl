#=
Precompute heavy EKP inputs for TransformInversion/Inversion processes.
Run this on a compute node (e.g., cpudev) before running run_calibration.jl from tmux.

NOTE: TransformUnscented does NOT need this - it works directly from tmux!
      This precompute is ONLY needed for TransformInversion/Inversion.

This script:
1. Loads the observation vector
2. Extracts observation data (y) and noise covariance (Γ)  
3. Constructs the initial ensemble
4. Saves everything to ekp_inputs.jld2
=#

using Dates
import Random
import ClimaCoupler
import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions as PD
import JLD2

# Override JLD2's default_iotype to use IOStream instead of MmapIO
JLD2.default_iotype() = IOStream

# Load shared priors from single source of truth
include(joinpath(@__DIR__, "calibration_priors.jl"))
prior = CALIBRATION_PRIOR
ensemble_size = CALIBRATION_ENSEMBLE_SIZE
rng_seed = CALIBRATION_RNG_SEED

@info "Loading observation vector..."
observation_vector = JLD2.load_object(
    joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal_weekly/obs_vec.jld2"),
)
@info "Observation vector loaded: $(length(observation_vector)) samples"

@info "Extracting observation data..."
obs = first(observation_vector)
y = EKP.get_obs(obs)
@info "Observation vector size: $(length(y))"

# NOTE: We store noise_scalar instead of full covariance matrix (saves ~5GB!)
# The noise covariance will be reconstructed as noise_scalar * I (UniformScaling)
noise_scalar = CALIBRATION_NOISE_SCALAR
@info "Noise scalar: $noise_scalar (will use UniformScaling in run_calibration.jl)"

@info "Constructing initial ensemble (ensemble_size=$ensemble_size)..."
rng = Random.MersenneTwister(rng_seed)
initial_ensemble = EKP.construct_initial_ensemble(rng, prior, ensemble_size)
@info "Initial ensemble shape: $(size(initial_ensemble))"

# Save to file (compact - no huge covariance matrix!)
output_path = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal_weekly/ekp_inputs.jld2")
@info "Saving to $output_path..."
JLD2.jldsave(output_path; 
    y = y,
    noise_scalar = noise_scalar,
    initial_ensemble = initial_ensemble,
    ensemble_size = ensemble_size,
    rng_seed = rng_seed,
)

@info "Done! EKP inputs saved to ekp_inputs.jld2"
@info "You can now run run_calibration.jl from tmux"
