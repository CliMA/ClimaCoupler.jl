#=
Shared calibration configuration - priors and ensemble settings.
This is the SINGLE SOURCE OF TRUTH for calibration parameters.

Included by both run_calibration.jl and precompute_ekp_inputs.jl
=#

import EnsembleKalmanProcesses.ParameterDistributions as PD
import EnsembleKalmanProcesses as EKP

# ==========================================================================
# PRIORS - Define your calibration parameters here
# ==========================================================================
const CALIBRATION_PRIORS = [
    # Atmospheric parameters
    PD.constrained_gaussian("entr_inv_tau", 0.002, 0.0015, 0.0, 0.01),
    # PD.constrained_gaussian("precipitation_timescale", 600, 300, 100, 1000),
    
    # Land parameters  
    PD.constrained_gaussian("leaf_Cd", 0.01, 0.005, 0.0, 0.1),
]

const CALIBRATION_PRIOR = EKP.combine_distributions(CALIBRATION_PRIORS)

# ==========================================================================
# ENSEMBLE SETTINGS
# ==========================================================================
# For TransformInversion/Inversion: set ensemble_size freely (typically 5-20)
# For TransformUnscented: this is IGNORED (uses 2*n_params + 1 automatically)
const CALIBRATION_ENSEMBLE_SIZE = 5

# Random seed for reproducibility
const CALIBRATION_RNG_SEED = 42

# Noise scalar for observation covariance (must match generate_observations.jl)
# This represents the noise level relative to normalized data variance
const CALIBRATION_NOISE_SCALAR = 0.5
