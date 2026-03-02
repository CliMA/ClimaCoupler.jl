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
    # Gravity wave parameter (primary)
    PD.constrained_gaussian("nogw_Bt_0", 0.0043, 0.002, 0.001, 0.01),

    # Optional: uncomment to also calibrate OGW parameter
    # PD.constrained_gaussian("ogw_mountain_height_width_exponent", 0.4, 0.2, 0.0, 1.0),
]

const CALIBRATION_PRIOR = EKP.combine_distributions(CALIBRATION_PRIORS)

# ==========================================================================
# ENSEMBLE SETTINGS
# ==========================================================================
# For TransformInversion/Inversion: set ensemble_size freely (typically 5-20)
# For TransformUnscented: this is IGNORED (uses 2*n_params + 1 automatically)
const CALIBRATION_ENSEMBLE_SIZE = 7

# Random seed for reproducibility
const CALIBRATION_RNG_SEED = 42

# Noise scalar for observation covariance
# Used by both generate_observations.jl and precompute_ekp_inputs.jl
# For normalized data (unit variance): 0.5 = 50% of std, 2.0 = 200% of std
const CALIBRATION_NOISE_SCALAR = 0.5
