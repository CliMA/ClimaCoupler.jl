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
    # PD.constrained_gaussian("entr_inv_tau", 0.002, 0.0015, 0.0, 0.01),

    # Linear regression coefficients (6-element vector)
    # PD.VectorOfParameterized([PD.Normal(0.0, 5.0), PD.Normal(0.0, 5.0), PD.Normal(0.0, 5.0), PD.Normal(0.0, 5.0), PD.Normal(0.0, 5.0), PD.Normal(0.4, 0.2)]),
    PD.ParameterDistribution(
        PD.VectorOfParameterized([PD.Normal(0.0, 5.0), PD.Normal(0.0, 5.0), PD.Normal(0.4, 0.2)]),
        repeat([PD.no_constraint()], 3),
        "entr_param_vec",
    ),

    # PD.constrained_gaussian("detr_buoy_coeff", 0.12, 0.06, 0.0, 1.0),
    PD.constrained_gaussian("detr_vertdiv_coeff", 0.6, 0.25, 0.0, 5.0),

    # PD.constrained_gaussian("precipitation_timescale", 600, 300, 100, 2000),
    PD.constrained_gaussian("precipitation_timescale", 1200, 300, 300, 2400),
    PD.constrained_gaussian("diagnostic_covariance_coeff", 2.1, 0.5, 0, 10),
    PD.constrained_gaussian("Tq_correlation_coefficient", 0, 0.5, -1, 1),

    
    # PD.constrained_gaussian("mixing_length_eddy_viscosity_coefficient", 0.2, 0.1, 0, 1.0),
    PD.constrained_gaussian("mixing_length_diss_coeff", 4.2, 2.1, 0, 10.0),
    # PD.constrained_gaussian("mixing_length_tke_surf_flux_coeff", 8.0, 4.0, 0, 100.0),

    # PD.constrained_gaussian("EDMF_surface_area", 0.1, 0.03, 0, 1),
    # Land parameters 
    # PD.constrained_gaussian("pmodel_cstar", 0.30, 0.15, 0.0, 1.0),
    # PD.constrained_gaussian("leaf_Cd", 0.01, 0.006, 0.0, 0.1),

    # sea ice parameters
    # PD.constrained_gaussian("ice_albedo", 0.7, 0.05, 0.4, 0.9), # not supported 

    # gravity wave parameters
    # PD.constrained_gaussian("nogw_Bt_0", 0.0043, 0.003, 0.001, 0.01),
    # PD.constrained_gaussian("ogw_mountain_height_width_exponent", 0.4, 0.3, 0.0, 1.0),
]

const CALIBRATION_PRIOR = EKP.combine_distributions(CALIBRATION_PRIORS)

# ==========================================================================
# ENSEMBLE SETTINGS
# ==========================================================================
# For TransformInversion/Inversion: set ensemble_size freely (typically 5-20)
# For TransformUnscented: this is IGNORED (uses 2*n_params + 1 automatically)
const CALIBRATION_ENSEMBLE_SIZE = 17  # For TransformUnscented: 2*n_params+1 = 2*8+1 = 17

# Random seed for reproducibility
const CALIBRATION_RNG_SEED = 42

# Noise scalar for observation covariance
# Used by both generate_observations.jl and precompute_ekp_inputs.jl
# For normalized data (unit variance): this is the VARIANCE, not std
# 0.1 = 10% variance (std ≈ 0.32), 0.5 = 50% variance (std ≈ 0.71)
# Lower values make EKP more aggressive in fitting observations
const CALIBRATION_NOISE_SCALAR = 3.0

