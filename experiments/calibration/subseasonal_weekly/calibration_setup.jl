import EnsembleKalmanProcesses.ParameterDistributions as PD
import EnsembleKalmanProcesses as EKP


const CALIBRATION_PRIORS = [
    # Atmospheric parameters
    PD.constrained_gaussian("entr_inv_tau", 0.002, 0.001, 0.0, 0.01),

    # PiGroup linear regression coefficients
    # PD.ParameterDistribution(
    #     PD.VectorOfParameterized([PD.Normal(0.0, 5.0), PD.Normal(0.0, 5.0), PD.Normal(0.4, 0.2)]),
    #     repeat([PD.no_constraint()], 3),
    #     "entr_param_vec",
    # ),

    # PD.constrained_gaussian("detr_buoy_coeff", 0.12, 0.06, 0.0, 1.0),
    PD.constrained_gaussian("detr_vertdiv_coeff", 0.8, 0.25, 0.0, 5.0),
    # PD.constrained_gaussian("EDMF_surface_area", 0.1, 0.03, 0, 1),

    PD.constrained_gaussian("precipitation_timescale", 1200, 300, 300, 2400),
    PD.constrained_gaussian("diagnostic_covariance_coeff", 2.1, 0.5, 0.0, 10.0),
    PD.constrained_gaussian("Tq_correlation_coefficient", 0.4, 0.4, -1.0, 1.0),

    PD.constrained_gaussian("mixing_length_eddy_viscosity_coefficient", 0.2, 0.1, 0, 1.0),
    PD.constrained_gaussian("mixing_length_diss_coeff", 0.22, 0.15, 0.0, 10.0),
    PD.constrained_gaussian("mixing_length_tke_surf_flux_coeff", 8.0, 4.0, 0, 100.0),

    # Land parameters 
    # PD.constrained_gaussian("pmodel_cstar", 0.30, 0.15, 0.0, 1.0),
    # PD.constrained_gaussian("leaf_Cd", 0.01, 0.006, 0.0, 0.1),

    # gravity wave parameters
    # PD.constrained_gaussian("nogw_Bt_0", 0.0043, 0.003, 0.001, 0.01),
    # PD.constrained_gaussian("ogw_mountain_height_width_exponent", 0.4, 0.3, 0.0, 1.0),
]

const CALIBRATION_PRIOR = EKP.combine_distributions(CALIBRATION_PRIORS)

const CALIBRATION_ENSEMBLE_SIZE = 17  # For TransformUnscented: 2*n_params+1

const CALIBRATION_RNG_SEED = 42

# Noise scalar for observation covariance
# For normalized data (unit variance): this is the variance, not std
# 0.1 = 10% variance (std ≈ 0.32), 0.5 = 50% variance (std ≈ 0.71)
const CALIBRATION_NOISE_SCALAR = 5.0

# If true, normalize each variable to zero mean and unit variance
const NORMALIZE_VARIABLES = false


# ----- DATA SOURCE SETTINGS -----

# Variables to load from CERES instead of ERA5 (radiation variables)
# Set to empty vector [] to use ERA5 for all variables
const CERES_VARIABLES = ["rsut", "rlut", "rsutcs", "rlutcs", "rsds", "rsus", "rlds", "rlus", "swcre", "lwcre"]
# Note: CERES data is monthly

# ----- PRESSURE LEVEL SETTINGS -----

# Variables that require pressure-level selection from 3D data
# Format: "varname_XXXhPa" where XXX is the pressure level
const PRESSURE_LEVEL_VARIABLES = Dict{String, Vector{Float64}}(
    "ta" => [850.0, 500.0, 200.0],   # Temperature at 850, 500, 200 hPa
    "hur" => [850.0, 500.0, 200.0],
)

# Generate flat list of pressure-level variable names (e.g., "ta_850hPa")
function get_pressure_level_short_names()
    names = String[]
    for (var, levels) in PRESSURE_LEVEL_VARIABLES
        for level in levels
            push!(names, "$(var)_$(Int(level))hPa")
        end
    end
    return names
end
