# Define which coupler file to use
config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Calibration target: a single fixed observation (Oct 2010) used every iteration.
# This is a standard fixed-target EKI setup and keeps the error trajectory clean
# and comparable across iterations. The list is repeated per minibatch slot.
sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:6]

# Dates used ONLY to estimate the observational + internal-variability covariance
# (see generate_observations.jl). Kept separate from the calibration target so
# that broadening the covariance sample does not change what is being calibrated.
# The same calendar month across several years gives the interannual spread the
# SVDplusD covariance needs. Requirements: all these dates must exist in the
# observational time series, and the calibration target date(s) above must be
# included here (SVDplusD requires the sampled date to be one of these).
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2005:2010
]

# On Derecho, it is preferable to save the calibration output to the scratch
# directory (e.g. "/glade/derecho/scratch")
output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    # Note: Pressure-level variables require model output with
    # pressure_coordinates: true in config
    short_names = ["lwp", "ta", "hur"],
    minibatch_size = 1,
    n_iterations = 5,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Used in generate_observations.jl and observation_map.jl
# Units: Pa (not hPa)
const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]

# Normalization is DISABLED: we now use a data-informed (SVDplusD) covariance
# that already carries each variable's physical scale, so per-variable
# normalization is both unnecessary and unsupported with that covariance (see
# generate_observations.jl and preprocess_sim_vars in observation_map.jl, which
# both skip normalization when this file is absent). Delete any stale
# normalization_stats.jld2 in output_dir before rerunning.
const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# Parameters chosen to actually be informed by the observation vector
# (lwp, ta, hur). The previous single knob, rain_autoconversion_timescale, only
# affects lwp/precip and left the ta/hur misfit (6 of 7 fields) untouched, so the
# error could not drop. We keep it for lwp and add mixing-length parameters that
# control turbulent heat/moisture transport and therefore the ta/hur profiles.
#
# NOTE ON COST: TransformUnscented uses 2p+1 members, so this 3-parameter set
# runs 7 forward models per iteration (vs 3 before). Add/remove parameters below
# to trade cost for the number of degrees of freedom being constrained.
const CALIBRATION_PRIORS = [
    # lwp / warm-rain
    PD.constrained_gaussian("rain_autoconversion_timescale", 1800, 300, 300, 3600),
    # ta / hur via turbulent mixing
    PD.constrained_gaussian("mixing_length_eddy_viscosity_coefficient", 0.2, 0.1, 0, 1.0),
    PD.constrained_gaussian("mixing_length_diss_coeff", 0.22, 0.15, 0.0, 10.0),
    # Additional candidates (uncomment to constrain more DOF, at 2 more members each):
    # PD.constrained_gaussian("mixing_length_tke_surf_flux_coeff", 8.0, 4.0, 0, 100.0),
    # PD.constrained_gaussian("Tq_correlation_coefficient", 0.4, 0.4, -1.0, 1.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
