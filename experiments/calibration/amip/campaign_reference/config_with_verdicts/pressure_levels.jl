# Define which coupler file to use
config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Calibration target: a single fixed observation (Oct 2010) used every iteration.
# This is a standard fixed-target EKI setup and keeps the error trajectory clean
# and comparable across iterations. The list is repeated per minibatch slot.
# Length must be >= n_iterations + 1, because forward_model indexes
# sample_date_ranges[iter + 1] (iter runs 1:n_iterations).
sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:11]

# Dates used ONLY to estimate the observational + internal-variability covariance
# (see generate_observations.jl). Kept separate from the calibration target so
# that broadening the covariance sample does not change what is being calibrated.
# The same calendar month across several years gives the interannual spread the
# SVDplusD covariance needs. Requirements: all these dates must exist in the
# observational time series, and the calibration target date(s) above must be
# included here (SVDplusD requires the sampled date to be one of these).
# Note: CALIPSO/CloudSat cloud-fraction coverage starts 2006-06, so the earliest
# usable October is 2006 (the binding constraint across lwp+cl).
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

# On Derecho, it is preferable to save the calibration output to the scratch
# directory (e.g. "/glade/derecho/scratch")
output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    # Calibrating liquid water path (MAC) + cloud fraction (CALIPSO/CloudSat `cl`,
    # level-resolved). `cl` is compared on altitude levels (see ALTITUDE_LEVELS).
    short_names = ["lwp", "cl"],
    minibatch_size = 1,
    # With the fixed Δt = 0.1 scheduler in run_calibration.jl, accumulated
    # pseudo-time T reaches ~1 (full posterior) over ~10 iterations.
    n_iterations = 10,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Used in generate_observations.jl and observation_map.jl
# Units: Pa (not hPa). Retained for the pressure-level preprocessing no-op; not
# used now that the vector is lwp + cl.
const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]

# Altitude levels (meters) at which cloud fraction `cl` is compared to CALIPSO.
# low / mid / high. Both obs (CALIPSO 240 m grid) and sim (model z-levels) are
# put onto these via nearest-level selection (see select_altitude_levels).
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

# Normalization is DISABLED: we now use a data-informed (SVDplusD) covariance
# that already carries each variable's physical scale, so per-variable
# normalization is both unnecessary and unsupported with that covariance (see
# generate_observations.jl and preprocess_sim_vars in observation_map.jl, which
# both skip normalization when this file is absent). Delete any stale
# normalization_stats.jld2 in output_dir before rerunning.
const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# Parameters targeted at LWP + cloud fraction, based on the sweep findings:
#   * cloud_liquid_water_specific_humidity_autoconversion_threshold — the warm-rain
#     autoconversion threshold. Currently 0 in the config (clouds rain at the first
#     drop of condensate → LWP ~3x low); the sweeps showed ~5-6e-4 recovers MAC LWP
#     without changing precip. Primary LWP knob.
#   * rain_autoconversion_timescale — bounded to its physical range [300, 3600] s so
#     it can't run to unphysical values; shapes the LWP field alongside the threshold.
#   * cloud_fraction_eps_rel — the saturation-variability floor in the SGS quadrature
#     cloud-fraction closure. Primary cloud-fraction knob.
#
# NOTE ON COST: TransformUnscented uses 2p+1 members, so 3 parameters -> 7 members.
const CALIBRATION_PRIORS = [
    # LWP / warm-rain
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        5e-4, 3e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1800, 700, 300, 3600),
    # Cloud fraction
    PD.constrained_gaussian("cloud_fraction_eps_rel", 0.03, 0.02, 0.0, 0.15),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
