# LWP + precipitation (pr) config — first parameter-space expansion beyond LWP.
#
# Rationale: the warm-rain parameters q_liq_threshold and rain_autoconversion_
# timescale BOTH raise LWP by suppressing autoconversion, so LWP alone cannot
# separate them (rain_tau stayed loosely constrained ~1450-1670 in the LWP-only
# runs). Precipitation is the autoconversion *sink* — it is expected to respond
# to rain_tau differently than LWP (the reservoir) does, which is what would
# break the q_liq/rain_tau degeneracy. Before committing this to a calibration
# we sweep both parameters against {lwp, pr} (run_parameter_sweep.jl) to confirm
# pr carries signal above the internal-variability noise floor AND that its
# response is distinct from LWP's.
#
# Same clean single-year Oct-2010 platform as lwp_only.jl (7-day spin-up), so the
# only change from the trusted LWP baseline is the added pr target.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Fixed calibration target: Oct 2010, repeated per minibatch slot. Length must be
# >= n_iterations + 1.
sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:11]

# Dates used ONLY to estimate the SVDplusD covariance (interannual spread). The
# target (Oct 2010) must be included. GPCP pr and MAC lwp both cover these years.
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

# Separate output_dir so this run's observations/output stay independent.
# (Symlinked to scratch on Derecho — see the submit script.)
output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_pr")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "pr"],
    minibatch_size = 1,
    n_iterations = 10,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Retained for the pressure/altitude preprocessing no-ops (lwp/pr have neither).
const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

# Normalization disabled (SVDplusD carries the physical scale). Delete any stale
# normalization_stats.jld2 in output_dir before rerunning.
const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# Same two warm-rain parameters as lwp_only.jl. The point of adding pr is to
# constrain rain_autoconversion_timescale, which LWP alone leaves loose.
const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        5e-4, 3e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1800, 700, 300, 3600),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
