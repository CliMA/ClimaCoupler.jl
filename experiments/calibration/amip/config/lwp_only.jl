# LWP-ONLY calibration config.
#
# Companion to pressure_levels.jl (lwp+cl). The noise investigation showed cl
# carries almost no usable signal for these parameters — its model–obs misfit
# (~0.083) is smaller than its own interannual noise (~0.10), because cloud
# fraction has enormous year-to-year variability. LWP, by contrast, has a
# misfit (~0.031) well above its noise (~0.018) and a strong physical knob (the
# warm-rain autoconversion threshold). This config isolates LWP so the inverse
# is driven only by the informative target.
#
# Select this config by setting ENV["CALIBRATION_CONFIG"] to this file's path
# when running generate_observations.jl / run_calibration.jl.

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Fixed calibration target: Oct 2010, repeated per minibatch slot. Length must be
# >= n_iterations + 1 (forward_model indexes sample_date_ranges[iter + 1]).
sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:11]

# Dates used ONLY to estimate the SVDplusD covariance (interannual spread). LWP
# (MAC) has coverage back further than CALIPSO cl, but we reuse the same 5
# Octobers known to work; the target (Oct 2010) must be included.
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

# Separate output_dir so this run's observations/output stay independent from the
# lwp+cl run. (Symlinked to scratch on Derecho — see the submit script.)
output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp"],
    minibatch_size = 1,
    n_iterations = 10,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Retained for the pressure/altitude preprocessing no-ops (lwp has neither).
const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

# Normalization disabled (SVDplusD carries the physical scale). Delete any stale
# normalization_stats.jld2 in output_dir before rerunning.
const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# LWP-only parameters: the warm-rain autoconversion threshold (primary LWP knob;
# currently 0 in the config, so clouds rain at the first drop of condensate →
# LWP ~3x low) and the rain autoconversion timescale (bounded to its physical
# range). cloud_fraction_eps_rel is dropped since there is no cl target here.
# TransformUnscented uses 2p+1 members, so 2 parameters -> 5 members.
const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        5e-4, 3e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1800, 700, 300, 3600),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
