# PHASE-4: phase3's parameters and observables, extended in time.
#
# Each member simulates a full OND season (Oct 1 to Jan 1, ~122 days from
# the Sep 1 son_ic initial condition, 30 day spinup), and the minibatch
# rotates the target year across iterations: 2010, 2008, 2007, 2006, then
# repeats. Three monthly means per variable enter each sample, so one
# iteration sees 3x the constraints of phase3 AND a seasonal-cycle signal
# (Oct -> Dec trend) that no single-month run contains. The SON spin-up
# study validated the Sep 1 start: the 30 day and 7 day spinup arms
# converged to consistent October states.
#
# 2009 is excluded from samples and covariance: the CALIPSO/CloudSat
# monthly record has no December 2009 (gap from 2009-11-30 to 2010-01-17).
# The covariance therefore uses 4 OND realizations (rank <= 3), one fewer
# than phase3's 5 Octobers.
#
# Cost: ~3.3x phase3 per iteration (122 vs 38 day members). 13 workers on
# 4 packed GPU nodes run one batch per iteration (~4 h); a 12 h job fits
# 2-3 iterations, so the run needs 1-2 resubmissions. Resume only at
# iteration boundaries and clean partial members first (restart bug).
#
# Predictions (falsifiable):
# 1. No COLLAPSED flag; residuals reach their floors as in phase3.
# 2. The rotating target does not destabilize the update: parameter drift
#    stays below 1x spread/iter from iteration 3 on, despite each
#    iteration fitting a different year.
# 3. q_liq lands in the 2.2-2.5e-4 cluster (zonal + cl family). A
#    seasonally-robust q_liq is the point of this run.
# 4. More data identifies rain_tau better: final spread below phase2's
#    0.78x prior.
# 5. December does not break the fit: per-month residuals of the final
#    iteration show no month with RMS above 2 sigma (a Dec blowup means
#    the closures carry a seasonal bias the single-month runs never saw).
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_3mo.yml",
)

# OND of one year per iteration, rotating over the CALIPSO-complete years.
const SAMPLE_YEARS = [2010, 2008, 2007, 2006, 2010, 2008]

sample_date_ranges =
    [(Dates.DateTime(y, 10, 1), Dates.DateTime(y, 12, 1)) for y in SAMPLE_YEARS]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 12, 1)) for
    y in (2006, 2007, 2008, 2010)
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_phase4_ond")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "pr", "swcre", "lwcre", "cl"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(30),
    output_dir,
    rng_seed = 42,
)

# Floors carried over from phase3 (validated per observable).
const OBS_NOISE_GROUPS = [
    (short_names = ["lwp", "pr"], model_error_scale = 0.2),
    (short_names = ["swcre"], model_error_scale = 0.5),
    (short_names = ["lwcre"], model_error_scale = 0.25),
    (short_names = ["cl"], model_error_scale = 0.5),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# Same priors as phase3: the time extension should re-derive the answer,
# not inherit it.
const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        3e-4, 1.5e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1450, 500, 300, 3600),
    PD.constrained_gaussian("detr_massflux_vertdiv_coeff", 0.35, 0.15, 0.0, 0.8),
    PD.constrained_gaussian("detr_buoy_coeff", 1.0, 0.6, 0.05, 5.0),
    PD.constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    PD.constrained_gaussian("cloud_fraction_steepness_scale", 0.7, 0.2, 0.1, 2.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
