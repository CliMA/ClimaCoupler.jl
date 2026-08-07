# Terminal-velocity calibration: same observations as
# lwp_clt_swcre_5d_relax (subseasonal, ERA5 init 2010-09-24, 7-day spinup,
# October 2010 monthly means, 5 degree grid coarsened 2x), but a different
# parameter set. Five of the relax parameters are frozen at its posterior
# in toml/calibration_relax_posterior.toml (snow_tau 1360, cond_evap 100,
# margin 1, eps_rel 0.05, steepness 0.55) and the three fixed terminal
# velocities are opened instead. 6 parameters -> 13 members, so this run
# costs about a quarter less than relax per iteration.
#
# The velocities are live only because ClimaAtmos defaults
# fixed_terminal_velocity to true (0.42.2 default_config.yml); they set the
# constant fall speeds in FixedTerminalVelocity. ClimaParams baselines are
# ice 0.01, rain 5.0, snow 1.0 m/s. The ice prior below is centered at
# 0.1 m/s, near observed crystal fall speeds and 1.8 sigma above that
# baseline, so it deliberately asserts the baseline is too slow. Cloud
# liquid keeps its 0.01 m/s default and is not calibrated.
#
# Predictions (falsifiable):
# 1. Rain fall speed is identifiable from lwp and swcre (it sets how fast
#    condensate leaves the column): spread contracts below 0.5x prior and
#    the posterior lands in 3.5-6.5 m/s. Staying near 1x prior means
#    column-integrated observables cannot see it at this resolution.
# 2. Ice fall speed is NOT identifiable: spread stays above 0.5x prior.
#    clt at 5 degrees carried little ice information in the relax pair
#    (sublimation-deposition drifted 278 -> 491 s across runs with wide
#    spread). Contraction here would be the first ice constraint from
#    these three observables.
# 3. Freezing the five cloud-fraction and phase parameters does not cost
#    fit: final residuals land within 0.05 sigma of relax (lwp 0.67,
#    clt 0.59, swcre 0.68). A higher clt residual means eps_rel, steepness
#    or the margin were doing work the velocities cannot replace.
# 4. q_liq is robust to the parameter set: posterior in 1.3-2.0e-4, which
#    spans relax (1.93e-4, broken wide prior) and icabl_sub (1.34e-4,
#    checked prior, Sep+Oct). Outside that window means q_liq was
#    compensating for the now-frozen cloud-fraction shape.
# 5. q_liq and rain fall speed both control rain removal, so they are the
#    degeneracy risk: |correlation| stays below 0.9 in the final ensemble.
#    Above it, the pair is not jointly identifiable from lwp alone.
#
# Verdict (run 2026-08-04, jobs 7003518 + workers, 6 iterations):
# Posterior: q_liq 1.493e-4, rain_tau 1512, subl_dep 608, v_ice 0.0931,
# v_rain 5.421, v_snow 1.689. Residuals lwp 0.647, clt 0.654, swcre 0.611.
# 1. PASS. Rain fall speed 5.42 m/s, spread 0.055x prior.
# 2. FALSIFIED as stated, upheld in substance. Ice spread contracted to
#    0.046x prior, so the stated metric fails. But spread contracted to
#    0.042-0.068x for EVERY parameter (relax: 0.03-0.04x for all eight),
#    so the spread ratio tracks the scheduler's accumulated pseudo-time,
#    not per-parameter information. Mean displacement in prior-sigma units
#    discriminates: v_snow +1.45, rain_tau -0.71, q_liq -0.46, v_rain
#    +0.42, subl_dep +0.39, v_ice -0.06. Ice is the only parameter the
#    data left alone, so lwp/clt/swcre carry no ice-velocity information.
#    USE MEAN DISPLACEMENT, NOT SPREAD RATIO, as the identifiability test.
# 3. PARTIAL. lwp 0.647 (within 0.05 of relax's 0.67) and swcre 0.611
#    (0.069 better), but clt 0.654 is 0.064 WORSE, and clt is the only
#    observable that degrades across iterations (0.620 -> 0.654). Total
#    fit still beats relax (loss minimum 0.181 vs 0.215 on identical
#    observations with 6 parameters instead of 8). Freezing eps_rel,
#    steepness and the margin costs clt specifically, which is direct
#    evidence that those three were serving clt in relax.
# 4. PASS. q_liq 1.493e-4, near icabl_sub's 1.34e-4.
# 5. PASS decisively. No parameter pair exceeds |r| = 0.5 in the final
#    ensemble, so the fall speeds do not trade against autoconversion.
# Loss (EKP.get_error): 0.0388, 0.0362, 0.0302, 0.0323, 0.0391, 0.0386.
# The minimum is at iteration 3 and the loss then rises. relax shows the
# same shape on the same observations, so it is a property of the fixed
# Delta t = 0.1 scheduler against this over-informative observation, not
# of the parameter set.
# CAVEAT: the ice prior asserted 0.1 m/s against a ClimaParams baseline of
# 0.01. Since the data did not move it, this run gives NO evidence for
# changing that baseline; the 0.093 posterior is the prior talking.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_vterm.yml",
)

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir =
    joinpath(pkgdir(ClimaCoupler), "amip_calibration_lwp_clt_swcre_5d_vterm")

const COARSEN_FACTOR = 2

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "clt", "swcre"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# lwp/swcre keep the 1.5x-scaled 10 degree floors; clt drops to 0.35
# (1.4x the 10d companion's 0.25) to make it informative.
const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [
    (
        short_names = ["lwp"],
        model_error_scale = 0.6,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["clt"],
        model_error_scale = 0.35,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["swcre"],
        model_error_scale = 0.75,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# checked_constrained_gaussian: EKP's constrained_gaussian silently ignores
# targets with magnitudes below ~1e-2 (see prior_tools.jl). Every earlier
# q_liq prior was in fact a unit normal over its bounds (mean = bounds
# center 5e-4), which is why iteration-1 q_liq always started there.
include(joinpath(@__DIR__, "..", "prior_tools.jl"))

const CALIBRATION_PRIORS = [
    checked_constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        # Tightened around the reproduced 2D answer (1.6-1.9e-4 across four
        # runs). The zonal cluster (2.2-2.55e-4) stays within ~1 sigma.
        1.8e-4, 0.6e-4, 0.0, 1e-3,
    ),
    checked_constrained_gaussian("rain_autoconversion_timescale", 1800, 400, 300, 3600),
    checked_constrained_gaussian("sublimation_deposition_timescale", 500, 400, 20, 3600),
    checked_constrained_gaussian("fixed_cloud_ice_terminal_velocity", 0.1, 0.05, 0, 1),
    checked_constrained_gaussian("fixed_rain_terminal_velocity", 5, 1, 0, 10),
    checked_constrained_gaussian("fixed_snow_terminal_velocity", 1, 0.5, 0, 5),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
