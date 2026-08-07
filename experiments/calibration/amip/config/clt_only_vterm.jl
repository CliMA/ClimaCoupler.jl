# clt-only diagnostic for lwp_clt_swcre_5d_vterm. Identical in every way
# except the observable set: same 6 parameters, same frozen 5 in
# toml/calibration_relax_posterior.toml, same subseasonal October 2010
# target, same covariance years, same floor for clt (0.35) and the same
# 800 km correlated floor, same seed and iteration count.
#
# Purpose. In the vterm run clt degraded monotonically (0.620 at iteration
# 1 to 0.654 at iteration 6) while lwp and swcre improved. Two
# explanations fit that: competition (the parameters can help clt but the
# joint fit prefers spending them on lwp/swcre) or structure (these 6
# parameters cannot move clt at all, because the parameters clt responds
# to - eps_rel, steepness, the floor-release margin - are frozen). Running
# clt alone separates them, and weighting experiments cannot.
#
# Predictions (falsifiable):
# 1. THE TEST. clt falls below 0.620, the best it reached in the joint run,
#    and keeps falling. That means competition, and the fix for the joint
#    run is reweighting or a longer schedule. If clt instead stalls near
#    0.62 or degrades again with nothing to compete against, the cause is
#    structural and clt needs eps_rel or steepness back, not a weight.
# 2. v_snow decomposes. It moved +1.45 prior-sigma (0.96 -> 1.69 m/s) in
#    the joint run. If clt alone also pushes it above ~1.5, the signal is
#    shared; if it stays near 1.0, the rise was bought by lwp and swcre.
# 3. q_liq loosens. Without lwp, condensate amount is only weakly
#    constrained, so expect q_liq above 1.49e-4 with a wider final spread.
#    Staying at 1.49e-4 would mean clt alone pins the autoconversion
#    threshold, which no earlier run has shown.
# 4. Information floor. At least one parameter should exceed 0.5 prior
#    sigma of mean displacement. If none do, clt at 5 degrees carries
#    almost no parameter information by itself, which would independently
#    explain why it loses every trade in the joint fit.
#
# Read the verdict with mean displacement in prior-sigma units, not the
# spread ratio: the vterm run showed the spread ratio contracts uniformly
# to 0.04-0.07x for every parameter regardless of information content.
#
# Verdict (run 2026-08-04, job 7009893, CANCELLED after 2 iterations: the
# answer was already unambiguous and the remaining 4 iterations would have
# cost ~65 GPU-hours to confirm a flat line):
# clt 0.620 at iteration 1 (identical to vterm, as it must be: same prior,
# seed and sigma points give identical members) and 0.623 at iteration 2.
# It got WORSE with nothing to compete against. Loss flat, 0.0204 ->
# 0.0206. Mean displacement after the first clt-only update, against the
# same step in vterm: v_snow +0.11 vs +1.23 sigma, v_rain +0.22 vs +0.49,
# rain_tau -0.13 vs -0.48, q_liq -0.03 vs -0.37.
# 1. FALSIFIED, informatively. clt never beats 0.620, so the joint run's
#    degradation was NOT competition. But the cause is not a missing
#    control either (see below), so the config-header framing of this
#    prediction was itself wrong.
# 2. ANSWERED. v_snow moves +0.11 sigma here vs +1.23 in the joint run, so
#    its +1.45 sigma posterior shift was bought entirely by lwp and swcre.
#    clt does not ask for faster snow.
# 3. and 4. UNRESOLVED by cancellation, though the direction is clear: the
#    largest displacement at iteration 2 was v_rain at +0.22 sigma, so
#    prediction 4's 0.5 sigma bar was unlikely to be cleared.
# ROOT CAUSE. clt is already fit to well inside its assumed noise. At
# iteration 1 the prior misfit is exactly zero, so bayes_loss = misfit / p
# gives a whitened misfit of 0.0204 * 6 = 0.122, i.e. an RMS residual of
# 0.35 of the assumed noise. There is nothing left to gain and the
# optimizer correctly does almost nothing. Note 0.35 whitened vs 0.62 on
# the covariance diagonal: the gap is the 800 km correlated floor
# absorbing the smooth part of clt's error, which is exactly what a
# large-scale cloud-cover deficit looks like.
# IMPLICATION for making clt informative, in preference order:
#   1. Shorten clt's decorrelation length, or make its floor diagonal. The
#      0.35-vs-0.62 gap localizes clt's lost weight in the correlated
#      part, and a near-global bias should not be nearly free.
#   2. Lower model_error_scale from 0.35 (misfit scales as 1/f^2, so a
#      misfit near 1 needs roughly 0.12). Bounded below by the low-rank
#      interannual part of Gamma, which is real climate variability.
# Reweighting alone is defensible here, contrary to what the joint run's
# monotone clt degradation first suggested.
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

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_clt_only_vterm")

const COARSEN_FACTOR = 2

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["clt"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# clt keeps the 0.35 floor of the joint run. Changing the floor and the
# observable set at once would confound the diagnostic.
const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [
    (
        short_names = ["clt"],
        model_error_scale = 0.35,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

include(joinpath(@__DIR__, "..", "prior_tools.jl"))

const CALIBRATION_PRIORS = [
    checked_constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        1.8e-4, 0.6e-4, 0.0, 1e-3,
    ),
    checked_constrained_gaussian("rain_autoconversion_timescale", 1800, 400, 300, 3600),
    checked_constrained_gaussian("sublimation_deposition_timescale", 500, 400, 20, 3600),
    checked_constrained_gaussian("fixed_cloud_ice_terminal_velocity", 0.1, 0.05, 0, 1),
    checked_constrained_gaussian("fixed_rain_terminal_velocity", 5, 1, 0, 10),
    checked_constrained_gaussian("fixed_snow_terminal_velocity", 1, 0.5, 0, 5),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
