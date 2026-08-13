# Multi-year September calibration: the two fixes that the vterm, clt-only
# and Southern Ocean band runs pointed to, applied together.
#
# FIX 1, sampling. Five DISTINCT Septembers (2006-2010) instead of seven
# copies of one month. Every earlier run in this campaign graded against a
# single month, and the vterm spread-contraction test showed that member
# spread at 38 days is mostly internal variability: signal-to-noise was
# lwp 1.25, clt 0.70, swcre 0.57, so two of three observables carried more
# weather than parameter response. Distinct years with minibatch_size = 1
# is EKP minibatching: each iteration targets a different year, so the
# parameters cannot be tuned to one year's storm track. The weather floor
# should fall roughly as sqrt(5) across the run.
#
# FIX 2, parameters. cloud_fraction_eps_rel and
# cloud_fraction_steepness_scale are FREE again. They are the only
# parameters that ever moved clt: relax improved clt with them free (ratio
# 2.7), while vterm and the band both degraded clt with them frozen (3.7
# and 5.7). fixed_cloud_ice_terminal_velocity is DROPPED; the data ignored
# it in four runs (-0.06 prior sigma globally, -0.17 in the band). q_liq
# stays free because fixing it in the band run cost lwp its largest lever
# (lwp spread fell 0.456 -> 0.145). 6 parameters, 13 members.
#
# WHY SEPTEMBER, not October. son_ic holds ERA5 initial conditions for
# 2006-2010 only at Aug 25 and Sep 1, and the default
# wxquest_initial_conditions artifact has late-September dates for 2010
# only. A 7-day spinup into October is therefore impossible for the earlier
# years. Starting Aug 25 and targeting September keeps vterm's protocol
# exactly (7-day spinup, 37-day members, 5 degree grid, same floors, 800 km
# correlated floor), so the leverage comparison is controlled apart from the
# month itself. The alternative, October from a Sep 1 start, would have
# meant a 30-day spinup and 61-day members: a different protocol and 1.6x
# the cost.
#
# Compare against these vterm iteration-1 baselines (leverage.jl):
#   lwp 0.696 / 0.456 = 1.5,  clt 0.620 / 0.166 = 3.7,  swcre 0.701 / 0.270 = 2.6
#
# Predictions (falsifiable):
# 1. THE GO/NO-GO, checked at iteration 1. clt's ratio drops below 3.0,
#    beating vterm's 3.7 and approaching relax's 2.7, because eps_rel and
#    steepness are free again. If it stays near 3.7 with them free, then
#    relax's clt improvement came from something else and we have been
#    chasing the wrong lever all along.
# 2. The weather floor falls. Run the spread-contraction test at the end
#    (parameter spread contracts about 12x over 6 iterations; a linear
#    response would contract G spread by the same factor). vterm's G spread
#    fell only 1.2-1.6x, giving floors of lwp 0.28, clt 0.136, swcre 0.235
#    sigma. Expect lower floors here. If they are unchanged, multi-year
#    sampling does not help and the fix is longer members.
# 3. v_snow is re-tested on independent data. It moved +1.45 prior sigma on
#    one October (0.96 -> 1.69 m/s) and -0.24 in the Southern Ocean. If it
#    again lands well above 1.4 m/s across five Septembers, the result is
#    real and not one month's weather. If it collapses toward 1.0, the
#    +1.45 was overfitting a single realization, which would be the most
#    important negative result of the campaign.
# 4. q_liq reproduces near 1.5e-4 (vterm gave 1.493e-4 on Oct 2010). A
#    posterior far from that across five Septembers would mean the
#    "low q_liq" answer, reproduced in four earlier runs, is seasonal or
#    weather-driven rather than a property of the model.
#
# Verdict (run 2026-08-05, job 7024987, 6 iterations, clean finish):
# Posterior: q_liq 1.329e-4, rain_tau 1807, eps_rel 0.0734, steepness
# 0.591, v_rain 5.436, v_snow 1.215.
# Displacements in prior sigma: eps_rel +2.07, v_snow +0.60, steepness
# -0.63, v_rain +0.54, q_liq -0.89, rain_tau +0.04.
# Loss 0.0433 -> 0.0282 at iteration 3 -> 0.0331, the same non-monotone
# shape as relax and vterm.
#
# 1. PASS. clt's ratio at iteration 1 was 2.2 against 3.7 with the levers
#    frozen and 2.7 in relax, and its spread rose 0.166 -> 0.259.
#    Per-parameter attribution was unambiguous: eps_rel gives clt 0.568
#    while every microphysics parameter gives 0.23-0.29.
# 2. FALSIFIED. The weather floor did NOT fall. Parameter spread contracted
#    14x while G spread fell 1.6x, leaving floors of lwp 0.287, clt 0.129,
#    swcre 0.188 sigma against vterm's 0.28, 0.136, 0.235 on ONE October.
#    MECHANISM: minibatch_size = 1 means each iteration still grades against
#    a single year. Minibatching over years removes the BIAS of overfitting
#    one realization; it does nothing to the VARIANCE within an iteration.
#    To lower the per-iteration floor, use minibatch_size > 1 (averaging G
#    over several years per update, which needs a multi-date forward model
#    and about 5x the workers) or longer members.
# 3. FALSIFIED, and it corrects a headline. v_snow moved +0.60 sigma to
#    1.215 m/s, not the +1.45 sigma to 1.69 seen on one October. Same sign,
#    less than half the size, and the October Southern Ocean gave -0.24.
#    The +1.45 was substantially one month's weather. QUOTE 1.22 m/s AND
#    +0.60 SIGMA, NOT THE OCTOBER NUMBER.
# 4. PASS. q_liq 1.329e-4, essentially icabl_sub's 1.34e-4 and 11 percent
#    below vterm's 1.493e-4. The low-q_liq answer now replicates across
#    season, across five years, and across three parameter sets.
#
# NEW RESULT, and read it as a WARNING rather than a discovery. eps_rel
# moved +2.07 prior sigma (0.0385 -> 0.0734), the largest displacement in
# the campaign. Do NOT read this as "the data wants more residual
# saturation variability". The g_vs_obs plots show what actually happened:
#   iteration 1  model clt tracks the observation across 0.3-0.95
#   iteration 6  model clt collapses to a flat 0.45-0.6 while the
#                observation still spans 0.05-1.0
# eps_rel is a FLOOR on subgrid saturation variability,
# sigma_S_floor^2 = (eps_rel * q_sat)^2 + sigma_abs^2, so raising it widens
# the assumed distribution and pushes cloud fraction TOWARD 0.5 from both
# sides: down where the box was near-saturated and cloudy, up where it was
# dry. A +2 sigma increase therefore HOMOGENISES the cloud field.
# Why the optimizer wants that: swcre started too reflective (its
# iteration-1 residuals are negative through the mid-latitudes, so mean(G)
# was more negative than observed). Flattening the cloud field removes
# cloud from the reflective storm tracks and improves swcre by 35 percent,
# and it pays for that with clt's spatial pattern.
# So an OPTICAL bias is being fixed with an AREAL knob. The right fix
# changes cloud optical properties (droplet number, effective radius; see
# the row 3 optics sweep), not cloud amount. eps_rel needs a tighter upper
# bound, or clt needs enough weight to veto the flattening. Do not quote
# eps_rel 0.0734 as a posterior without this caveat.
#
# CAVEAT ON THE FIX. Freeing eps_rel and steepness gave clt LEVERAGE but
# not a better FIT: clt's residual still degraded, 0.580 -> 0.749, while
# swcre improved 35 percent (0.858 -> 0.556) and lwp stayed flat. The
# optimizer spent the cloud-fraction parameters on swcre. So clt's
# remaining problem is PRIORITY, not reachability, and up-weighting clt is
# now the right lever where before it provably was not (a run with clt
# unreachable cannot be helped by weighting; this one can).
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_sep_multi.yml",
)

# Distinct-year September samples, ordered to span the range early so the
# first iterations see diverse years, then cycled to fill n_iterations + 1.
const _SEP_YEARS = [2006, 2010, 2008, 2007, 2009]
sample_date_ranges = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 9, 1)) for
    y in _SEP_YEARS[1 .+ ((0:6) .% length(_SEP_YEARS))]
]

# Covariance from the interannual spread of September 2006-2010.
const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 9, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_sep_multi")

# Same 5 degree grid as vterm, so the leverage numbers are comparable.
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

# Unchanged from vterm: same floors, same 800 km correlated floor.
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

include(joinpath(@__DIR__, "..", "prior_tools.jl"))

const CALIBRATION_PRIORS = [
    checked_constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        1.8e-4, 0.6e-4, 0.0, 1e-3,
    ),
    checked_constrained_gaussian("rain_autoconversion_timescale", 1800, 400, 300, 3600),
    checked_constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    checked_constrained_gaussian("cloud_fraction_steepness_scale", 0.7, 0.2, 0.1, 2.0),
    checked_constrained_gaussian("fixed_rain_terminal_velocity", 5, 1, 0, 10),
    checked_constrained_gaussian("fixed_snow_terminal_velocity", 1, 0.5, 0, 5),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
