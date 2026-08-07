# January Southern Ocean, lwp + swcre, ZONAL reduction, lwp up-weighted.
# Follow-up to lwp_clt_swcre_so_band.jl, which failed on leverage in
# October. Four things change, each for a measured reason.
#
# 1. ZONAL REDUCTION instead of a 2-D grid. No COARSEN_FACTOR, so
#    reduce_spatial falls back to zonal_average: 8 latitude values per
#    observable in the band. Measured on vterm's existing G ensembles, in
#    physical units, collapsing longitude improves signal-to-noise by 1.7x
#    for lwp, 2.7x for clt and 2.9x for swcre, and lifts clt and swcre above
#    1 for the first time in this campaign:
#
#      obs          weather floor   parameter signal   S/N
#      lwp 2-D            0.0168             0.0231   1.37
#      lwp zonal          0.0049             0.0116   2.38
#      swcre 2-D            6.11               5.85   0.96
#      swcre zonal          1.25               3.53   2.83
#
#    The floor comes from the iteration-6 spread (parameter spread has
#    contracted about 12x by then, so what remains is weather) and the
#    signal by subtracting it in quadrature from iteration 1. Averaging
#    keeps the coherent parameter response and removes incoherent weather.
#
#    WARNING when reading the results: plain residual/spread gets WORSE
#    under averaging (lwp 1.6 -> 2.4) purely because averaging removes
#    weather from the denominator. Divide by the parameter signal, not the
#    total spread, whenever the reduction changes. Noise-corrected, zonal
#    is flat against 2-D: lwp 1.9 -> 2.6, clt 7.4 -> 7.5, swcre 4.1 -> 4.4.
#
# 2. JANUARY, not October. MAC lwp is passive microwave over ocean and sea
#    ice breaks the retrieval. October sits near the Antarctic sea-ice
#    maximum, January near the minimum. Valid lwp points in 45S-65S at
#    2.5 degrees: October 762, January 931, February 1119, with 65-60S
#    going from 12 percent in October to 41 in January. February would be
#    better still but needs a 31-day spinup from the Jan 1 IC and 59-day
#    members; Dec 25 2009 exists in wxquest_initial_conditions, so January
#    keeps the 7-day spinup and 38-day members of every earlier run.
#    Southern summer also means high insolation, so swcre should carry real
#    information here, unlike September (noise 92.8 percent of signal).
#
# 3. clt REMOVED. Three runs showed these parameters cannot move it
#    (ratio 3.7 globally, 5.7 in the band, and a clt-only run went nowhere
#    in two iterations). Its levers, eps_rel and steepness, are fixed here
#    because with clt gone nothing would constrain them.
#
# 4. lwp UP-WEIGHTED, model_error_scale 0.6 -> 0.3, so 4x the weight per
#    point. Weighting cannot make an observable reachable, since
#    residual/spread has sigma in both terms. It DOES decide which
#    observable the update serves when two conflict, which is exactly how
#    clt lost its fit in vterm. lwp earns the priority: it had the best
#    ratio everywhere, 1.5 globally and 3.7 in the band against swcre's
#    2.6 and 3.7.
#
# Only 16 constraints (8 zonal values x 2 observables) for 4 parameters.
# Formally over-determined 4:1, but thin. If the update looks
# prior-dominated, the fix is longitude sectors (say 4 sectors x 8
# latitudes = 32 values), which keeps most of the averaging gain and some
# east-west contrast, and needs a small addition to preprocessing.jl.
#
# Predictions (falsifiable):
# 1. THE GO/NO-GO at iteration 1. lwp's NOISE-CORRECTED ratio is below 3.
#    The October band gave 3.7 uncorrected on 2-D data. If January plus
#    zonal still leaves lwp above 4, the Southern Ocean is unresponsive to
#    warm-rain parameters year round and the region should be abandoned for
#    this parameter set.
# 2. S/N above 1 for both observables, measured with the spread-contraction
#    test at the end. The zonal gain should carry over: rings in the band
#    are full 144-longitude rings, so the averaging is the same as in the
#    global test. If S/N stays below 1, zonal reduction does not deliver
#    here and the noise-limited problem needs longer members instead.
# 3. q_liq lands near 1.5e-4 again. That value has now been reproduced on
#    October 2010 with 6 free parameters (1.493e-4) and in four earlier 2-D
#    runs. A very different answer here would mean it is seasonal or
#    regional rather than a property of the model.
# 4. v_snow barely moves. It went +1.45 prior sigma globally in October and
#    -0.24 in the October Southern Ocean. Midsummer Southern Ocean has
#    little snow, so expect little movement. Strong movement would mean the
#    snow signal is not seasonal after all.
#
# Verdict (run 2026-08-06, job 7026239, 6 iterations, clean finish):
# NULL RESULT. The run learned nothing. Parameter spread contracted to
# 0.9926 of iteration 1, that is, not at all. Displacements: q_liq +0.03,
# rain_tau +0.01, v_rain +0.01, v_snow -0.01 prior sigma. Loss 0.0244 ->
# 0.0327 with the best value at iteration 1. The posterior is the prior.
#
# CAUSE. The observation carries no information under its own noise model.
# Loss 0.0244 times 4 parameters gives a whitened misfit of 0.098, an RMS
# residual of 0.31 of the assumed noise, the same "already fit" condition
# that made clt inert. Two choices combined to produce it:
#   1. Zonal averaging cancels opposing biases within each ring, shrinking
#      the residual, while model_error_scale keeps scaling the floor with
#      field magnitude. lwp's residual is 0.66 sigma on the diagonal but
#      0.31 whitened.
#   2. A 400 km correlated floor over 8 latitudes spanning 2200 km leaves
#      few effective constraints. Adjacent rings are 278 km apart, so
#      exp(-278/400) = 0.50, giving about 2200/400 = 5 independent
#      constraints per observable, roughly 11 against 4 parameters.
#
# THE FIX for a rerun: with zonal data the correlated floor DOUBLE COUNTS,
# because averaging already removed the small-scale noise that the
# correlation exists to model. Use a DIAGONAL floor (omit
# decorrelation_length) and drop model_error_scale, lwp to about 0.2.
#
# Predictions 1-4 are all unevaluable rather than passed or failed: with no
# ensemble movement there is nothing to judge. v_snow's -0.01 sigma
# "confirms" prediction 4 only vacuously.
#
# LIMIT OF THE DIAGNOSTIC, worth knowing before trusting it elsewhere. The
# spread-contraction estimate of the weather floor assumes the parameter
# spread shrinks enough that the final spread is pure weather. Here it did
# not shrink, so sqrt(spread1^2 - spreadN^2) returns about zero signal by
# construction and the reported S/N of 0.70 and 0.59, and the corrected
# ratios of 5.5 and 11.8, are ARTEFACTS. Require contraction below about
# 0.3 before believing that estimate.
#
# The zonal-versus-2D question is therefore NOT settled by this run. The
# measurement on vterm's own ensembles (see the header note above) remains
# the evidence, and it still favours reduction.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_so_jan.yml",
)

sample_date_ranges =
    [(Dates.DateTime(2010, 1, 1), Dates.DateTime(2010, 1, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 1, 1), Dates.DateTime(y, 1, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_so_jan")

# Southern Ocean only. Both the observation and the simulation read this.
const LAT_WINDOW = (-65, -45)

# The covariance regularization takes a quantile of the per-point model
# error scale and needs at least 1/qtl values per variable. The default
# 0.05 needs 20; a zonal mean over this 20 degree band leaves 8, which
# fails with "Insufficient samples for computing quantile". 0.2 needs 5.
# With 8 values a 20 percent quantile still sits near the low end of the
# floor magnitudes, which is what the regularization wants.
const OBS_REGULARIZATION_QUANTILE = 0.2

# NO COARSEN_FACTOR on purpose: reduce_spatial then zonal-averages, which is
# the whole point of this run. See the note above.

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "swcre"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# lwp at half its usual floor: 4x the weight per point. 400 km because the
# band is only about 2200 km wide, so an 800 km floor would treat a
# band-wide bias as expected error.
const OBS_DECORRELATION_LENGTH = 4.0e5
const OBS_NOISE_GROUPS = [
    (
        short_names = ["lwp"],
        model_error_scale = 0.3,
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
    checked_constrained_gaussian("fixed_rain_terminal_velocity", 5, 1, 0, 10),
    checked_constrained_gaussian("fixed_snow_terminal_velocity", 1, 0.5, 0, 5),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
