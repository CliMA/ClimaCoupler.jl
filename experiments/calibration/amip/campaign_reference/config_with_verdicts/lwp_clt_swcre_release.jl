# Floor-release calibration: give the optimizer a way to fix the shortwave
# bias WITHOUT removing cloud. Designed from reading the cloud-fraction
# closure (ClimaAtmos 0.42.2 cloud_fraction.jl, byte-identical on main).
#
# What sep_multi showed. The optimizer fixed swcre (too reflective, -35
# percent residual) by raising eps_rel +2.07 prior sigma, which lowers CF
# monotonically wherever the floor is active: sigma_aug rises, C = q_c /
# sigma_aug falls, CF = Phi(z(C)) falls. Model clt collapsed from tracking
# 0.3-0.95 to a flat 0.45-0.6. An OPTICAL bias was fixed with an AREAL knob.
#
# Why eps_rel was the only lever it had: eps_rel appears in BOTH the floor
# magnitude (sigma_floor ~ eps_rel*q_sat) and the release width
# (w ~ abs_margin*eps_rel*q_sat), so raising it also raises the saturation
# margin needed to release the floor, capping CF even in equilibrated
# overcast decks - exactly the failure the closure docstring warns about.
# With the release shape frozen the coupling cannot be broken.
#
# The fix, two-pronged:
#   1. Free the four release-shape parameters (margin, abs_margin,
#      sharpness, residual) so patchiness in subsaturated air no longer
#      forces a CF cap in saturated decks.
#   2. Free prescribed_cloud_droplet_number_concentration, the true optical
#      knob: it is N0 in ml_N_cloud_liquid_droplets, multiplying the
#      aerosol-modulated droplet number that feeds the Liu-Hallett effective
#      radius in RRTMGP. Lower Nd = larger droplets = dimmer clouds at
#      identical LWP and CF. It does not touch 1M microphysics.
# And clt is up-weighted (floor 0.35 -> 0.25) so it can veto cloud removal.
# sep_multi established clt HAS leverage with eps_rel free (ratio 2.2) and
# loses on priority; weighting is the right lever now, where in vterm it
# provably was not.
#
# steepness_scale is deliberately FIXED at 0.55 despite past calibration:
# alpha = 1/steepness multiplies the whole augmented sigma (degenerate with
# eps_rel; sep_multi moved both the same direction, doubling the effective
# floor alpha*eps_rel) and leaks into the microphysics condensate closure
# via sigma_S_eff. See toml/calibration_release_fixed.toml.
#
# Everything else copies sep_multi, the best-behaved protocol so far: five
# distinct Septembers 2006-2010 (minibatch over years kills the overfit-
# one-year bias; it does NOT lower the per-iteration weather floor, that
# was measured), 5 degree grid, 800 km correlated floor, 7-day spinup from
# the son_ic Aug 25 ICs. 6 free parameters -> 13 members.
#
# Predictions (falsifiable):
# 1. THE POINT. swcre improves without clt degrading: swcre residual ends
#    at or below sep_multi's 0.556 AND clt's ends at or below its 0.580
#    iteration-1 value, instead of the 0.749 it degraded to. If clt still
#    degrades with the release free, Nd available, and clt up-weighted,
#    the cloud-removal pressure is not a parameterization artifact and
#    points at the EDMF state itself.
# 2. The Nd pathway is used: Nd moves at least 0.3 prior sigma (its
#    attribution row at iteration 1 will already show whether swcre
#    responds to it more than to eps_rel).
# 3. eps_rel stays within 1 sigma of its prior (0.02-0.06) now that the
#    release shape can decouple deck-capping from patchiness. Another +2
#    sigma excursion means the degeneracy was not broken.
# 4. The release parameters are identifiable as a group: at least one of
#    margin/abs_margin/sharpness/residual moves 0.5 prior sigma. All four
#    flat means clt cannot see the release shape at 5 degrees and the
#    ClimaParams promotion of the margin should be reconsidered.
#
# Iteration-1 leverage (go/no-go, PASSED): ratios lwp 2.1, clt 2.6,
# swcre 2.4 (sep_multi: 1.5, 2.2, 2.5). clt residual 0.758 sigma is NOT
# worse than sep_multi's 0.580: the clt noise scale here is 0.25 against
# 0.35, and 0.758 * 0.25/0.35 = 0.54 in sep_multi units, the same start.
# Apply that conversion before judging prediction 1. Attribution: eps_rel
# is still the top clt lever (0.751) but Nd shows the wanted asymmetry,
# strong on lwp 0.488 and swcre 0.410, weak on clt 0.274, so the optical
# pathway can fix swcre without spending cloud area. Each release-shape
# parameter responds 0.27-0.42 on all three observables.
#
# Verdict (run 2026-08-06/07, tmux driver, 6 iterations, clean finish):
# Posterior: eps_rel 0.0583, release_margin 0.848, release_abs_margin
# 1.29, release_sharpness 0.826, floor_residual 0.0333, Nd 4.11e7.
# Displacements in prior sigma: eps_rel +1.17, Nd -0.94, abs_margin
# +0.56, margin -0.19, residual -0.19, sharpness -0.16.
# Error: lwp 0.65 -> 0.66, clt 0.76 -> 0.93, swcre 0.69 -> 0.55 (this
# run's sigma units). Weather floor validly measured (contraction 14x
# params, 1.76x G).
#
# 1. HALF PASS. swcre 0.553 <= 0.556. clt degraded from its own start,
#    0.76 -> 0.93 (+22 percent), against sep_multi's +29. In ABSOLUTE
#    units the release posterior ends with clt error 11 percent BELOW
#    sep_multi's end at equal swcre, so the trade got cheaper but did
#    not vanish. The stated consequence applies: with the release shape
#    free, Nd available, and clt up-weighted, the leftover cloud-removal
#    pressure is not a closure artifact and points at the EDMF state.
#    Confirmed by the posterior leverage: clt ratio 5.6, spread 0.168
#    sigma, OUT OF REACH of all six parameters at the posterior.
# 2. PASS, decisively. Nd fell -0.94 sigma (9.1e7 -> 4.1e7 per m^3),
#    the second-largest move. The optical pathway was used: dimmer
#    clouds through larger droplets instead of fewer clouds.
# 3. MARGINAL. eps_rel +1.17 sigma (0.0385 -> 0.0583), above the 1
#    sigma line but nearly half of sep_multi's +2.07. Freeing the
#    release shape relieved most of the excursion; the coupling is
#    weakened, not fully broken.
# 4. PASS. abs_margin, the parameter added to break the eps_rel
#    coupling, moved +0.56 sigma to 1.29. The others moved 0.16-0.19
#    sigma. The release shape is identifiable at 5 degrees, led by
#    exactly the term the closure reading predicted.
#
# Follow-up pointers: quote Nd 41 per cm^3 with the caveat that it is
# marine-pristine low and optics may be absorbing an EDMF cloud-state
# error; next clt lever must come from EDMF parameters (entrainment /
# detrainment), not microphysics or the CF closure.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_release.yml",
)

const _SEP_YEARS = [2006, 2010, 2008, 2007, 2009]
sample_date_ranges = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 9, 1)) for
    y in _SEP_YEARS[1 .+ ((0:6) .% length(_SEP_YEARS))]
]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 9, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_release")

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

# clt 0.35 -> 0.25: with leverage established (sep_multi ratio 2.2), clt's
# remaining problem is priority against swcre, and this is the veto.
const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [
    (
        short_names = ["lwp"],
        model_error_scale = 0.6,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["clt"],
        model_error_scale = 0.25,
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
    checked_constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    # Release shape (see _compute_cloud_fraction): margin scales the
    # equilibrium-width term of the release width, abs_margin the
    # eps_rel*q_sat term (freeing it breaks the eps_rel coupling),
    # sharpness the transition exponent. Defaults are all 1.
    checked_constrained_gaussian("cloud_fraction_floor_release_margin", 1, 0.5, 0, 3),
    checked_constrained_gaussian("cloud_fraction_floor_release_abs_margin", 1, 0.5, 0, 2),
    checked_constrained_gaussian("cloud_fraction_floor_release_sharpness", 1, 0.5, 0.2, 4),
    # residual caps CF below 1 inside fully-released saturated decks, so
    # the prior keeps its mass at small values and the bound stays well
    # under where it visibly caps stratocumulus. sigma 0.04 not 0.05:
    # EKP rejects mu - sigma at the lower bound.
    checked_constrained_gaussian("cloud_fraction_floor_residual", 0.05, 0.04, 0, 0.3),
    # The optical knob: N0 multiplying the aerosol-modulated droplet number
    # behind the Liu-Hallett effective radius. Prior from row3_optics.
    checked_constrained_gaussian(
        "prescribed_cloud_droplet_number_concentration",
        1e8, 8e7, 1e6, 1e9,
    ),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
