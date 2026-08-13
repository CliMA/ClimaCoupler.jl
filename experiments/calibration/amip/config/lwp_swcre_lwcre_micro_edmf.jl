# First science calibration on the post-rebase physics: 1M microphysics +
# ice nucleation + EDMF against lwp / swcre / lwcre.
#
# Context. The release_repro check (2026-08-12) proved the migrated pipeline
# is bit-exact but the forward model changed: cloud_ice_formation:
# TemperatureDependent (commit 674ccf36) plus ClimaAtmos 0.42.3 -> 0.42.4
# moved G by 23-41% rel RMS at identical parameters. Campaign posteriors are
# therefore priors-to-reverify, and the new ice-formation pathway has knobs
# (Frostenberg2023) the campaign never saw. Sigma units and residuals are
# NOT comparable to campaign-era runs.
#
# MODEL (decided 2026-08-13): ClimaAtmos#main, pinned by the manifest
# (repo-rev main, tree 5ee6833a, declares 0.42.4), with ClimaCore 0.15.1
# and ClimaParams 1.1.4 - chosen over the repro-validated registered
# 0.42.4 to calibrate the current development head. The environment IS
# experiments/AMIP/Manifest-v1.11.toml: commit it alongside this config,
# and do not update packages while the run or its analysis is in flight.
# This stack has never run a member on this branch; preflight (wiring
# check included) must be green on it before launch.
#
# Loss: lwp (MAC), swcre and lwcre (CERES). clt is dropped on the campaign's
# closing result (its pattern error is directionally orthogonal to the
# parameter space at this rung; edmf + panom verdicts). lwcre is new to the
# 2-D era: it constrains the high-cloud/ice pathway that the freed ice and
# snow parameters touch, where lwp/swcre are mostly liquid-side.
#
# Parameters (8 -> 17 TransformUnscented members): warm/cold 1M microphysics
# (rain & snow autoconversion timescales, sublimation/deposition timescale,
# snow terminal velocity), ice nucleation (Frostenberg2023 a), and EDMF
# (entrainment, massflux-vertdiv detrainment, eddy viscosity) - a direct
# retest of the edmf-run null (all EDMF params flat) under a different loss
# and different physics.
#
# TOML layering (member last, later file wins): toml/amip_progedmf_1m.toml
# (the PRODUCTION parameter set, same file amip.yml uses; content-identical
# to the campaign's wxquest_progedmf.toml) then the member's
# parameters.toml. No calibration fixed TOML at all: the old-physics
# posteriors of calibration_release_fixed.toml are deliberately dropped, so
# the baseline is exactly the shipped model and posteriors upstream 1:1
# into the same file. Net non-calibrated operating point: q_liq threshold
# 2e-4, cond_evap 100 s, steepness 0.55 (all in the production TOML);
# v_ice 0.01 and v_rain 5.0 (ClimaParams defaults - a DELIBERATE departure
# from the campaign-era pins 0.1 / 5.4; v_ice is 10x lower than any
# campaign run ever used, in the ice pathway lwcre grades - watch the
# iteration-1 lwcre residual). Prior means re-center at the production
# operating point / ClimaParams 1.1.4 defaults, not at campaign posteriors
# (those belong to the old physics). NOTE entr_coeff prior mean 1 vs the
# production pin 0.1 and dmfvd 1 vs 0.3: ClimaParams 1.1.4 defaults are
# 1 - confirm the current entrainment closure expects O(1) coefficients
# before launch (do not assume campaign-era scaling).
#
# Predictions (DRAFT - edit before launch, then record a verdict after):
# 1. Go/no-go at iteration 1: every observable's leverage ratio
#    (residual / ensemble response spread) is <= ~3. lwcre in particular is
#    informative at its 0.25 floor (phase2 measured 0.5 leaves it inert).
# 2. The new ice pathway is identifiable: Frostenberg2023_a or subl_dep
#    displaces >= 0.5 prior sigma by the final iteration. Both flat means
#    TemperatureDependent ice formation is invisible to lwp/swcre/lwcre at
#    5 degrees and SCM/profile observables are needed for it.
# 3. The edmf null is retested: at least one of entr / dmfvd / eddy_visc
#    displaces >= 0.5 prior sigma. All three flat reproduces the edmf-run
#    conclusion under the new loss and physics.
# 4. No trade-away: no observable ends more than the weather floor
#    (~0.13 sigma) above its iteration-1 residual.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_micro_edmf.yml",
)

# Five distinct Septembers, minibatched one per iteration (sep_multi
# protocol): kills the overfit-one-year bias without changing the
# per-iteration weather floor. 7-day spinup from the son_ic Aug 25 ICs,
# September graded (extend = 1 month past the Sep 1 sample date).
const _SEP_YEARS = [2006, 2010, 2008, 2007, 2009]
sample_date_ranges = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 9, 1)) for
    y in _SEP_YEARS[1 .+ ((0:6) .% length(_SEP_YEARS))]
]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 9, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_micro_edmf")

const COARSEN_FACTOR = 2

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "swcre", "lwcre"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Floors: lwp 0.6 and swcre 0.75 are the release-era 2-D values. swcre is
# full-field again (the zsw zonal reduction existed to stop swcre pattern
# error financing clt removal; clt is not in this loss). lwcre 0.25 is the
# phase2_lwcre025 value, kept through phase5: 0.5 measured lwcre inert,
# 0.25 made it a real constraint - same transfer clt's scale made between
# representations. All three keep the 800 km correlated floor; diagonal
# floors on 2-D grids over-inform and collapse the ensemble (lwp_pr_2d10/15).
const OBS_DECORRELATION_LENGTH = 8.0e5
const OBS_NOISE_GROUPS = [
    (
        short_names = ["lwp"],
        model_error_scale = 0.6,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["swcre"],
        model_error_scale = 0.75,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
    (
        short_names = ["lwcre"],
        model_error_scale = 0.25,
        decorrelation_length = OBS_DECORRELATION_LENGTH,
    ),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

include(joinpath(@__DIR__, "..", "prior_tools.jl"))

# checked_constrained_gaussian, not the plain EKP constructor: it verifies
# the achieved constrained mean/std against the requested ones and errors on
# silent drift (the q_liq wide-prior incident). Values as requested.
const CALIBRATION_PRIORS = [
    checked_constrained_gaussian("rain_autoconversion_timescale", 1200, 300, 100, 7200),
    checked_constrained_gaussian("snow_autoconversion_timescale", 1200, 300, 100, 7200),
    checked_constrained_gaussian("Frostenberg2023_a_coefficient", 1, 0.3, 0, 10),
    checked_constrained_gaussian("sublimation_deposition_timescale", 300, 100, 0, 7200),
    checked_constrained_gaussian("fixed_snow_terminal_velocity", 1, 0.3, 0.1, 3),
    checked_constrained_gaussian(
        "mixing_length_eddy_viscosity_coefficient",
        0.2, 0.1, 0, 1.0,
    ),
    checked_constrained_gaussian("entr_coeff", 1, 0.3, 0, 10.0),
    checked_constrained_gaussian("detr_massflux_vertdiv_coeff", 1, 0.3, 0, 10.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
