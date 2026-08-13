# EDMF calibration: the first parameter set that can move the clt
# PATTERN, against the zsw loss.
#
# Why. The release run proved clt is out of reach of every cloud-closure
# and optics scalar at its posterior (ratio 5.6, spread 0.168 sigma): a
# global scalar fixes amplitude, and clt's remaining error is spatial
# pattern. Mining the row-2 EDMF sweeps (validated against the release
# ensemble: max-overlap cover proxies clt responses at correlation 0.81,
# and the anchor row reproduces the leverage attribution) showed all
# EDMF candidates move column cover 7.8-10.3 percentage-point RMS at
# plausible ranges, comparable to eps_rel's 12.4. clt is reachable
# through EDMF.
#
# Free parameters, ranked by mined cover response:
#   detr_massflux_vertdiv_coeff 10.3, static_stab 8.6, entr_coeff 8.4,
#   eddy_viscosity 8.3, detr_buoy 7.8, EDMF_max_area 7.8.
# Ri_crit duplicates static_stab (8.5) and is left fixed.
# mixing_length_diss_coeff is EXCLUDED: swept 0.073-0.66 in row 2 with
# exactly zero field response; suspect it is not wired. Do not free it
# without a wiring check.
#
# Priors: means at the campaign defaults (wxquest_progedmf.toml), sigma
# about half the sweep half-range, bounds generous but physical.
#
# The loss is the zsw loss UNCHANGED (lwp and clt full-field at 5
# degrees, swcre zonal with diagonal floor): one change at a time. The
# mining cannot separate pattern from amplitude in the EDMF swcre
# responses, so restoring full-field swcre could re-fund the
# cloud-removal trade; it returns only if this run shows the EDMF
# parameters fix clt pattern without spending swcre.
#
# Fixed parameters: toml/calibration_edmf_fixed.toml freezes the six
# closure/optics parameters at the zsw posterior and the settled
# microphysics (see that file's header).
#
# Predictions (falsifiable):
# 1. THE POINT. clt improves for the first time in the campaign: its
#    final residual ends BELOW its iteration-1 value, and by more than
#    the weather-floor noise (about 0.13 sigma). Every previous run
#    ended with clt flat or degraded.
# 2. The improvement is pattern, not amplitude: in the g_vs_obs strips
#    the iteration-6 model clt tracks the observation across its range
#    instead of flattening. (Judged from the plots; the run cannot fake
#    this with a mean shift because the noise floor already absorbs it.)
# 3. swcre's zonal amplitude does not degrade: final zonal swcre
#    residual within the weather floor of its iteration-1 value.
# 4. At least two EDMF parameters are identifiable (displacement at
#    least 0.5 prior sigma). All six flat means September AMIP at 5
#    degrees cannot see the EDMF state and the SCM rung is the only
#    place to calibrate it.
#
# VERDICT (2026-08-08, 6 iterations, driver exited 0, ClimaAtmos 0.42.3;
# iterations 1 and 6 both grade Sep 2006, so endpoints are same-weather).
# Residuals FLAT: lwp 0.65 -> 0.68, clt 0.84 -> 0.86, swcre 0.76 -> 0.78,
# all moves inside the weather floor. Posterior displacements in prior
# sigma: detr_buoy +0.48, static_stab +0.39, entr -0.27, dmfvd +0.12,
# max_area +0.05, eddy_visc +0.04.
# P1 FAIL: clt did not improve. P2 moot. P3 PASS but trivially: nothing
#    moved, so nothing was traded. P4 FAIL: no displacement reaches 0.5.
# The iteration-1 leverage explains it: clt response spread 0.23 sigma
# against a 0.84 residual, ratio 3.7 = OUT OF REACH, worse than the
# closure/optics scalars had in release (2.6). The mining ranking was
# right in direction (dmfvd is the family's strongest clt lever, 0.41)
# and wrong in sufficiency: single-parameter responses at sweep ranges
# do not add up to reach at prior ranges under this noise model.
# CONCLUSION: September AMIP at 5 degrees cannot see the EDMF state.
# With warm rain, closure, optics, and now EDMF all exhausted, clt sits
# at its assumed floor with no parameter family able to move it: the
# AMIP-rung scalar campaign is CONVERGED at this noise model. EDMF
# calibration moves to the SCM rung (profile observations, per-regime
# forcing); the AMIP rung re-enters only with a sharper observable
# (pattern metrics, more months, or a lower clt floor justified by a
# better simulator).
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_edmf.yml",
)

const _SEP_YEARS = [2006, 2010, 2008, 2007, 2009]
sample_date_ranges = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 9, 1)) for
    y in _SEP_YEARS[1 .+ ((0:6) .% length(_SEP_YEARS))]
]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 9, 1), Dates.DateTime(y, 9, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_edmf")

const COARSEN_FACTOR = 2
const ZONAL_SHORT_NAMES = ["swcre"]

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
    (short_names = ["swcre"], model_error_scale = 0.3),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

include(joinpath(@__DIR__, "..", "prior_tools.jl"))

const CALIBRATION_PRIORS = [
    checked_constrained_gaussian("entr_coeff", 0.1, 0.06, 0.005, 0.5),
    checked_constrained_gaussian("detr_buoy_coeff", 1, 0.6, 0.05, 5),
    checked_constrained_gaussian("detr_massflux_vertdiv_coeff", 0.3, 0.2, 0.02, 1.5),
    checked_constrained_gaussian("EDMF_max_area", 0.7, 0.12, 0.3, 0.99),
    checked_constrained_gaussian(
        "mixing_length_eddy_viscosity_coefficient",
        0.14, 0.08, 0.01, 0.6,
    ),
    checked_constrained_gaussian("mixing_length_static_stab_coeff", 0.4, 0.25, 0.02, 2.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
