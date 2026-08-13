# PHASE-3 merged calibration: 6 parameters x 5 observables, zonal.
#
# Merges the two independently validated sets from the identifiability map:
#   phase2       {q_liq_threshold, rain_autoconversion_timescale,
#                 detr_massflux_vertdiv_coeff, detr_buoy_coeff}
#   cls zonal    {cloud_fraction_eps_rel, cloud_fraction_steepness_scale}
# snow_autoconversion_timescale is dropped: it barely contracted in the cls
# runs (0.86x prior). Observables add cl to the phase2 four; the lwcre floor
# comes from the phase2_lwcre025 verdict (0.25).
#
# The one unmeasured risk is cross-set degeneracy: dmfvd and the
# cloud-fraction pair both touch swcre and cl but were never swept together.
# The steering degen line monitors this live.
#
# Predictions (falsifiable):
# 1. No collapse; all five observables at their floors by iteration ~5
#    (lwp/pr/swcre near 1 sigma, lwcre 0.6 to 1.1, cl 0.8 to 1).
# 2. Posteriors agree with the donor runs within one final spread:
#    q_liq 2 to 3e-4, dmfvd 0.4 to 0.6, detr_buoy at or above 1,
#    eps_rel 0.04 to 0.06, steepness 0.55 to 0.75.
# 3. All parameter pairs stay separable, in particular dmfvd vs eps_rel and
#    steepness. If a pair goes degenerate, stop and design a targeted
#    two-parameter sweep instead of rerunning.
# 4. The added lwcre + cl information tightens the detrainment posteriors
#    below phase2's final spreads (0.58 and 0.73 x prior).
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:9]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_phase3")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "pr", "swcre", "lwcre", "cl"],
    minibatch_size = 1,
    n_iterations = 8,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Floors from the validated runs: lwp/pr 0.2 (phase2), swcre 0.5 (irreducible
# brightness bias), lwcre 0.25 (phase2_lwcre025), cl 0.5 (cls zonal).
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

# Priors reuse the donor runs' (not their posteriors): the merged run should
# re-derive the answer from the data, not inherit confidence.
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
