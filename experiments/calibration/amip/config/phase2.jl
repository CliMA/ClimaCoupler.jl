# PHASE-2 joint calibration: 4 parameters x 4 observables.
#
# Design comes from the identifiability map (identifiability_map.md):
#   parameters  q_liq_threshold, rain_autoconversion_timescale (row 1),
#               detr_massflux_vertdiv_coeff, detr_buoy_coeff (rows 2, 2b)
#   observables lwp (MAC), pr (GPCP), swcre + lwcre (CERES)
# Each parameter has a designated identifying observable: lwp -> q_liq,
# pr -> rain_tau, {lwp, swcre} jointly -> dmfvd, swcre -> detr_buoy.
#
# Predictions stated before launch, to compare against the outcome:
#   1. q_liq stays near 3e-4 (the replicated warm-rain result).
#   2. dmfvd is pulled up toward 0.4-0.6 by swcre, trading a small lwp cost.
#   3. detr_buoy stays at or above 1 (values below 1 worsen swcre).
#   4. lwp and pr stay at their noise floors.
#   5. swcre improves by 0.3-0.6 sigma but plateaus above 1 sigma. The
#      remaining gap is structural and is covered by the larger noise floor.
#   6. lwcre stays near 1 sigma.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file =
    joinpath(pkgdir(ClimaCoupler), "config", "amip_configs", "amip_calibration.yml")

# Fixed calibration target: Oct 2010, repeated per minibatch slot. Length must be
# >= n_iterations + 1.
sample_date_ranges =
    [(Dates.DateTime(2010, 10, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:9]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 10, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_phase2")

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "pr", "swcre", "lwcre"],
    minibatch_size = 1,
    # 9 members run in parallel across 9 workers, so an iteration costs the same
    # wall time as before (~1.1 h). 8 iterations fit inside the 12 h walltime
    # with margin, which avoids the mid-iteration restart failure mode.
    n_iterations = 8,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Day(7),
    output_dir,
    rng_seed = 42,
)

# Per-group structural-error floors. lwp and pr keep the 20% floor that their
# calibrations converged under. swcre carries a 13-28 W/m^2 irreducible bias
# (30-50% of the field mean) that no swept parameter closes, so it gets a 50%
# floor. Without it the calibration would distort the reachable parameters to
# chase swcre. lwcre shares the radiation group.
const OBS_NOISE_GROUPS = [
    (short_names = ["lwp", "pr"], model_error_scale = 0.2),
    (short_names = ["swcre", "lwcre"], model_error_scale = 0.5),
]

const PRESSURE_LEVELS = 100.0 .* [200.0, 500.0, 850.0]
const ALTITUDE_LEVELS = [2000.0, 5000.0, 10000.0]

const NORMALIZATION_STATS_FP =
    joinpath(CALIBRATE_CONFIG.output_dir, "normalization_stats.jld2")

# Priors. Warm rain centers on the replicated posterior. dmfvd uses the row-2b
# result: the swcre gain per lwp cost is favorable below ~0.7, so the prior
# covers (0, 0.8) with its mass near the base value. detr_buoy gets room on
# both sides of 1 because its response saturates above 1.
const CALIBRATION_PRIORS = [
    PD.constrained_gaussian(
        "cloud_liquid_water_specific_humidity_autoconversion_threshold",
        3e-4, 1.5e-4, 0.0, 1.5e-3,
    ),
    PD.constrained_gaussian("rain_autoconversion_timescale", 1450, 500, 300, 3600),
    PD.constrained_gaussian("detr_massflux_vertdiv_coeff", 0.35, 0.15, 0.0, 0.8),
    PD.constrained_gaussian("detr_buoy_coeff", 1.0, 0.6, 0.05, 5.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
