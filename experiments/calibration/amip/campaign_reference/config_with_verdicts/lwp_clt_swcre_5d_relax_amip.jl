# AMIP-mode rerun of lwp_clt_swcre_5d_relax: same 8 parameters and floors,
# but free-running AMIP members instead of subseasonal ERA5-initialized
# ones. Members run 4 months from Jul 1 2010 (atmos default IC, AMIP
# SST/SIC); the first month is spinup, so the samples are the Aug/Sep/Oct
# 2010 monthly means. This run also actually uses the tightened q_liq
# prior; the subseasonal relax run accidentally ran with the old wide one
# (see the incident note in lwp_clt_swcre_5d_relax.jl). VERIFY priors.png
# at driver init before letting workers start members.
#
# Predictions (falsifiable):
# 1. With the tightened prior the error trace shows learning, not the
#    iteration-1 jump: q_liq drift stays below 2x spread/iter at iter 1
#    and the posterior lands in 1.6-2.4e-4.
# 2. IC-mode robustness: the posterior lands in the subseasonal relax
#    run's neighborhood (margin ~1.0-1.3, cond_evap ~100 s, subl_dep
#    ~280 s, eps_rel 0.04-0.06). Large shifts mean the subseasonal
#    posteriors were tuned to ERA5 initialization, not climate.
# 3. Per-month residuals stay flat Aug -> Oct (a free-running AMIP has no
#    reanalysis state to drift from, unlike phase5's Sep -> Nov growth).
#    Growth here instead suggests the parameters trade transient vs
#    equilibrated cloud states.
# 4. Margin identifiability replicates: spread below 0.5x prior by the
#    final iteration.
#
# Select via ENV["CALIBRATION_CONFIG"].

config_file = joinpath(
    pkgdir(ClimaCoupler),
    "config",
    "amip_configs",
    "amip_calibration_amip_mode.yml",
)

sample_date_ranges =
    [(Dates.DateTime(2010, 8, 1), Dates.DateTime(2010, 10, 1)) for _ in 1:7]

const COVARIANCE_DATE_RANGES = [
    (Dates.DateTime(y, 8, 1), Dates.DateTime(y, 10, 1)) for y in 2006:2010
]

output_dir = joinpath(
    pkgdir(ClimaCoupler),
    "amip_calibration_lwp_clt_swcre_5d_relax_amip",
)

const COARSEN_FACTOR = 2

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file,
    short_names = ["lwp", "clt", "swcre"],
    minibatch_size = 1,
    n_iterations = 6,
    sample_date_ranges,
    extend = Dates.Month(1),
    spinup = Dates.Month(1),
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
    checked_constrained_gaussian("snow_autoconversion_timescale", 1800, 400, 300, 3600),
    checked_constrained_gaussian("condensation_evaporation_timescale", 150, 100, 20, 3600),
    checked_constrained_gaussian("sublimation_deposition_timescale", 500, 400, 20, 3600),
    checked_constrained_gaussian("cloud_fraction_floor_release_abs_margin", 1, 0.5, 0, 2),
    checked_constrained_gaussian("cloud_fraction_eps_rel", 0.04, 0.02, 0.001, 0.2),
    checked_constrained_gaussian("cloud_fraction_steepness_scale", 0.7, 0.2, 0.1, 2.0),
]

const PRIORS = EKP.combine_distributions(CALIBRATION_PRIORS)
