# SON 2010 cloud radiative effects as zonal means, plus the global-mean net downward
# radiation at the top of the atmosphere (rtmt). Without rtmt this loss moved net TOA from
# +1.0 to +3.5 W m^-2 (CERES +1.4). Members start 2010-07-01 (July and August are spinup) and
# end 2010-12-01. The yml uses "amip" mode, whose SST lasts beyond the two months
# "subseasonal" mode covers, and a 120 s step (180 s went unstable on day 57).

# Seasonal means of the zonal mean and of the global mean.
seasonal_zonal = seasonal_average ∘ zonal_average
seasonal_global = seasonal_average ∘ global_mean

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config",
        "amip_configs",
        "amip_calibration.yml",
    ),
    # One line per loss field. The model errors (32% and 19% of the observed mean) were set
    # for 10-degree cells, where the flagship scores about 1.3 per point; re-derive them for
    # zonal means before relying on this config. rtmt is one value, so its noise sets its
    # weight: 0.2 W m^-2 is about one run's spread, and a 2 W m^-2 miss costs 100.
    observables = [
        (name = "swcre", reduce = seasonal_zonal, model_error_fraction = 0.316),
        (name = "lwcre", reduce = seasonal_zonal, model_error_fraction = 0.187),
        (name = "rtmt", reduce = seasonal_global, noise_sd = 0.2),
    ],
    # The SON 2010 calibration converged after two updates, and two 5.5 h iterations fit the
    # default 12 h worker walltime. A season is dated by its first month.
    n_iterations = 2,
    sample_date_ranges = [(DateTime(2010, 9, 1), DateTime(2010, 9, 1))],
    covariance_date_ranges = [(DateTime(y, 9, 1), DateTime(y, 9, 1)) for y in 2001:2010],
    extend = Month(3),
    spinup = Month(1),
    # On Derecho, make this a symlink to scratch: member output fills the home quota.
    output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_son2010_cre_toa"),
    priors = EKP.combine_distributions([
        PD.constrained_gaussian(
            "cloud_liquid_water_specific_humidity_autoconversion_threshold_stratiform_scale",
            3.0,
            1.0,
            1.0,
            10.0,
        ),
        PD.constrained_gaussian(
            "cloud_ice_specific_humidity_autoconversion_threshold",
            5e-6,
            1e-6,
            1e-6,
            1e-5,
        ),
        PD.constrained_gaussian("cloud_fraction_steepness_scale", 1.2, 0.3, 0.3, 3.0),
        PD.constrained_gaussian("fixed_cloud_ice_terminal_velocity", 0.1, 0.02, 0.0, 0.3),
        PD.constrained_gaussian(
            "prescribed_cloud_droplet_number_concentration",
            1e8,
            5e7,
            1e7,
            5e8,
        ),
    ]),
)
