# One-year calibration on rsut, rlut, and pr against CERES and GPCP, December 2009 to November
# 2010. Each field enters the loss through `mean_pattern_seasonal_observables` as three
# observables that add up to its seasonal maps and do not overlap, so one cannot trade
# against another: its global annual mean (`annual_global_mean`), its annual map minus that
# mean (`annual_pattern`), and each season's map minus the annual map (`seasonal_anomaly`).
#
# Members start 2009-11-01 (one month of spinup) and end 2010-12-01.

const CALIBRATE_CONFIG = CalibrationTools.CalibrateConfig(;
    config_file = joinpath(
        pkgdir(ClimaCoupler),
        "config",
        "amip_configs",
        "amip_calibration.yml",
    ),
    # Each observable's per-entry error scale (`global_mean_sd` in physical units for the
    # 1-entry global mean, `annual_pattern_fraction` and `seasonal_cycle_fraction` of the
    # reduced observable's time mean for the two maps) is set from the baseline h_elem = 12
    # run's root-mean-square error (RMSE) against CERES and GPCP so every part starts near
    # one whitened misfit unit per entry (fraction = baseline_RMSE / observed_RMS):
    # - rsut: baseline global-mean bias is 12.6 W m^-2, and `baseline_RMSE / observed_RMS`
    #   on the reduced maps gives 0.42 for the annual pattern and 0.10 for the seasonal
    #   cycle.
    # - rlut: baseline global-mean bias is 4.8 W m^-2 (0.187 of the 25.5 W m^-2 observed
    #   longwave cloud effect), and 0.187 of the observed RMS for both maps.
    # - pr: baseline global-mean error scale is 0.85 mm/day (0.316 of the 2.7 mm/day
    #   observed global mean), 0.316 for the annual pattern, and 0.534 for the seasonal
    #   cycle.
    # `mean_pattern_seasonal_observables` sets `effective_entries = 29` on all three parts,
    # which normalizes each part by its entry count and weights it as 29 effective
    # observations in both the loss and `EffectiveEntriesController` (giving equal 1/3 block
    # weight to the global mean, annual pattern, and seasonal cycle, equivalent to a
    # `global_mean_sd / sqrt(29)` effective tolerance on the 1-entry global mean).
    observables = [
        mean_pattern_seasonal_observables(
            "rsut";
            global_mean_sd = 12.6,
            annual_pattern_fraction = 0.42,
            seasonal_cycle_fraction = 0.10,
        )
        mean_pattern_seasonal_observables(
            "rlut";
            global_mean_sd = 4.8,
            annual_pattern_fraction = 0.187,
            seasonal_cycle_fraction = 0.187,
        )
        mean_pattern_seasonal_observables(
            "pr";
            global_mean_sd = 0.85,
            annual_pattern_fraction = 0.316,
            seasonal_cycle_fraction = 0.534,
        )
    ],
    # A 13-month member takes about 10.7 h on one GPU. Each worker job requests 11 h (within
    # Derecho's 12 h queue limit), and `observation_map` retires workers that lack walltime
    # for another iteration so ClimaCalibrate 0.7 relaunches fresh workers for iterations 2
    # and 3.
    n_iterations = 3,
    worker_walltime = 11 * 60,
    # DJF 2009-10 to SON 2010. The noise is from the December-to-November years that CERES and
    # GPCP both cover.
    sample_date_ranges = [(DateTime(2009, 12, 1), DateTime(2010, 11, 1))],
    covariance_date_ranges = [
        (DateTime(y, 12, 1), DateTime(y + 1, 11, 1)) for y in 2000:2022
    ],
    extend = Month(1),
    spinup = Month(1),
    # On Derecho, make this a symlink to scratch: member output fills the home quota.
    output_dir = joinpath(pkgdir(ClimaCoupler), "amip_calibration_year2010_rsut_rlut_pr"),
    # Prior mean, standard deviation, and bounds in physical units. The sgs_variance
    # parameters' released defaults (1.24, and 1.08 m/s) are fitted for h_elem 16 to 60 and
    # make rsut far too bright at h_elem 12, so their prior means use the coarse-grid
    # values.
    priors = PD.combine_distributions([
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
            3.6e-6,
            1e-6,
            5e-5,
        ),
        PD.constrained_gaussian(
            "ice_cloud_effective_radius",
            2.5e-5,
            7.5e-6,
            1e-5,
            5e-5,
        ),
        PD.constrained_gaussian("entr_coeff", 0.30, 0.09, 0.0, 0.55),
        PD.constrained_gaussian(
            "rain_autoconversion_timescale",
            1200.0,
            430.0,
            100.0,
            10000.0,
        ),
        PD.constrained_gaussian(
            "sgs_variance_horizontal_scale_factor",
            0.219,
            0.27,
            0.0,
            1.8,
        ),
        PD.constrained_gaussian(
            "sgs_variance_vertical_velocity_floor",
            0.2,
            0.2,
            0.0,
            2.5,
        ),
    ]),
)
