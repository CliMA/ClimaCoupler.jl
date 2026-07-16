import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import JLD2

include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "run_calibration.jl",
    ),
)
include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "preprocessing.jl",
    ),
)

"""
    make_data_informed_observation_vector(
        vars,
        sample_date_ranges,
        covariance_date_ranges;
        model_error_scale = 0.1,
        regularization = ClimaCalibrate.ObservationRecipe.QuantileRegularization(0.05),
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
        rank = nothing,
    )

Make one `EKP.Observation` per calibration target in `sample_date_ranges`, all
sharing a single data-informed `SVDplusDCovariance` estimated from the
(independent) `covariance_date_ranges`.

The two date lists are deliberately separate: `sample_date_ranges` defines *what
is being calibrated against* (the target samples), while `covariance_date_ranges`
defines the realizations used to *estimate the noise*. This lets us broaden the
covariance sample (e.g. the same month across many years) without changing the
calibration target. Each `sample_date_ranges` entry must be contained in
`covariance_date_ranges` (SVDplusD requires the sampled date to be one of the
covariance dates).

Unlike `ScalarCovariance(1.0)`, which asserts a uniform, physically meaningless
unit noise, `SVDplusDCovariance` estimates the observational + internal
variability directly from the interannual spread of `vars`. This is what tells
EKP how much of the model–obs mismatch is signal versus noise; with the scalar
covariance the parameter signal was drowned and the error stayed flat.

Keyword arguments
=================
- `model_error_scale`: structural model-error term added to the diagonal, as
  `(model_error_scale * mean_sample)^2`. This sets the irreducible error floor a
  perfect parameter set is expected to leave (e.g. 0.1 ⇒ 10% of the field mean).
- `regularization`: floor added to the covariance for conditioning. A
  `QuantileRegularization` is used because it scales per-variable, which matters
  now that variables are in physical (unnormalized) units of very different
  magnitudes (ta ~ 100s K vs lwp ~ 0.1 kg m⁻²).
- `rank`: SVD rank; `nothing` infers it (≤ number of covariance dates − 1).
"""
function make_data_informed_observation_vector(
    vars,
    sample_date_ranges,
    covariance_date_ranges;
    model_error_scale = 0.1,
    regularization = ClimaCalibrate.ObservationRecipe.QuantileRegularization(0.05),
    use_latitude_weights = true,
    min_cosd_lat = 0.1,
    rank = nothing,
)
    @info "Using SVDplusD data-informed covariance with"
    @info "Model error scale: $model_error_scale"
    @info "Latitude weighting: $use_latitude_weights"
    @info "Covariance estimated from $(length(covariance_date_ranges)) dates"
    # The covariance is built once from the interannual spread across the
    # covariance dates and reused for every calibration target.
    covar_estimator = ClimaCalibrate.ObservationRecipe.SVDplusDCovariance(
        covariance_date_ranges;
        model_error_scale,
        regularization,
        use_latitude_weights,
        min_cosd_lat,
        rank,
    )
    obs_vec = map(sample_date_ranges) do sample_date_range
        start_date = first(sample_date_range)
        end_date = last(sample_date_range)
        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            vars,
            start_date,
            end_date,
        )
    end
    return obs_vec
end

if abspath(PROGRAM_FILE) == @__FILE__
    # Prevent MPI from being used which is not needed for generating
    # observations
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"

    # Create data loaders (constructing these are relatively cheap)
    era5_pl_data_loader = CalibrationTools.ERA5PressureLevelDataLoader()
    ceres_data_loader = CalibrationTools.CERESDataLoader()
    modis_data_loader = CalibrationTools.ModisDataLoader()
    mac_data_loader = CalibrationTools.MACDataLoader()
    # Both MODIS and MAC provide `lwp`, so disambiguate to get `lwp` from MAC.
    # MODIS is kept for its ice water path (`clivi`).
    data_loader = CalibrationTools.CompositeDataLoader(
        era5_pl_data_loader,
        ceres_data_loader,
        modis_data_loader,
        mac_data_loader;
        varname_to_loader = Dict("lwp" => mac_data_loader),
    )

    (; short_names) = CALIBRATE_CONFIG

    vars = map(short_names) do short_name
        source_data_loader = CalibrationTools.find_source_loader(data_loader, short_name)
        @info "Retrieving $(short_name) from $(typeof(source_data_loader))"
        var = get(source_data_loader, short_name)
    end

    # For now, we apply the preprocessing to all the variables if possible
    # If the preprocessing does not apply, then it is a no-op.
    # In the future, if we want to do specific preprocessing, this needs to
    # change
    vars = select_pressure_levels.(vars, Ref(PRESSURE_LEVELS))
    lonlat_regridder = get_lonlat_regridder(config_file)
    vars = lonlat_regridder.(vars)
    lat_left = -90
    lat_right = 90
    vars = apply_lat_window.(vars, lat_left, lat_right)

    # Harmonize each variable's NaN pattern across the covariance dates so all
    # SVDplusD samples have equal length. Satellite lwp coverage varies by year;
    # without this the interannual samples drop different numbers of points and
    # SVDplusDCovariance errors ("Length of all the samples are not the same").
    foreach(v -> harmonize_nan_mask_over_dates!(v, COVARIANCE_DATE_RANGES), vars)

    # NOTE: Normalization is intentionally NOT applied. The SVDplusD covariance
    # below carries each variable's physical scale, and normalization is
    # unsupported with it. `preprocess_sim_vars` skips normalization when
    # NORMALIZATION_STATS_FP is absent, so the sim side stays consistent as long
    # as no stale normalization_stats.jld2 is left in output_dir.
    (; output_dir) = CALIBRATE_CONFIG
    isfile(NORMALIZATION_STATS_FP) && rm(NORMALIZATION_STATS_FP)

    # Create observation vector: calibration targets from sample_date_ranges, with
    # a data-informed covariance estimated from the independent COVARIANCE_DATE_RANGES.
    (; sample_date_ranges) = CALIBRATE_CONFIG
    observation_vec = make_data_informed_observation_vector(
        vars,
        sample_date_ranges,
        COVARIANCE_DATE_RANGES;
        model_error_scale = 0.1,
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
    )

    # Save observation vector
    output_path = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip")
    JLD2.save_object(joinpath(output_path, "observation_vec.jld2"), observation_vec)

    # Reconstruct the variables from the observation and show them for debugging
    for (i, obs) in enumerate(observation_vec)
        @info "Observation $i"
        @info ClimaCalibrate.ObservationRecipe.reconstruct_vars(observation_vec[i])
    end
end
