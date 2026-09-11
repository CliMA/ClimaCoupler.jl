import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCalibrate: ObservationRecipe
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import EnsembleKalmanProcesses as EKP
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
include(
    joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "calibration",
        "amip",
        "correlated_noise.jl",
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

Make one `EKP.Observation` per entry of `sample_date_ranges`, all sharing one
`SVDplusDCovariance` estimated from the interannual spread over
`covariance_date_ranges`.

`sample_date_ranges` defines the date range of the sample data, while
`covariance_date_ranges` defines the date range used to estimate the
covariance.. Every sample entry must be contained in the covariance entries.

Keyword arguments:
- `model_error_scale`: adds `(model_error_scale * mean_sample)^2` to the
  covariance diagonal, the error a perfect parameter set can not reduce.
- `regularization`: conditioning floor, quantile-based so it scales per variable
  in physical units.
- `rank`: SVD rank, `nothing` infers it (at most n_covariance_dates - 1).
- `decorrelation_length`: if set (meters), splits the diagonal floor D into
  `f * D + (1 - f) * sqrt(D) * K * sqrt(D)` with `K_ij = exp(-d_ij / L)`, so a
  patch of correlated points counts as one constraint rather than one per point.
  See correlated_noise.jl.
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
    decorrelation_length = nothing,
)
    @info "Using SVDplusD data-informed covariance with"
    @info "Model error scale: $model_error_scale"
    @info "Latitude weighting: $use_latitude_weights"
    @info "Covariance estimated from $(length(covariance_date_ranges)) dates"
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
    isnothing(decorrelation_length) && return obs_vec
    @info "Correlating the noise floor with decorrelation length $(decorrelation_length) m"
    return apply_floor_correlation(obs_vec, decorrelation_length)
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
    # SVDplusD samples have equal length.
    foreach(v -> harmonize_nan_mask_over_dates!(v, COVARIANCE_DATE_RANGES), vars)

    # Reduction
    vars = reduce_spatial.(vars, COARSEN_FACTOR)

    # No normalization: the SVDplusD covariance carries physical scales.
    (; output_dir) = CALIBRATE_CONFIG
    isfile(NORMALIZATION_STATS_FP) && rm(NORMALIZATION_STATS_FP)

    (; sample_date_ranges) = CALIBRATE_CONFIG

    # A config may define OBS_NOISE_GROUPS to give variable groups their own
    # model error floors. A floor should reflect the bias the model keeps at the
    # best parameter values, which differs per observable.
    noise_groups =
        @isdefined(OBS_NOISE_GROUPS) ? OBS_NOISE_GROUPS :
        [(short_names = short_names, model_error_scale = 0.2)]

    grouped_obs_vecs = map(noise_groups) do group
        group_vars =
            filter(v -> ClimaAnalysis.short_name(v) in group.short_names, collect(vars))
        isempty(group_vars) && error(
            "OBS_NOISE_GROUPS entry $(group.short_names) matches no loaded variables",
        )
        make_data_informed_observation_vector(
            group_vars,
            sample_date_ranges,
            COVARIANCE_DATE_RANGES;
            model_error_scale = group.model_error_scale,
            use_latitude_weights = true,
            min_cosd_lat = 0.1,
            decorrelation_length = hasproperty(group, :decorrelation_length) ?
                                   group.decorrelation_length : nothing,
        )
    end
    observation_vec = [
        length(grouped_obs_vecs) == 1 ? grouped_obs_vecs[1][k] :
        EKP.combine_observations([gv[k] for gv in grouped_obs_vecs]) for
        k in eachindex(sample_date_ranges)
    ]

    # Saved observations and coverage masks 
    coverage_masks = Dict{String, Any}(
        ClimaAnalysis.short_name(v) => coverage_mask(v, COVARIANCE_DATE_RANGES) for
        v in vars
    )
    JLD2.save_object(joinpath(output_dir, "observation_vec.jld2"), observation_vec)
    JLD2.save_object(coverage_mask_path(output_dir), coverage_masks)
    for (name, (dims, mask)) in coverage_masks
        @info "Coverage mask for $name" dims missing_fraction =
            round(count(mask) / length(mask); digits = 3)
    end

    # Reconstruct the variables from the observation and show them for debugging
    for (i, obs) in enumerate(observation_vec)
        @info "Observation $i"
        @info ObservationRecipe.reconstruct_vars(obs)
    end
end
