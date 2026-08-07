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
- `decorrelation_length`: if set (meters), the diagonal floor of the covariance
  is replaced with a spatially correlated one, `sqrt(D) * K * sqrt(D)` with
  `K_ij = exp(-d_ij / L)`. Use this for 2-D (lon×lat) observation vectors:
  a diagonal floor treats every grid point as independent, over-informs the
  inverse, and collapses the ensemble. See correlated_noise.jl.
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
    # CALIPSO/CloudSat provides level-resolved cloud fraction `cl`.
    calipso_data_loader = CalibrationTools.CalipsoDataLoader()
    # GPCP provides monthly precipitation `pr` (mm/day; converted to match the
    # model's kg m^-2 s^-1 on the sim side in observation_map.jl).
    gpcp_data_loader = CalibrationTools.GPCPDataLoader()
    # Both MODIS and MAC provide `lwp`, so disambiguate to get `lwp` from MAC.
    # MODIS is kept for its ice water path (`clivi`).
    data_loader = CalibrationTools.CompositeDataLoader(
        era5_pl_data_loader,
        ceres_data_loader,
        modis_data_loader,
        mac_data_loader,
        calipso_data_loader,
        gpcp_data_loader;
        varname_to_loader = Dict(
            "lwp" => mac_data_loader,
            "cl" => calipso_data_loader,
            "pr" => gpcp_data_loader,
        ),
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
    vars = select_altitude_levels.(vars, Ref(ALTITUDE_LEVELS))
    lonlat_regridder = get_lonlat_regridder(config_file)
    vars = lonlat_regridder.(vars)
    lat_left, lat_right = lat_window()
    vars = apply_lat_window.(vars, lat_left, lat_right)

    # Harmonize each variable's NaN pattern across the covariance dates so all
    # SVDplusD samples have equal length. Satellite lwp coverage varies by year;
    # without this the interannual samples drop different numbers of points and
    # SVDplusDCovariance errors ("Length of all the samples are not the same").
    #
    # NOTE: this runs BEFORE the zonal mean (it used to run after). Harmonizing the
    # 2-D field means every date's zonal mean is taken over an identical set of
    # grid points, so year-to-year changes in satellite coverage can no longer leak
    # into the interannual spread the covariance is estimated from. It also still
    # guarantees equal sample lengths, since a latitude is now dropped only when all
    # of its longitudes are missing — consistently for every date.
    foreach(v -> harmonize_nan_mask_over_dates!(v, COVARIANCE_DATE_RANGES), vars)

    # Save where each observation actually has data, so the simulation can be
    # restricted to the SAME points before its own zonal mean. Without this the
    # observed zonal mean is ocean-only (MAC lwp is NaN over land, ~54% of points)
    # while the simulated one covers all longitudes — two different spatial
    # samples, which for lwp flips the sign of the area-weighted bias. Must be
    # written before zonal_average, while longitude still exists.
    coverage_masks = Dict{String, Any}(
        ClimaAnalysis.short_name(v) => coverage_mask(v, COVARIANCE_DATE_RANGES) for
        v in vars
    )

    # Reduce the spatial dimensions to approximately independent constraints.
    # Native grid points are strongly correlated, and treating them as
    # independent over-informs the inverse and collapses the ensemble. The
    # default is a zonal mean. A config that defines COARSEN_FACTOR keeps two
    # spatial dimensions on a block-averaged grid instead (see coarsen_lonlat).
    vars = reduce_spatial.(vars)

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
    # model_error_scale adds `(model_error_scale * mean(field))^2` to the diagonal
    # as a structural-model-error floor. The floor states how much mismatch a
    # perfect parameter set is expected to leave. 0.2 (20% of the field mean) is
    # comparable to the internal variability of lwp and pr. An observable with a
    # larger irreducible bias needs a larger floor, otherwise the calibration
    # distorts reachable parameters to chase it.
    #
    # A config may define OBS_NOISE_GROUPS to give variable groups their own
    # floors. Each group becomes its own covariance block and the per-sample
    # observations are combined with EKP.combine_observations, so the flattened
    # layout is group order then variable order within the group.
    noise_groups = @isdefined(OBS_NOISE_GROUPS) ? OBS_NOISE_GROUPS :
                   [(short_names = short_names, model_error_scale = 0.2)]

    grouped_obs_vecs = map(noise_groups) do group
        group_vars = filter(
            v -> ClimaAnalysis.short_name(v) in group.short_names, collect(vars),
        )
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
            # The regularization quantile needs at least 1/qtl values per
            # variable, so 0.05 requires 20. A reduced observation can have
            # far fewer: a zonal mean over a 20 degree band leaves 8. Set
            # OBS_REGULARIZATION_QUANTILE above 1/n_values in that case.
            regularization = ClimaCalibrate.ObservationRecipe.QuantileRegularization(
                @isdefined(OBS_REGULARIZATION_QUANTILE) ?
                OBS_REGULARIZATION_QUANTILE : 0.05,
            ),
        )
    end
    observation_vec = [
        length(grouped_obs_vecs) == 1 ? grouped_obs_vecs[1][k] :
        EKP.combine_observations([gv[k] for gv in grouped_obs_vecs]) for
        k in eachindex(sample_date_ranges)
    ]

    # Save observation vector into the config's output_dir so that different
    # calibration setups (e.g. lwp+cl vs LWP-only) keep independent observations.
    JLD2.save_object(joinpath(output_dir, "observation_vec.jld2"), observation_vec)

    # Coverage masks live next to the observation vector; preprocess_sim_vars reads
    # them so the simulation is sampled at the same points as the observation.
    JLD2.save_object(coverage_mask_path(output_dir), coverage_masks)
    for (name, (dims, mask)) in coverage_masks
        @info "Coverage mask for $name" dims missing_fraction =
            round(count(mask) / length(mask); digits = 3)
    end

    # Reconstruct the variables from the observation and show them for debugging
    for (i, obs) in enumerate(observation_vec)
        @info "Observation $i"
        @info ClimaCalibrate.ObservationRecipe.reconstruct_vars(observation_vec[i])
    end
end
