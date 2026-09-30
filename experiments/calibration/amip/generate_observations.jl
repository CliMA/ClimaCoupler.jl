import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCalibrate: ObservationRecipe, SampleBuilder
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
    joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip", "noise_model.jl"),
)

"""
    make_svdplusd_observation_vector(
        vars,
        sample_date_ranges;
        covariance_date_ranges = sample_date_ranges,
        beta = 0.05^2,
        rank = 2,
        sigma2 = 1e-6,
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
        FT = Float32,
    )

Make a vector of `EKP.Observation`s with an `SVDplusD` covariance matrix, one for each
sample corresponding to the dates in `sample_date_ranges`.

The covariance is the `Gamma` of `noise_model.jl`, estimated from one sample per date
range in `covariance_date_ranges`, so it is the interannual spread of the observation
across those dates. Give it enough date ranges to estimate that spread: `rank` modes need
appreciably more than `rank` samples, and identical date ranges produce a singular
covariance. The calibration targets in `sample_date_ranges` may repeat, and may be a
subset of the covariance dates, since the model only needs initial conditions for those.

`beta` accepts one value for all variables or one per variable, in the order of `vars`.
"""
function make_svdplusd_observation_vector(
    vars,
    sample_date_ranges;
    covariance_date_ranges = sample_date_ranges,
    beta = 0.05^2,
    rank = 2,
    sigma2 = 1e-6,
    use_latitude_weights = true,
    min_cosd_lat = 0.1,
    FT = Float32,
)
    @info "Using SVDplusD covariance matrix with" beta rank sigma2 use_latitude_weights min_cosd_lat
    covar_estimator =
        noise_covariance_estimator(; beta, rank, sigma2, use_latitude_weights, min_cosd_lat)

    covariance_samples =
        SampleBuilder.build_samples_by_times(vars, covariance_date_ranges; FT)
    @info "Built covariance samples" covariance_samples
    covar = ObservationRecipe.covariance(covar_estimator, covariance_samples)

    target_samples = SampleBuilder.build_samples_by_times(vars, sample_date_ranges; FT)
    @info "Built target samples" target_samples
    n_cov = size(SampleBuilder.get_samples(covariance_samples), 1)
    n_target = size(SampleBuilder.get_samples(target_samples), 1)
    n_cov == n_target || error(
        "The covariance samples ($n_cov values) and the targets ($n_target values) have different lengths",
    )

    # The same assembly as ObservationRecipe.observation, with the covariance supplied.
    obs_vec = map(1:SampleBuilder.num_samples(target_samples)) do i
        sample = collect(view(SampleBuilder.get_samples(target_samples), :, i))
        metadata = collect(view(SampleBuilder.get_metadata(target_samples), :, i))
        name = join(ClimaAnalysis.short_name.(metadata), ";")
        EKP.Observation(
            Dict(
                "samples" => sample,
                "covariances" => covar,
                "names" => name,
                "metadata" => metadata,
            ),
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
        get(source_data_loader, short_name)
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

    # Keep this in step with the zonal average in the other file.
    vars = zonal_average.(vars)

    # Give every sample the same NaN mask, or `build_samples_by_times` rejects a product
    # whose coverage varies by year.
    vars = ClimaAnalysis.propagate_nans.(vars; dims = ("time",))

    # Normalize data
    normalization_stats = Dict()
    compute_normalization!.(Ref(normalization_stats), vars)
    apply_normalization!.(Ref(normalization_stats), vars)
    (; output_dir) = CALIBRATE_CONFIG
    JLD2.save_object(NORMALIZATION_STATS_FP, normalization_stats)

    # Create observation vector
    (; sample_date_ranges) = CALIBRATE_CONFIG
    observation_vec = make_svdplusd_observation_vector(
        vars,
        sample_date_ranges;
        covariance_date_ranges,
        beta = 0.05^2,
        rank = 2,
        sigma2 = 1e-6,
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
    )

    # Save observation vector
    output_path = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip")
    JLD2.save_object(joinpath(output_path, "observation_vec.jld2"), observation_vec)

    # Reconstruct the variables from the observation and show them for debugging
    for (i, obs) in enumerate(observation_vec)
        @info "Observation $i"
        @info ObservationRecipe.reconstruct_vars(obs)
    end
end
