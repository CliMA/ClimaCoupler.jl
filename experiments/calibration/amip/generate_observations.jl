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

"""
    make_scalar_covariance_observation_vector(
        vars,
        sample_date_ranges;
        scalar = 1.0,
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
        FT = Float32,
    )

Make a vector of `EKP.Observation`s with a scalar covariance matrix, one for
each sample corresponding to the dates in `sample_date_ranges`.

The `OutputVar`s in `vars` are windowed by the date ranges in
`sample_date_ranges` to build the samples, and each sample in turn is used as
the observation. The matrix of samples has element type `FT`.
"""
function make_scalar_covariance_observation_vector(
    vars,
    sample_date_ranges;
    scalar = 1.0,
    use_latitude_weights = true,
    min_cosd_lat = 0.1,
    FT = Float32,
)
    @info "Using scalar covariance matrix with"
    @info "Scalar: $scalar"
    @info "Latitude weighting: $use_latitude_weights"
    @info "Min cosd lat: $min_cosd_lat"
    covar_estimator =
        ObservationRecipe.ScalarCovariance(; scalar, use_latitude_weights, min_cosd_lat)

    # Each date range becomes one sample (one column of the sample collection)
    sample_collection = SampleBuilder.build_samples_by_times(vars, sample_date_ranges; FT)
    @info "Built samples" sample_collection

    obs_vec = map(1:SampleBuilder.num_samples(sample_collection)) do i
        ObservationRecipe.observation(covar_estimator, sample_collection, i)
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

    # Normalize data
    normalization_stats = Dict()
    compute_normalization!.(Ref(normalization_stats), vars)
    apply_normalization!.(Ref(normalization_stats), vars)
    (; output_dir) = CALIBRATE_CONFIG
    JLD2.save_object(NORMALIZATION_STATS_FP, normalization_stats)

    # Create observation vector
    (; sample_date_ranges) = CALIBRATE_CONFIG
    observation_vec = make_scalar_covariance_observation_vector(
        vars,
        sample_date_ranges;
        scalar = 1.0,
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
