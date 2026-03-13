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
    make_scalar_covariance_observation_vector(
        vars,
        sample_date_ranges;
        scalar = 1.0,
        use_latitude_weights = true,
        min_cosd_lat = 0.1,
    )

Make a scalar covariance matrix using `vars` for each sample corresponding to
the dates in `sample_date_ranges`.
"""
function make_scalar_covariance_observation_vector(
    vars,
    sample_date_ranges;
    scalar = 1.0,
    use_latitude_weights = true,
    min_cosd_lat = 0.1,
)
    obs_vec = map(sample_date_ranges) do sample_date_range
        start_date = first(sample_date_range)
        end_date = last(sample_date_range)
        @info "Using scalar covariance matrix with"
        @info "Scalar: $scalar"
        @info "Latitude weighting: $use_latitude_weights"
        @info "Min cosd lat: $min_cosd_lat"
        covar_estimator = ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
            scalar,
            use_latitude_weights,
            min_cosd_lat,
        )
        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            vars,
            start_date,
            end_date,
        )
    end
    return obs_vec
end

function make_svd_observation_vector(vars, sample_date_ranges)
    available_sample_dates = intersect(ClimaAnalysis.dates.(vars)...)
    min_year = minimum(Dates.year.(available_sample_dates))
    max_year = maximum(Dates.year.(available_sample_dates))

    obs_vec = map(sample_date_ranges) do sample_date_range
        start_date = first(sample_date_range)
        end_date = last(sample_date_range)
        # Note: This assumes that the data is monthly
        # Find the starting and ending dates for the entire years
        min_start_date = Date(min_year, month(start_date), day(start_date))
        max_start_date = Date(max_year, month(start_date), day(start_date))
        min_end_date = Date(min_year, month(end_date), day(end_date))
        max_end_date = Date(max_year, month(end_date), day(end_date))

        monthly_start_dates = collect(min_start_date:Dates.Year(1):max_start_date)
        monthly_end_dates = collect(min_end_date:Dates.Year(1):max_end_date)

        monthly_sample_date_ranges = [
            (monthly_start_date, monthly_end_date) for
            (monthly_start_date, monthly_end_date) in
            zip(monthly_start_dates, monthly_end_dates)
        ]
        @info "Samples used for $monthly_sample_date_ranges for generating SVDplusDCovariance matrix"
        covar_estimator = ClimaCalibrate.ObservationRecipe.SVDplusDCovariance(
            monthly_sample_date_ranges;
            model_error_scale = 0.05,
            regularization = ClimaCalibrate.ObservationRecipe.QuantileRegularization(0.01),
            use_latitude_weights = true,
            min_cosd_lat = 0.1,
            rank = 5,
        )
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
    data_loaders = [era5_pl_data_loader, ceres_data_loader, modis_data_loader]

    (; short_names) = CALIBRATE_CONFIG
    # Map short name to the corresponding data loader
    loader_registry = Dict()

    # Populate loader registry, mapping short names to loaders
    for short_name in short_names
        idx = findfirst(l -> short_name in l.available_vars, data_loaders)
        !isnothing(idx) && (loader_registry[short_name] = data_loaders[idx])
    end
    data_loader_non_unique_names =
        intersect(ClimaCoupler.CalibrationTools.available_vars.(data_loaders)...)
    any(x -> x in data_loader_non_unique_names, data_loader_non_unique_names) &&
        error("Data loader variable names are not unique: $data_loaders")

    Set(short_names) == keys(loader_registry) || error("Not all short names are being used")

    vars = map(short_names) do short_name
        data_loader = loader_registry[short_name]
        @info "Retrieving $(short_name) from $(typeof(data_loader))"
        var = get(data_loader, short_name)
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
    observation_vec = make_svd_observation_vector(
        vars,
        sample_date_ranges
    )

    # Save observation vector
    output_path = joinpath(pkgdir(ClimaCoupler), "experiments", "calibration", "amip")
    JLD2.save_object(joinpath(output_path, "observation_vec.jld2"), observation_vec)
    # Analyze noise covariance structure immediately upon construction
    # for (i, obs) in enumerate(observation_vec)
    #     noise_analysis = analyze_noise_covariance(obs)
    #     @info "Observation $i noise covariance" eigvalues = noise_analysis.eigvalues effective_rank =
    #         noise_analysis.effective_rank condition_number = noise_analysis.condition_number
    # end
    # Reconstruct the variables from the observation and show them for debugging
    for (i, obs) in enumerate(observation_vec)
        @info "Observation $i"
        @info ClimaCalibrate.ObservationRecipe.reconstruct_vars(observation_vec[i])
    end
end
