import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCoupler
import ClimaCoupler: CalibrationTools
import Dates
import JLD2
import NCDatasets
import Statistics
using OrderedCollections: OrderedDict

# Access CalibrateConfig (this also includes observation_utils.jl via observation_map.jl)
include(joinpath(@__DIR__, "run_calibration.jl"))

# Include preprocessing utils
include(joinpath(@__DIR__, "preprocessing_utils.jl"))

# Mapping from ERA5 file prefixes to short names used in calibration
const ERA5_FILE_PREFIX_TO_SHORT_NAME = Dict(
    "2m_temperature" => "tas",
    "mean_sea_level_pressure" => "mslp",
    "total_precipitation" => "pr",
    "sea_surface_temperature" => "ts",
    "surface_pressure" => "sp",
    "mean_top_upward_short_wave_radiation_flux" => "rsut",
    "mean_top_upward_long_wave_radiation_flux" => "rlut",
)

# Mapping from ERA5 NetCDF variable names to short names
const ERA5_VARNAME_TO_SHORT_NAME = Dict(
    "t2m" => "tas",
    "msl" => "mslp",
    "tp" => "pr",
    "sst" => "ts",
    "sp" => "sp",
    "mean_top_upward_short_wave_radiation_flux" => "rsut",
    "mean_top_upward_long_wave_radiation_flux" => "rlut",
)

"""
    get_weekly_filename(obs_dir, short_name, start_date, end_date)

Construct the filename for a weekly ERA5 file based on variable and date range.
Files follow the pattern: {era5_var}_weekly_mean_{start}_{end}.nc
"""
function get_weekly_filename(obs_dir, short_name, start_date, end_date)
    era5_prefix = short_name_to_era5_prefix(short_name)
    start_str = Dates.format(start_date, "yyyymmdd")
    end_str = Dates.format(end_date, "yyyymmdd")
    filename = "$(era5_prefix)_weekly_mean_$(start_str)_$(end_str).nc"
    return joinpath(obs_dir, filename)
end

"""
    get_daily_filename(obs_dir, short_name, target_date)

Find a daily-mean ERA5 file that contains the target_date.
Files follow the pattern: {era5_var}_daily-mean_{start}_{end}.nc
Returns the filepath and the (start_date, end_date) range of the file.
"""
function get_daily_filename(obs_dir, short_name, target_date)
    era5_prefix = short_name_to_era5_prefix(short_name)
    parent_dir = dirname(obs_dir)
    
    pattern = Regex("$(era5_prefix)_daily-mean_(\\d{8})_(\\d{8})\\.nc")
    
    for filename in readdir(parent_dir)
        m = match(pattern, filename)
        if !isnothing(m)
            file_start = Dates.DateTime(m.captures[1], "yyyymmdd")
            file_end = Dates.DateTime(m.captures[2], "yyyymmdd")
            if file_start <= target_date <= file_end
                return joinpath(parent_dir, filename), (file_start, file_end)
            end
        end
    end
    
    error(
        "No daily-mean file found containing date $target_date for variable $short_name in $parent_dir",
    )
end

"""
    short_name_to_era5_prefix(short_name)

Convert calibration short_name to ERA5 file prefix.
"""
function short_name_to_era5_prefix(short_name)
    for (prefix, sn) in ERA5_FILE_PREFIX_TO_SHORT_NAME
        if sn == short_name
            return prefix
        end
    end
    error("Unknown short name: $short_name")
end

"""
    load_daily_var(filepath, short_name, target_date)

Load a specific day from a daily-mean ERA5 file.
Uses ClimaAnalysis to load and select the specific date.
Note: Returns 2D data (lon, lat) - time dimension is added in preprocess_vars.
"""
function load_daily_var(filepath, short_name, target_date)
    era5_varname = short_name_to_era5_varname(short_name)
    @info "Loading daily $short_name for $target_date from $filepath"
    
    var = ClimaAnalysis.OutputVar(filepath, era5_varname)

    time_name = ClimaAnalysis.time_name(var)
    var = ClimaAnalysis.slice(var; (Symbol(time_name) => target_date,)...)

    if !issorted(ClimaAnalysis.latitudes(var))
        var = ClimaAnalysis.reverse_dim(var, ClimaAnalysis.latitude_name(var))
    end
    new_attribs = copy(var.attributes)
    new_attribs["short_name"] = short_name
    var = ClimaAnalysis.OutputVar(new_attribs, var.dims, var.dim_attributes, var.data)
    
    return var
end

"""
    short_name_to_era5_varname(short_name)

Convert calibration short_name to ERA5 NetCDF variable name.
"""
function short_name_to_era5_varname(short_name)
    for (varname, sn) in ERA5_VARNAME_TO_SHORT_NAME
        if sn == short_name
            return varname
        end
    end
    error("Unknown short name: $short_name")
end

"""
    load_weekly_var(filepath, short_name, start_date, end_date)

Load a weekly ERA5 variable from a NetCDF file and set appropriate attributes.
"""
function load_weekly_var(filepath, short_name, start_date, end_date)
    era5_varname = short_name_to_era5_varname(short_name)
    @info "Loading weekly $short_name from $filepath"
    var = OutputVar(filepath, era5_varname)
    
    var.attributes["short_name"] = short_name
    var.attributes["start_date"] = start_date
    
    if !issorted(ClimaAnalysis.latitudes(var))
        var = ClimaAnalysis.reverse_dim(var, ClimaAnalysis.latitude_name(var))
    end

    # Note: We do NOT add time dimension here - it's added in preprocess_vars after resampling
    # This avoids issues with ClimaAnalysis.resampled_as which doesn't handle 3D data
    
    return var
end

# """
#     add_time_dimension(var, start_date)

# Add a time dimension to a 2D OutputVar for ObservationRecipe compatibility.
# Time is stored as seconds since start_date, so ClimaAnalysis.dates() can 
# convert it back to DateTime using the start_date attribute.
# """
# function add_time_dimension(var, start_date)
#     if haskey(var.dims, "time")
#         return var  # Already has time dimension
#     end
    
#     # Time value = 0.0 means "at start_date"
#     # ClimaAnalysis.dates() will compute: start_date + time_value seconds
#     time_val = 0.0
    
#     # Create new dims with time, preserving order (lon, lat, time)
#     # and converting to Float64 vectors (ClimaCalibrate expects plain Float64)
#     new_dims = OrderedDict{String, Vector{Float64}}()
#     for k in keys(var.dims)
#         new_dims[k] = Float64.(collect(skipmissing(var.dims[k])))
#     end
#     new_dims["time"] = [time_val]
    
#     # Create new dim_attributes with time
#     # ClimaCalibrate expects "s" as the unit string for seconds
#     new_dim_attribs = copy(var.dim_attributes)
#     new_dim_attribs["time"] = Dict{String, Any}(
#         "units" => "s",
#         "calendar" => "standard",
#     )
    
#     # Reshape data to add time dimension, ensure Float64
#     new_data = reshape(Float64.(var.data), size(var.data)..., 1)
    
#     return ClimaAnalysis.OutputVar(
#         var.attributes,
#         new_dims,
#         new_dim_attribs,
#         new_data,
#     )
# end


##### CERES Data Loading Functions

"""
    load_ceres_var(short_name, config_file)

Load a variable from CERES monthly data for the month containing start_date.
Returns a 2D OutputVar (lon, lat) representing the monthly mean.

Note: CERES data is monthly, so we use the monthly value for the month
containing the calibration period. Time dimension is added later in preprocess_vars.
"""
function load_ceres_var(short_name, config_file)
    loader = CalibrationTools.CERESDataLoader()
    
    # Load the full time series
    var = Base.get(loader, short_name)
    resampler = resampled_lonlat(config_file)
    var = resampler(var)

    @info "Loaded CERES $short_name"
    
    return var
end

"""
    is_ceres_variable(short_name, ceres_variables)

Check if a variable should be loaded from CERES instead of ERA5.
"""
is_ceres_variable(short_name, ceres_variables) = short_name in ceres_variables


##### ERA5 Pressure-Level Data Loading Functions

"""
    load_era5_pressure_level_var(short_name)

Load a variable from ERA5 pressure-level data for the month containing start_date.
Handles parsing of pressure-level variable names (e.g., "ta_850hPa") and selects
the appropriate pressure level from the 3D data.

Returns a 2D OutputVar (lon, lat) representing the monthly mean at that pressure level.
"""
function load_era5_pressure_level_var(short_name, config_file)
    # Parse the short_name to get base variable and pressure level
    base_name, pressure_hPa = parse_pressure_level_variable(short_name)

    loader = CalibrationTools.ERA5PressureLevelDataLoader()
    # (lon, lat, pressure_level, time)
    var = Base.get(loader, String(base_name))

    # ERA5 artifact uses "pressure_level" dimension in Pa, convert hPa to Pa for slicing
    pressure_Pa = pressure_hPa * 100.0
    var = ClimaAnalysis.slice(var, pressure_level = pressure_Pa)

    resampler = resampled_lonlat(config_file)
    var = resampler(var)

    # Set the short_name to include pressure level (e.g., "ta_850hPa")
    new_attribs = copy(var.attributes)
    new_attribs["short_name"] = short_name
    var = ClimaAnalysis.OutputVar(new_attribs, var.dims, var.dim_attributes, var.data)
    
    @info "Loaded ERA5 $short_name ($(base_name) at $(Int(pressure_hPa)) hPa)"
    
    return var
end

"""
    load_vars(obs_dir, short_names, sample_date_ranges; ceres_variables = String[])

Load observation data for the specified short_names and date ranges.
- Pressure-level variables (e.g., "ta_850hPa") are loaded from ERA5 pressure-level artifact
- Variables in `ceres_variables` are loaded from CERES monthly data
- Other variables are loaded from ERA5 (daily or weekly files)

Returns a Dict mapping (short_name, date_range) -> OutputVar
"""
function load_vars(obs_dir, short_names, config_file; ceres_variables = String[])
    vars =[]

    pressure_level_vars = filter(is_pressure_level_variable, short_names)
    ceres_vars_in_names = filter(sn -> is_ceres_variable(sn, ceres_variables), short_names)
    # TODO: simplify these filters
    era5_surface_vars = filter(
        sn -> !is_pressure_level_variable(sn) && !is_ceres_variable(sn, ceres_variables),
        short_names,
    )

    # Create loaders once — artifact-backed loaders return full time series regardless of date
    era5_pl_loader = isempty(pressure_level_vars) ? nothing : CalibrationTools.ERA5PressureLevelDataLoader()
    ceres_loader = isempty(ceres_vars_in_names) ? nothing : CalibrationTools.CERESDataLoader()

    # Load ERA5 pressure-level vars once per variable
    for short_name in pressure_level_vars
        var = load_era5_pressure_level_var(short_name, config_file)
        push!(vars, var)
    end

    # Load CERES vars once per variable
    for short_name in ceres_vars_in_names
        var = load_ceres_var(short_name, config_file)
        for dr in sample_date_ranges
            push!(vars, var)
        end
    end

    # Load ERA5 surface vars from weekly/daily files — these are date-specific.
        for short_name in era5_surface_vars
            try
                if is_single_day
                    filepath, _ = get_daily_filename(obs_dir, short_name, start_date)
                    var = load_daily_var(filepath, short_name, start_date)
                else
                    # Load from ERA5 weekly files
                    filepath =
                        get_weekly_filename(obs_dir, short_name, start_date, end_date)
                    if !isfile(filepath)
                        @warn "Weekly file not found: $filepath"
                        continue
                    end
                    var = load_weekly_var(filepath, short_name, start_date, end_date)
                end
                push!(vars, var)
            catch e
                @warn "Failed to load $short_name for $start_date - $end_date: $e"
                rethrow(e)
            end
        end
    return vars
end

"""
    preprocess_vars(vars_by_date, config_file)

Preprocess each OutputVar by resampling to the model grid, setting units,
and adding time dimension for ObservationRecipe compatibility.

Precipitation (pr) requires special handling:
ERA5 'tp' is in meters/hour, we convert to mm/day using ERA5_PR_CONVERSION.
"""
function preprocess_vars(vars)
    # processed = Dict{Tuple{String, NTuple{2, Dates.DateTime}}, OutputVar}()
    
    # for (key, var) in vars_by_date
    #     short_name = key[1]
    #     date_range = key[2]
    #     start_date = date_range[1]

    #     if haskey(var_units, short_name)
    #         var = ClimaAnalysis.set_units(var, get_var_units(short_name))
    #     end
    #     var.attributes["short_name"] = short_name
    #     var.attributes["start_date"] = string(start_date)

    #     processed[key] = var
    # end
    processed = []
    for var in vars
        if haskey(var_units, ClimaAnalysis.short_name(var))
            # var = ClimaAnalysis.set_units(var, get_var_units(short_name))
            # There isn't a set_units! function yet
            var.attributes["units"] = get_var_units(ClimaAnalysis.short_name(var))
        end
        push!(processed, var)
    end
    return processed
end

"""
    make_observation_vector(vars_by_date, sample_date_ranges, short_names)

Create the observation vector for EKP using ObservationRecipe.ScalarCovariance,
matching the subseasonal pipeline pattern.

ObservationRecipe handles latitude weighting and covariance estimation internally,
replacing the custom normalization and land masking that was previously done manually.

Note: For weekly averages, each OutputVar has a single time point at start_date.
We pass start_date for both start and end to ObservationRecipe since the data
represents the entire week as a single snapshot.
"""
function make_observation_vector(vars, sample_date_ranges)
    include(joinpath(@__DIR__, "calibration_setup.jl"))

    available_sample_dates = intersect(ClimaAnalysis.dates.(vars)...)
    min_year = minimum(Dates.year.(available_sample_dates))
    max_year = maximum(Dates.year.(available_sample_dates))
    # covar_estimator = ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
    #     scalar = Float64(CALIBRATION_NOISE_SCALAR),
    #     use_latitude_weights = true,
    # )

    obs_vec = map(sample_date_ranges) do sample_date_range
        start_date = first(sample_date_range)
        min_date = Date(min_year, month(start_date), day(start_date))
        max_date = Date(max_year, month(start_date), day(start_date))
        monthly_sample_dates = collect(min_date:Dates.Year(1):max_date)
        monthly_sample_date_ranges = [(monthly_date, monthly_date) for monthly_date in monthly_sample_dates]
        covar_estimator = ClimaCalibrate.ObservationRecipe.SVDplusDCovariance(
                   monthly_sample_date_ranges;
                   model_error_scale = 0.05,
                   regularization = 0.1,
                   use_latitude_weights = true,
                   min_cosd_lat = 0.1,
        )
        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            vars,
            start_date,
            start_date,
        )
    end
    return obs_vec
end

"""
    compute_normalization_stats(vars_by_date, short_names)

Compute global mean and std for each variable across all date ranges.
Returns a Dict mapping short_name -> (mean, std).
Uses latitude-weighted averaging for physically meaningful statistics.
"""
function compute_normalization_stats(vars_by_date, short_names)
    norm_stats = Dict{String, Tuple{Float64, Float64}}()
    
    for short_name in short_names
        # Collect all data for this variable across date ranges
        all_data = Float64[]
        all_weights = Float64[]
        
        for (key, var) in vars_by_date
            if key[1] == short_name
                lats = ClimaAnalysis.latitudes(var)
                # Compute latitude weights (cosine weighting)
                lat_weights = cosd.(lats)
                
                # Flatten spatial data and replicate weights
                data = vec(var.data)
                valid_mask = .!isnan.(data)
                
                # Create weight array matching data shape
                nlat = length(lats)
                nlon = size(var.data, 1)
                weights_2d = repeat(lat_weights', nlon, 1)
                weights_flat = vec(weights_2d)
                
                append!(all_data, data[valid_mask])
                append!(all_weights, weights_flat[valid_mask])
            end
        end
        
        if isempty(all_data)
            @warn "No data found for $short_name, using default normalization (0, 1)"
            norm_stats[short_name] = (0.0, 1.0)
            continue
        end
        
        # Compute weighted mean and std
        total_weight = sum(all_weights)
        weighted_mean = sum(all_data .* all_weights) / total_weight
        weighted_var = sum(all_weights .* (all_data .- weighted_mean).^2) / total_weight
        weighted_std = sqrt(weighted_var)
        
        # Avoid division by zero
        if weighted_std < 1e-10
            @warn "$short_name has near-zero std, using 1.0"
            weighted_std = 1.0
        end
        
        norm_stats[short_name] = (weighted_mean, weighted_std)
        @info "Normalization stats for $short_name: mean=$(round(weighted_mean, sigdigits=5)), std=$(round(weighted_std, sigdigits=5))"
    end
    
    return norm_stats
end

"""
    normalize_var(var, norm_stats)

Normalize an OutputVar using precomputed normalization statistics.
Returns a new OutputVar with normalized data: (data - mean) / std
"""
function normalize_var(var, norm_stats)
    short_name = var.attributes["short_name"]
    if !haskey(norm_stats, short_name)
        @warn "No normalization stats for $short_name, returning unchanged"
        return var
    end
    
    mean_val, std_val = norm_stats[short_name]
    normalized_data = (var.data .- mean_val) ./ std_val
    
    return ClimaAnalysis.OutputVar(
        var.attributes,
        var.dims,
        var.dim_attributes,
        normalized_data,
    )
end

"""
    normalize_vars(vars_by_date, norm_stats)

Apply normalization to all variables in vars_by_date.
Returns a new Dict with normalized OutputVars.
"""
function normalize_vars(vars_by_date, norm_stats)
    normalized = Dict{Tuple{String, NTuple{2, Dates.DateTime}}, OutputVar}()
    
    for (key, var) in vars_by_date
        normalized[key] = normalize_var(var, norm_stats)
    end
    
    return normalized
end

"""
    resampled_lonlat(config_file)

Return a function to resample longitude and latitudes according to the model
grid specified by `config_file`.

For spectral element grids, the default interpolation grid is computed from h_elem:
  nlon = h_elem * 4 * 3 (cubed-sphere panels × spectral element degree)
  nlat = nlon / 2
"""
function resampled_lonlat(config_file)
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)
    if !isnothing(get(config_dict, "netcdf_interpolation_num_points", nothing))
        (nlon, nlat, nlev) = tuple(config_dict["netcdf_interpolation_num_points"]...)
    else
        # Compute from h_elem (spectral element grid)
        h_elem = get(config_dict, "h_elem", 12)
        # Default formula: h_elem * 4 panels * 3 (spectral degree)
        nlon = h_elem * 4 * 3
        nlat = nlon ÷ 2
        @info "Using model grid from h_elem=$h_elem: $(nlon)×$(nlat)"
    end
    lon_vals = range(-180, 180, nlon)
    lat_vals = range(-90, 90, nlat)
    return var ->
        ClimaAnalysis.resampled_as(var; longitude = lon_vals, latitude = lat_vals)
end

if abspath(PROGRAM_FILE) == @__FILE__
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"

    # Load calibration setup for NORMALIZE_VARIABLES and CERES_VARIABLES settings
    include(joinpath(@__DIR__, "calibration_setup.jl"))

    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    config_file = CALIBRATE_CONFIG.config_file
    obs_dir = ERA5_OBS_DIR
    output_path = joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal_weekly")

    # Determine which variables come from CERES vs ERA5
    ceres_vars_to_load = filter(v -> v in short_names, CERES_VARIABLES)
    era5_vars_to_load = filter(v -> !(v in CERES_VARIABLES), short_names)

    @info "Generating observations for $short_names"
    @info "ERA5 variables: $era5_vars_to_load (from: $obs_dir)"
    @info "CERES variables: $ceres_vars_to_load (from artifact)"
    @info "Normalization enabled: $NORMALIZE_VARIABLES"
    @info "The number of samples is $(length(sample_date_ranges)) over $sample_date_ranges"

    # Load observation data (ERA5 for some vars, CERES for radiation vars)
    unprocessed_vars = load_vars(obs_dir, short_names, config_file;
                                  ceres_variables = CERES_VARIABLES)
    @info "Loaded $(length(unprocessed_vars)) variable(s)"

    # Preprocess (set units)
    preprocessed_vars = preprocess_vars(unprocessed_vars)

    # Compute and apply normalization if enabled
    if NORMALIZE_VARIABLES
        norm_stats = compute_normalization_stats(preprocessed_vars, short_names)

        JLD2.save_object(joinpath(output_path, "norm_stats.jld2"), norm_stats)
        @info "Saved normalization stats to norm_stats.jld2"

        # Apply normalization to observations
        preprocessed_vars = normalize_vars(preprocessed_vars, norm_stats)
        @info "Applied normalization to observations"
    end

    lon_left = -60
    lon_right = 60
    preprocessed_vars = apply_lat_window(preprocessed_vars, lon_left, lon_right)

    # Save preprocessed variables (for inspection/debugging)
    JLD2.save_object(joinpath(output_path, "preprocessed_vars.jld2"), preprocessed_vars)

    # Create observation vector using ObservationRecipe (like subseasonal pipeline)
    observation_vector =
        make_observation_vector(preprocessed_vars, sample_date_ranges)
    JLD2.save_object(joinpath(output_path, "obs_vec.jld2"), observation_vector)

    @info "Saved observation vector with $(length(observation_vector)) samples"
end
