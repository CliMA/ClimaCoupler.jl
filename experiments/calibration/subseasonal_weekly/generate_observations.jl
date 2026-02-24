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
    var = ClimaAnalysis.select(var, by = ClimaAnalysis.MatchValue(), time = target_date)

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

"""
    add_time_dimension(var, start_date)

Add a time dimension to a 2D OutputVar for ObservationRecipe compatibility.
Time is stored as seconds since start_date, so ClimaAnalysis.dates() can 
convert it back to DateTime using the start_date attribute.
"""
function add_time_dimension(var, start_date)
    if haskey(var.dims, "time")
        return var  # Already has time dimension
    end
    
    # Time value = 0.0 means "at start_date"
    # ClimaAnalysis.dates() will compute: start_date + time_value seconds
    time_val = 0.0
    
    # Create new dims with time, preserving order (lon, lat, time)
    # and converting to Float64 vectors (ClimaCalibrate expects plain Float64)
    new_dims = OrderedDict{String, Vector{Float64}}()
    for k in keys(var.dims)
        new_dims[k] = Float64.(collect(skipmissing(var.dims[k])))
    end
    new_dims["time"] = [time_val]
    
    # Create new dim_attributes with time
    # ClimaCalibrate expects "s" as the unit string for seconds
    new_dim_attribs = copy(var.dim_attributes)
    new_dim_attribs["time"] = Dict{String, Any}(
        "units" => "s",
        "calendar" => "standard",
    )
    
    # Reshape data to add time dimension, ensure Float64
    new_data = reshape(Float64.(var.data), size(var.data)..., 1)
    
    return ClimaAnalysis.OutputVar(
        var.attributes,
        new_dims,
        new_dim_attribs,
        new_data,
    )
end


##### CERES Data Loading Functions

"""
    load_ceres_var(short_name, start_date)

Load a variable from CERES monthly data for the month containing start_date.
Returns a 2D OutputVar (lon, lat) representing the monthly mean.

Note: CERES data is monthly, so we use the monthly value for the month
containing the calibration period. Time dimension is added later in preprocess_vars.
"""
function load_ceres_var(short_name, start_date)
    loader = CalibrationTools.CERESDataLoader()
    
    # Load the full time series
    var = Base.get(loader, short_name)
    # Select the month containing start_date (CERES dates are at start of month)
    month_start = Dates.firstdayofmonth(start_date)
    var = ClimaAnalysis.select(var, by = ClimaAnalysis.MatchValue(), time = month_start)
    @info "Loaded CERES $short_name for $(Dates.monthname(start_date)) $(Dates.year(start_date))"
    
    return var
end

"""
    is_ceres_variable(short_name, ceres_variables)

Check if a variable should be loaded from CERES instead of ERA5.
"""
is_ceres_variable(short_name, ceres_variables) = short_name in ceres_variables

"""
    load_vars(obs_dir, short_names, sample_date_ranges; ceres_variables = String[])

Load observation data for the specified short_names and date ranges.
- Variables in `ceres_variables` are loaded from CERES monthly data
- Other variables are loaded from ERA5 (daily or weekly files)

Returns a Dict mapping (short_name, date_range) -> OutputVar
"""
function load_vars(obs_dir, short_names, sample_date_ranges; ceres_variables = String[])
    vars_by_date = Dict{Tuple{String, NTuple{2, Dates.DateTime}}, OutputVar}()
    
    for (start_date, end_date) in sample_date_ranges
        is_single_day = (start_date == end_date)
        
        for short_name in short_names
            try
                # Check if this variable should come from CERES
                if is_ceres_variable(short_name, ceres_variables)
                    var = load_ceres_var(short_name, start_date)
                elseif is_single_day
                    # Load from ERA5 daily files
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
                vars_by_date[(short_name, (start_date, end_date))] = var
            catch e
                @warn "Failed to load $short_name for $start_date - $end_date: $e"
                rethrow(e)
            end
        end
    end
    
    return vars_by_date
end

"""
    preprocess_vars(vars_by_date, config_file)

Preprocess each OutputVar by resampling to the model grid, setting units,
and adding time dimension for ObservationRecipe compatibility.

Precipitation (pr) requires special handling:
ERA5 'tp' is in meters/hour, we convert to mm/day using ERA5_PR_CONVERSION.
"""
function preprocess_vars(vars_by_date, config_file)
    resample_var = resampled_lonlat(config_file)
    processed = Dict{Tuple{String, NTuple{2, Dates.DateTime}}, OutputVar}()
    
    for (key, var) in vars_by_date
        short_name = key[1]
        date_range = key[2]
        start_date = date_range[1]
        
        # Resample to model grid (must be done before adding time dimension)
        var = resample_var(var)
        
        # Apply precipitation unit conversion for ERA5 data
        if short_name == "pr"
            new_data = var.data .* ERA5_PR_CONVERSION
            @info "ERA5 pr conversion: factor=$ERA5_PR_CONVERSION → mean=$(round(Statistics.mean(filter(!isnan, vec(new_data))), sigdigits=4)) mm/day"
            var = ClimaAnalysis.OutputVar(
                var.attributes,
                var.dims,
                var.dim_attributes,
                new_data,
            )
        end
        
        if haskey(var_units, short_name)
            var = ClimaAnalysis.set_units(var, var_units[short_name])
        end
        var.attributes["short_name"] = short_name
        var.attributes["start_date"] = start_date
        
        # Add time dimension AFTER resampling for ObservationRecipe compatibility
        var = add_time_dimension(var, start_date)
        
        processed[key] = var
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
function make_observation_vector(vars_by_date, sample_date_ranges, short_names)
    include(joinpath(@__DIR__, "calibration_setup.jl"))

    covar_estimator = ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
        scalar = Float64(CALIBRATION_NOISE_SCALAR),
        use_latitude_weights = true,
    )

    obs_vec = map(sample_date_ranges) do sample_date_range
        # Collect OutputVars for this date range
        date_vars = OutputVar[]
        start_date = first(sample_date_range)
        for sn in short_names
            key = (sn, sample_date_range)
            if haskey(vars_by_date, key)
                push!(date_vars, vars_by_date[key])
            else
                @warn "Missing variable for $key"
            end
        end
        
        # Use start_date for both args since weekly data has single time point
        # The data represents the weekly average, stored at start_date
        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            date_vars,
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
    unprocessed_vars = load_vars(obs_dir, short_names, sample_date_ranges; 
                                  ceres_variables = CERES_VARIABLES)
    @info "Loaded $(length(unprocessed_vars)) variable(s)"

    # Preprocess (resample to model grid and set units)
    preprocessed_vars = preprocess_vars(unprocessed_vars, config_file)

    # Compute and apply normalization if enabled
    if NORMALIZE_VARIABLES
        norm_stats = compute_normalization_stats(preprocessed_vars, short_names)

        JLD2.save_object(joinpath(output_path, "norm_stats.jld2"), norm_stats)
        @info "Saved normalization stats to norm_stats.jld2"
        
        # Apply normalization to observations
        preprocessed_vars = normalize_vars(preprocessed_vars, norm_stats)
        @info "Applied normalization to observations"
    end

    # Save preprocessed variables (for inspection/debugging)
    JLD2.save_object(joinpath(output_path, "preprocessed_vars.jld2"), preprocessed_vars)
    
    # Create observation vector using ObservationRecipe (like subseasonal pipeline)
    observation_vector =
        make_observation_vector(preprocessed_vars, sample_date_ranges, short_names)
    JLD2.save_object(joinpath(output_path, "obs_vec.jld2"), observation_vector)
    
    @info "Saved observation vector with $(length(observation_vector)) samples"
end
