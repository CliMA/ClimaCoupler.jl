import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCoupler
import Dates
import EnsembleKalmanProcesses as EKP
import JLD2
import ClimaDiagnostics
import ClimaCore
import NCDatasets
import LinearAlgebra
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
)

# Mapping from ERA5 NetCDF variable names to short names
const ERA5_VARNAME_TO_SHORT_NAME = Dict(
    "t2m" => "tas",
    "msl" => "mslp",
    "tp" => "pr",
    "sst" => "ts",
    "sp" => "sp",
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
    # Daily files are in parent directory of weekly
    parent_dir = dirname(obs_dir)
    
    # Pattern to match: prefix_daily-mean_YYYYMMDD_YYYYMMDD.nc
    pattern = Regex("$(era5_prefix)_daily-mean_(\\d{8})_(\\d{8})\\.nc")
    
    # Scan directory for matching files
    for filename in readdir(parent_dir)
        m = match(pattern, filename)
        if !isnothing(m)
            file_start = Dates.DateTime(m.captures[1], "yyyymmdd")
            file_end = Dates.DateTime(m.captures[2], "yyyymmdd")
            # Check if target_date falls within this file's range
            if file_start <= target_date <= file_end
                return joinpath(parent_dir, filename), (file_start, file_end)
            end
        end
    end
    
    error("No daily-mean file found containing date $target_date for variable $short_name in $parent_dir")
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
    load_daily_var(filepath, short_name, target_date, file_date_range)

Load a specific day from a daily-mean ERA5 file.
The daily files have valid_time dimension with multiple days of data.
We extract just the day matching target_date.

Arguments:
- filepath: Path to the NetCDF file
- short_name: Variable short name (e.g., "tas")
- target_date: The specific date to extract
- file_date_range: (start_date, end_date) tuple indicating the file's date coverage
"""
function load_daily_var(filepath, short_name, target_date, file_date_range)
    era5_varname = short_name_to_era5_varname(short_name)
    file_start, file_end = file_date_range
    @info "Loading daily file $filepath, variable $era5_varname, extracting date $target_date"
    
    # Open the NetCDF file and extract the specific day
    NCDatasets.Dataset(filepath) do ds
        # Get the variable data - Julia/NCDatasets uses column-major order
        # NetCDF dimensions (valid_time, latitude, longitude) become Julia array (lon, lat, time)
        data_3d = Array(ds[era5_varname])
        lats = Array(ds["latitude"])
        lons = Array(ds["longitude"])
        n_times = size(data_3d, 3)  # Time is the 3rd dimension in Julia
        
        # Calculate day index based on file's actual date range
        day_idx = Dates.value(Dates.Day(target_date - file_start)) + 1
        
        if day_idx < 1 || day_idx > n_times
            error("target_date $target_date (index $day_idx) is outside file range $file_start to $file_end ($n_times days)")
        end
        
        @info "Extracting day index $day_idx for $target_date (file covers $file_start to $file_end)"
        
        # Extract the 2D slice for this day - data is (lon, lat, time)
        data_2d = data_3d[:, :, day_idx]
        
        # Create dims dict for OutputVar - must use OrderedDict with correct dimension names
        # and order matching the data layout (lon, lat)
        dims = OrderedDict{String, Vector{Union{Missing, Float64}}}(
            "longitude" => Float64.(lons),
            "latitude" => Float64.(lats),
        )
        
        # Create dim_attributes with required entries for ClimaAnalysis
        dim_attribs = OrderedDict{String, Dict{String, Any}}(
            "longitude" => Dict{String, Any}(
                "units" => "degrees_east",
                "long_name" => "longitude",
                "standard_name" => "longitude",
            ),
            "latitude" => Dict{String, Any}(
                "units" => "degrees_north",
                "long_name" => "latitude",
                "standard_name" => "latitude",
            ),
        )
        
        # Create attributes - must be Dict{String, Any} for ClimaAnalysis
        attribs = Dict{String, Any}(
            "short_name" => short_name,
            "start_date" => target_date,
            "units" => get(ds[era5_varname].attrib, "units", ""),
            "long_name" => get(ds[era5_varname].attrib, "long_name", short_name),
        )
        
        # Create OutputVar - data is already (lon, lat) from Julia's column-major ordering
        var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, Float64.(data_2d))
        
        # Ensure latitudes are sorted in ascending order
        if !issorted(ClimaAnalysis.latitudes(var))
            var = ClimaAnalysis.reverse_dim(var, ClimaAnalysis.latitude_name(var))
        end
        
        return var
    end
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
The weekly files have no time dimension - just (latitude, longitude).
We add the start_date attribute that ClimaCalibrate expects.
"""
function load_weekly_var(filepath, short_name, start_date, end_date)
    era5_varname = short_name_to_era5_varname(short_name)
    @info "Loading weekly file $filepath with variable $era5_varname"
    var = OutputVar(filepath, era5_varname)
    
    var.attributes["short_name"] = short_name
    var.attributes["start_date"] = start_date
    
    # Ensure latitudes are sorted in ascending order
    if !issorted(ClimaAnalysis.latitudes(var))
        var = ClimaAnalysis.reverse_dim(var, ClimaAnalysis.latitude_name(var))
    end
    
    return var
end

"""
    load_vars(obs_dir, short_names, sample_date_ranges)

Load weekly ERA5 NetCDF files for the specified short_names and date ranges.
Returns a Dict mapping (short_name, date_range) -> OutputVar
"""
function load_vars(obs_dir, short_names, sample_date_ranges)
    vars_by_date = Dict{Tuple{String, NTuple{2, Dates.DateTime}}, OutputVar}()
    
    for (start_date, end_date) in sample_date_ranges
        # Determine if this is a single-day or multi-day range
        is_single_day = (start_date == end_date)
        
        for short_name in short_names
            try
                if is_single_day
                    # Use daily-mean files for single-day ranges
                    filepath, file_date_range = get_daily_filename(obs_dir, short_name, start_date)
                    var = load_daily_var(filepath, short_name, start_date, file_date_range)
                else
                    # Use weekly files for multi-day ranges
                    filepath = get_weekly_filename(obs_dir, short_name, start_date, end_date)
                    if !isfile(filepath)
                        @warn "Weekly file not found: $filepath"
                        continue
                    end
                    var = load_weekly_var(filepath, short_name, start_date, end_date)
                end
                vars_by_date[(short_name, (start_date, end_date))] = var
            catch e
                @warn "Failed to load $short_name for $start_date - $end_date: $e"
            end
        end
    end
    
    return vars_by_date
end

"""
    preprocess_vars(vars_by_date, config_file)

Preprocess each OutputVar by resampling to the model grid and setting units.
"""
function preprocess_vars(vars_by_date, config_file)
    resample_var = resampled_lonlat(config_file)
    processed = Dict{Tuple{String, NTuple{2, Dates.DateTime}}, OutputVar}()
    
    for (key, var) in vars_by_date
        short_name = key[1]
        date_range = key[2]
        start_date = date_range[1]
        
        var = resample_var(var)
        if haskey(var_units, short_name)
            var = ClimaAnalysis.set_units(var, var_units[short_name])
        end
        # Ensure start_date attribute is preserved after resampling
        var.attributes["short_name"] = short_name
        var.attributes["start_date"] = start_date
        
        processed[key] = var
    end
    
    return processed
end

"""
    flatten_var_with_latitude_weights(var)

Flatten a 2D OutputVar to a vector, applying latitude (cosine) weights.
Returns (flattened_data, weights).
"""
function flatten_var_with_latitude_weights(var)
    lats = ClimaAnalysis.latitudes(var)
    lons = ClimaAnalysis.longitudes(var)
    
    # Cosine latitude weights
    weights = cos.(deg2rad.(lats))
    weights ./= sum(weights)  # Normalize
    
    # Get data array - should be (lon, lat) or (lat, lon)
    data = var.data
    
    # Flatten data
    flat_data = vec(data)
    
    # Create weight matrix matching the data shape
    # Assuming data is (lon, lat)
    nlon = length(lons)
    nlat = length(lats)
    weight_matrix = repeat(weights', nlon, 1)
    flat_weights = vec(weight_matrix)
    
    return flat_data, flat_weights
end

"""
    make_observation_vector(vars_by_date, sample_date_ranges, short_names)

Create the observation vector for EKP from the preprocessed variables.
Each sample in the vector corresponds to one date range with all short_names.

Data is normalized to zero mean and unit variance for stable EKP convergence.
Normalization constants are saved for use in observation_map.jl.
"""
function make_observation_vector(vars_by_date, sample_date_ranges, short_names)
    # First pass: collect all data to compute mean/std for normalization
    all_data_raw = Float64[]
    for sample_date_range in sample_date_ranges
        for short_name in short_names
            key = (short_name, sample_date_range)
            if haskey(vars_by_date, key)
                var = vars_by_date[key]
                flat_data, _ = flatten_var_with_latitude_weights(var)
                valid_data = filter(!isnan, collect(skipmissing(flat_data)))
                append!(all_data_raw, valid_data)
            end
        end
    end
    
    # Compute normalization constants
    data_mean = Statistics.mean(all_data_raw)
    data_std = Statistics.std(all_data_raw)
    
    @info "Normalizing observations: mean=$data_mean, std=$data_std"
    
    # Save normalization constants for use in observation_map
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/normalization.jld2"),
        (mean=data_mean, std=data_std)
    )
    
    # Noise scalar for normalized data (unit variance, so 1.0 = 100% of std)
    # Use smaller value for tighter fit, larger for slower convergence
    noise_scalar = 0.5  # 50% of unit variance - adjust to control convergence speed
    
    obs_vec = map(sample_date_ranges) do sample_date_range
        all_data = Float64[]
        all_names = String[]
        
        for short_name in short_names
            key = (short_name, sample_date_range)
            if haskey(vars_by_date, key)
                var = vars_by_date[key]
                flat_data, _ = flatten_var_with_latitude_weights(var)
                valid_data = filter(!isnan, collect(skipmissing(flat_data)))
                # Normalize the data
                normalized_data = (valid_data .- data_mean) ./ data_std
                append!(all_data, normalized_data)
                push!(all_names, short_name)
            else
                @warn "Missing variable for $key"
            end
        end
        
        obs_name = join(all_names, "_") * "_" * Dates.format(first(sample_date_range), "yyyymmdd")
        noise = noise_scalar * LinearAlgebra.I
        
        @info "Creating observation '$obs_name' with $(length(all_data)) normalized data points"
        EKP.Observation(all_data, noise, obs_name)
    end
    return obs_vec
end

# OLD VERSION (unnormalized) - kept for reference
# """
#     make_observation_vector_unnormalized(vars_by_date, sample_date_ranges, short_names)

# # """
# function make_observation_vector_unnormalized(vars_by_date, sample_date_ranges, short_names)
#     # Covariance scalar - adjust based on expected observation noise
#     noise_scalar = 5.0
    
#     obs_vec = map(sample_date_ranges) do sample_date_range
#         # Collect and flatten all variables for this date range
#         all_data = Float64[]
#         all_names = String[]
        
#         for short_name in short_names
#             key = (short_name, sample_date_range)
#             if haskey(vars_by_date, key)
#                 var = vars_by_date[key]
#                 flat_data, _ = flatten_var_with_latitude_weights(var)
#                 # Filter out NaN/missing values
#                 valid_data = collect(skipmissing(flat_data))
#                 valid_data = filter(!isnan, valid_data)
#                 append!(all_data, valid_data)
#                 push!(all_names, short_name)
#             else
#                 @warn "Missing variable for $key"
#             end
#         end
        
#         # Create observation with scalar covariance (diagonal matrix)
#         obs_name = join(all_names, "_") * "_" * Dates.format(first(sample_date_range), "yyyymmdd")
#         noise = noise_scalar * LinearAlgebra.I
        
#         @info "Creating observation '$obs_name' with $(length(all_data)) data points"
#         EKP.Observation(all_data, noise, obs_name)
#     end
#     return obs_vec
# end

"""
    resampled_lonlat(config_file)

Return a function to resample longitude and latitudes according to the model
grid specified by `config_file`.
"""
function resampled_lonlat(config_file)
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)
    if !isnothing(get(config_dict, "netcdf_interpolation_num_points", nothing))
        (nlon, nlat, nlev) = tuple(config_dict["netcdf_interpolation_num_points"]...)
    else
        # Default grid resolution if not specified
        nlon, nlat = 360, 180
    end
    lon_vals = range(-180, 180, nlon)
    lat_vals = range(-90, 90, nlat)
    # Use full dimension names (longitude/latitude) as ERA5 files use these
    return var -> ClimaAnalysis.resampled_as(var; longitude = lon_vals, latitude = lat_vals)
end

if abspath(PROGRAM_FILE) == @__FILE__
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"
    
    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    config_file = CALIBRATE_CONFIG.config_file
    obs_dir = CALIBRATE_CONFIG.obs_dir
    
    @info "Generating observations for $short_names"
    @info "Using ERA5 data from: $obs_dir"
    @info "The number of samples is $(length(sample_date_ranges)) over $sample_date_ranges"

    # Load weekly ERA5 files
    unprocessed_vars = load_vars(obs_dir, short_names, sample_date_ranges)
    @info "Loaded $(length(unprocessed_vars)) variable(s)"

    # Preprocess (resample and set units)
    preprocessed_vars = preprocess_vars(unprocessed_vars, config_file)

    # Save preprocessed variables
    JLD2.save_object(
        joinpath(
            pkgdir(ClimaCoupler),
            "experiments/calibration/subseasonal/preprocessed_vars.jld2",
        ),
        preprocessed_vars,
    )
    
    # Create observation vector for EKP
    observation_vector = make_observation_vector(preprocessed_vars, sample_date_ranges, short_names)
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/obs_vec.jld2"),
        observation_vector,
    )
    
    @info "Saved observation vector with $(length(observation_vector)) samples"
end
