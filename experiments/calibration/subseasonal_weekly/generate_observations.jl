import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCoupler
import ClimaComms
import Dates
import EnsembleKalmanProcesses as EKP
import JLD2
import ClimaDiagnostics
import ClimaCore
import NCDatasets
import LinearAlgebra
import Statistics
import Interpolations
using OrderedCollections: OrderedDict

# Access CalibrateConfig (this also includes observation_utils.jl via observation_map.jl)
include(joinpath(@__DIR__, "run_calibration.jl"))

"""
    load_land_mask(target_lons, target_lats; threshold=0.5)

Load ERA5 land fraction and resample to target grid, returning a binary land mask.
Points with land fraction > threshold are considered land (true).
"""
function load_land_mask(target_lons, target_lats; threshold=0.5)
    # Get ERA5 land fraction artifact
    artifact_path = ClimaCoupler.Input.@clima_artifact("era5_land_fraction", ClimaComms.context())
    land_frac_file = joinpath(artifact_path, "era5_land_fraction.nc")
    
    @info "Loading land mask from $land_frac_file"
    
    land_mask = NCDatasets.Dataset(land_frac_file) do ds
        # ERA5 land fraction: lsm(lat, lon)
        lsm = ds["lsm"][:, :]  # (lat, lon) in file, but Julia reads as (lon, lat)
        src_lons = ds["lon"][:]
        src_lats = ds["lat"][:]
        
        # NCDatasets reads as (lon, lat) due to column-major ordering
        # lsm is actually (lon, lat) in Julia memory
        
        # Create interpolation (need to handle lon wrapping)
        # ERA5 lons are typically 0-360 or -180 to 180
        # Make sure we handle the interpolation correctly
        
        # Simple bilinear interpolation
        itp = Interpolations.interpolate(
            (src_lons, src_lats), 
            lsm, 
            Interpolations.Gridded(Interpolations.Linear())
        )
        itp_extrap = Interpolations.extrapolate(itp, Interpolations.Flat())
        
        # Resample to target grid
        resampled = [itp_extrap(lon, lat) for lon in target_lons, lat in target_lats]
        
        return resampled
    end
    
    # Convert to binary mask (true = land)
    binary_mask = land_mask .> threshold
    
    n_land = count(binary_mask)
    n_total = length(binary_mask)
    @info "Land mask: $n_land / $n_total points are land ($(round(100*n_land/n_total, digits=1))%)"
    
    return binary_mask
end

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

# Variables that should be evaluated only over land (others are evaluated everywhere)
const LAND_ONLY_VARS = Set(["tas", "pr"])


# Mapping from ERA5 NetCDF variable names to short names
# Note: rsut/rlut variable names (mtuswrf/mtulwrf) may need adjustment based on actual file contents
const ERA5_VARNAME_TO_SHORT_NAME = Dict(
    "t2m" => "tas",
    "msl" => "mslp",
    "tp" => "pr",
    "sst" => "ts",
    "sp" => "sp",
    "mean_top_upward_short_wave_radiation_flux" => "rsut",  # Mean top upward short-wave radiation flux
    "mean_top_upward_long_wave_radiation_flux" => "rlut",  # Mean top upward long-wave radiation flux
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

Preprocess each OutputVar by resampling to the model grid, applying unit conversions,
and setting units metadata.

**IMPORTANT**: Precipitation (pr) requires special handling:
- ERA5 'tp' is in meters/hour, we convert to mm/day using ERA5_PR_CONVERSION
"""
function preprocess_vars(vars_by_date, config_file)
    resample_var = resampled_lonlat(config_file)
    processed = Dict{Tuple{String, NTuple{2, Dates.DateTime}}, OutputVar}()
    
    for (key, var) in vars_by_date
        short_name = key[1]
        date_range = key[2]
        start_date = date_range[1]
        
        var = resample_var(var)
        
        # Apply precipitation unit conversion for ERA5 data
        if short_name == "pr"
            # ERA5 'tp' is in meters/hour, convert to mm/day
            raw_mean = Statistics.mean(filter(!isnan, vec(var.data)))
            raw_min = minimum(filter(!isnan, vec(var.data)))
            raw_max = maximum(filter(!isnan, vec(var.data)))
            
            new_data = var.data .* ERA5_PR_CONVERSION
            
            conv_mean = Statistics.mean(filter(!isnan, vec(new_data)))
            conv_min = minimum(filter(!isnan, vec(new_data)))
            conv_max = maximum(filter(!isnan, vec(new_data)))
            
            @info "ERA5 pr conversion (m/hour → mm/day, factor=$ERA5_PR_CONVERSION):"
            @info "  Before: mean=$(round(raw_mean, sigdigits=4)), range=[$raw_min, $raw_max]"
            @info "  After:  mean=$(round(conv_mean, sigdigits=4)), range=[$conv_min, $conv_max] mm/day"
            
            var = ClimaAnalysis.OutputVar(var.attributes, var.dims, var.dim_attributes, new_data)
        end
        
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
    flatten_var(var; land_mask=nothing)

Flatten a 2D OutputVar to a vector.

If `land_mask` is provided (a 2D boolean array matching var dimensions),
only land points (where land_mask == true) are kept.
"""
function flatten_var(var; land_mask=nothing)
    data = var.data
    flat_data = vec(data)
    
    # Apply land mask if provided
    if !isnothing(land_mask)
        flat_mask = vec(land_mask)
        land_indices = findall(flat_mask)
        flat_data = flat_data[land_indices]
    end
    
    return flat_data
end

"""
    make_observation_vector(vars_by_date, sample_date_ranges, short_names)

Create the observation vector for EKP from the preprocessed variables.
Each sample in the vector corresponds to one date range with all short_names.

**PER-VARIABLE NORMALIZATION**: Each variable is normalized independently to zero mean
and unit variance. This ensures all variables (tas, mslp, pr) contribute equally to 
the calibration objective, regardless of their different physical scales.

Land masking is applied per-variable based on LAND_ONLY_VARS:
- Variables in LAND_ONLY_VARS (tas, pr): only land points included
- Other variables (mslp): all points included (evaluated everywhere)
"""
function make_observation_vector(vars_by_date, sample_date_ranges, short_names)
    # Load config from single source of truth
    include(joinpath(@__DIR__, "calibration_priors.jl"))
    noise_scalar = CALIBRATION_NOISE_SCALAR
    
    @info "Using noise_scalar=$noise_scalar from calibration_priors.jl"
    
    # Get grid from first variable to determine land mask resolution
    first_key = first(keys(vars_by_date))
    first_var = vars_by_date[first_key]
    target_lons = ClimaAnalysis.longitudes(first_var)
    target_lats = ClimaAnalysis.latitudes(first_var)
    
    # Load land mask for variables that need it
    land_mask = load_land_mask(target_lons, target_lats)
    
    # Save land mask and land_only_vars for use in observation_map
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/land_mask.jld2"),
        (mask=land_mask, lons=collect(target_lons), lats=collect(target_lats), 
         land_only_vars=collect(LAND_ONLY_VARS))
    )
    @info "Saved land mask to land_mask.jld2 (land-only vars: $(LAND_ONLY_VARS))"
    
    return make_observation_vector_gridpoints(vars_by_date, sample_date_ranges, short_names, land_mask, noise_scalar)
end

"""
    make_observation_vector_gridpoints(...)

Create observation vector using full gridpoint data (original behavior).
Returns 113k+ values per sample.
"""
function make_observation_vector_gridpoints(vars_by_date, sample_date_ranges, short_names, land_mask, noise_scalar)
    # ==========================================================================
    # PER-VARIABLE NORMALIZATION: Compute mean/std for EACH variable separately
    # This ensures each variable contributes equally to the objective function
    # ==========================================================================
    var_normalization = Dict{String, NamedTuple{(:mean, :std), Tuple{Float64, Float64}}}()
    
    for short_name in short_names
        # Collect all data for this variable across all date ranges
        var_data_raw = Float64[]
        for sample_date_range in sample_date_ranges
            key = (short_name, sample_date_range)
            if haskey(vars_by_date, key)
                var = vars_by_date[key]
                # Apply land mask only for land-only variables
                var_land_mask = short_name in LAND_ONLY_VARS ? land_mask : nothing
                flat_data = flatten_var(var; land_mask=var_land_mask)
                valid_data = filter(!isnan, collect(skipmissing(flat_data)))
                append!(var_data_raw, valid_data)
            end
        end
        
        # Compute per-variable normalization constants
        var_mean = Statistics.mean(var_data_raw)
        var_std = Statistics.std(var_data_raw)
        var_normalization[short_name] = (mean=var_mean, std=var_std)
        
        mask_status = short_name in LAND_ONLY_VARS ? "land only" : "all points"
        @info "  $short_name: mean=$(round(var_mean, sigdigits=6)), std=$(round(var_std, sigdigits=6)), $(length(var_data_raw)) points ($mask_status)"
    end
    
    # Save per-variable normalization constants for use in observation_map
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal/normalization.jld2"),
        var_normalization
    )
    @info "Saved per-variable normalization constants to normalization.jld2"
    
    obs_vec = map(sample_date_ranges) do sample_date_range
        all_data = Float64[]
        all_names = String[]
        
        for short_name in short_names
            key = (short_name, sample_date_range)
            if haskey(vars_by_date, key)
                var = vars_by_date[key]
                # Apply land mask only for land-only variables
                var_land_mask = short_name in LAND_ONLY_VARS ? land_mask : nothing
                flat_data = flatten_var(var; land_mask=var_land_mask)
                valid_data = filter(!isnan, collect(skipmissing(flat_data)))
                
                # Apply PER-VARIABLE normalization
                norm = var_normalization[short_name]
                normalized_data = (valid_data .- norm.mean) ./ norm.std
                append!(all_data, normalized_data)
                push!(all_names, short_name)
                
                mask_status = isnothing(var_land_mask) ? "all points" : "land only"
                @info "    $short_name: $(length(valid_data)) points ($mask_status), normalized with mean=$(round(norm.mean, sigdigits=4)), std=$(round(norm.std, sigdigits=4))"
            else
                @warn "Missing variable for $key"
            end
        end
        
        obs_name = join(all_names, "_") * "_" * Dates.format(first(sample_date_range), "yyyymmdd")
        noise = noise_scalar * LinearAlgebra.I
        
        @info "Creating observation '$obs_name' with $(length(all_data)) normalized data points (per-variable normalization)"
        EKP.Observation(all_data, noise, obs_name)
    end
    return obs_vec
end

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
