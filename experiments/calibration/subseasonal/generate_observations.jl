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
    # Reverse lookup: short_name -> ERA5 file prefix
    era5_prefix = nothing
    for (prefix, sn) in ERA5_FILE_PREFIX_TO_SHORT_NAME
        if sn == short_name
            era5_prefix = prefix
            break
        end
    end
    isnothing(era5_prefix) && error("Unknown short name: $short_name")
    
    start_str = Dates.format(start_date, "yyyymmdd")
    end_str = Dates.format(end_date, "yyyymmdd")
    filename = "$(era5_prefix)_weekly_mean_$(start_str)_$(end_str).nc"
    return joinpath(obs_dir, filename)
end

"""
    load_weekly_var(filepath, short_name, start_date, end_date)

Load a weekly ERA5 variable from a NetCDF file and set appropriate attributes.
The weekly files have no time dimension - just (latitude, longitude).
We add the start_date attribute that ClimaCalibrate expects.
"""
function load_weekly_var(filepath, short_name, start_date, end_date)
    # Get the ERA5 variable name from the short name
    era5_varname = nothing
    for (varname, sn) in ERA5_VARNAME_TO_SHORT_NAME
        if sn == short_name
            era5_varname = varname
            break
        end
    end
    isnothing(era5_varname) && error("Unknown short name: $short_name")
    
    @info "Loading $filepath with variable $era5_varname"
    var = OutputVar(filepath, era5_varname)
    
    var.attributes["short_name"] = short_name
    
    # Add start_date attribute that ClimaCalibrate needs
    # The weekly files don't have a time dimension, so we set it from the known date range
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
        for short_name in short_names
            filepath = get_weekly_filename(obs_dir, short_name, start_date, end_date)
            if !isfile(filepath)
                @warn "File not found: $filepath"
                continue
            end
            var = load_weekly_var(filepath, short_name, start_date, end_date)
            vars_by_date[(short_name, (start_date, end_date))] = var
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

For 2D static data (no time dimension), we create observations manually
by flattening the spatial data and using a scalar covariance.
"""
function make_observation_vector(vars_by_date, sample_date_ranges, short_names)
    # Covariance scalar - adjust based on expected observation noise
    noise_scalar = 5.0
    
    obs_vec = map(sample_date_ranges) do sample_date_range
        # Collect and flatten all variables for this date range
        all_data = Float64[]
        all_names = String[]
        
        for short_name in short_names
            key = (short_name, sample_date_range)
            if haskey(vars_by_date, key)
                var = vars_by_date[key]
                flat_data, _ = flatten_var_with_latitude_weights(var)
                # Filter out NaN/missing values
                valid_data = collect(skipmissing(flat_data))
                valid_data = filter(!isnan, valid_data)
                append!(all_data, valid_data)
                push!(all_names, short_name)
            else
                @warn "Missing variable for $key"
            end
        end
        
        # Create observation with scalar covariance (diagonal matrix)
        obs_name = join(all_names, "_") * "_" * Dates.format(first(sample_date_range), "yyyymmdd")
        noise = noise_scalar * LinearAlgebra.I
        
        @info "Creating observation '$obs_name' with $(length(all_data)) data points"
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
    lon = range(-180, 180, nlon)
    lat = range(-90, 90, nlat)
    return var -> ClimaAnalysis.resampled_as(var; lon, lat)
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
