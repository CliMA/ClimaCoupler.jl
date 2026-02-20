import ClimaAnalysis
import ClimaAnalysis: OutputVar
import ClimaCalibrate
import ClimaCoupler
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
    load_daily_var(filepath, short_name, target_date, file_date_range)

Load a specific day from a daily-mean ERA5 file.
The daily files have valid_time dimension with multiple days of data.
We extract just the day matching target_date.
"""
function load_daily_var(filepath, short_name, target_date, file_date_range)
    era5_varname = short_name_to_era5_varname(short_name)
    file_start, file_end = file_date_range
    @info "Loading daily $short_name for $target_date from $filepath"

    NCDatasets.Dataset(filepath) do ds
        data_3d = Array(ds[era5_varname])
        lats = Array(ds["latitude"])
        lons = Array(ds["longitude"])
        n_times = size(data_3d, 3)

        day_idx = Dates.value(Dates.Day(target_date - file_start)) + 1
        if day_idx < 1 || day_idx > n_times
            error(
                "target_date $target_date (index $day_idx) is outside file range $file_start to $file_end ($n_times days)",
            )
        end

        data_2d = data_3d[:, :, day_idx]

        # Create 3D var with time dimension for ObservationRecipe compatibility
        time_val = Float64(Dates.datetime2unix(target_date))

        dims = OrderedDict{String, Vector{Union{Missing, Float64}}}(
            "longitude" => Float64.(lons),
            "latitude" => Float64.(lats),
            "time" => [time_val],
        )

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
            "time" => Dict{String, Any}(
                "units" => "seconds since 1970-01-01",
                "calendar" => "standard",
            ),
        )

        attribs = Dict{String, Any}(
            "short_name" => short_name,
            "start_date" => target_date,
            "units" => get(ds[era5_varname].attrib, "units", ""),
            "long_name" => get(ds[era5_varname].attrib, "long_name", short_name),
        )

        data_3d_reshaped = reshape(Float64.(data_2d), size(data_2d)..., 1)
        var = ClimaAnalysis.OutputVar(attribs, dims, dim_attribs, data_3d_reshaped)

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

    # Ensure time dimension exists for ObservationRecipe compatibility
    if !haskey(var.dims, "time")
        time_val = Float64(Dates.datetime2unix(start_date))
        var.dims["time"] = [time_val]
        var.dim_attributes["time"] = Dict{String, Any}(
            "units" => "seconds since 1970-01-01",
            "calendar" => "standard",
        )
        var = ClimaAnalysis.OutputVar(
            var.attributes,
            var.dims,
            var.dim_attributes,
            reshape(var.data, size(var.data)..., 1),
        )
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
        is_single_day = (start_date == end_date)

        for short_name in short_names
            try
                if is_single_day
                    filepath, file_date_range =
                        get_daily_filename(obs_dir, short_name, start_date)
                    var =
                        load_daily_var(filepath, short_name, start_date, file_date_range)
                else
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
            end
        end
    end

    return vars_by_date
end

"""
    preprocess_vars(vars_by_date, config_file)

Preprocess each OutputVar by resampling to the model grid and setting units.

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
"""
function make_observation_vector(vars_by_date, sample_date_ranges, short_names)
    include(joinpath(@__DIR__, "calibration_priors.jl"))

    covar_estimator = ClimaCalibrate.ObservationRecipe.ScalarCovariance(;
        scalar = Float64(CALIBRATION_NOISE_SCALAR),
        use_latitude_weights = true,
    )

    obs_vec = map(sample_date_ranges) do sample_date_range
        # Collect OutputVars for this date range
        date_vars = OutputVar[]
        for sn in short_names
            key = (sn, sample_date_range)
            if haskey(vars_by_date, key)
                push!(date_vars, vars_by_date[key])
            else
                @warn "Missing variable for $key"
            end
        end

        ClimaCalibrate.ObservationRecipe.observation(
            covar_estimator,
            date_vars,
            first(sample_date_range),
            last(sample_date_range),
        )
    end
    return obs_vec
end

"""
    resampled_lonlat(config_file)

Return a function to resample longitude and latitudes according to the model
grid specified by `config_file`.
(Same as subseasonal pipeline)
"""
function resampled_lonlat(config_file)
    config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)
    if !isnothing(get(config_dict, "netcdf_interpolation_num_points", nothing))
        (nlon, nlat, nlev) = tuple(config_dict["netcdf_interpolation_num_points"]...)
    else
        nlon, nlat = 360, 180
    end
    lon_vals = range(-180, 180, nlon)
    lat_vals = range(-90, 90, nlat)
    return var ->
        ClimaAnalysis.resampled_as(var; longitude = lon_vals, latitude = lat_vals)
end

if abspath(PROGRAM_FILE) == @__FILE__
    ENV["CLIMACOMMS_CONTEXT"] = "SINGLETON"

    sample_date_ranges = CALIBRATE_CONFIG.sample_date_ranges
    short_names = CALIBRATE_CONFIG.short_names
    config_file = CALIBRATE_CONFIG.config_file
    obs_dir = ERA5_OBS_DIR

    @info "Generating observations for $short_names"
    @info "Using ERA5 data from: $obs_dir"
    @info "The number of samples is $(length(sample_date_ranges)) over $sample_date_ranges"

    # Load weekly ERA5 files
    unprocessed_vars = load_vars(obs_dir, short_names, sample_date_ranges)
    @info "Loaded $(length(unprocessed_vars)) variable(s)"

    # Preprocess (resample to model grid and set units)
    preprocessed_vars = preprocess_vars(unprocessed_vars, config_file)

    # Save preprocessed variables (for inspection/debugging)
    JLD2.save_object(
        joinpath(
            pkgdir(ClimaCoupler),
            "experiments/calibration/subseasonal_weekly/preprocessed_vars.jld2",
        ),
        preprocessed_vars,
    )

    # Create observation vector using ObservationRecipe (like subseasonal pipeline)
    observation_vector =
        make_observation_vector(preprocessed_vars, sample_date_ranges, short_names)
    JLD2.save_object(
        joinpath(pkgdir(ClimaCoupler), "experiments/calibration/subseasonal_weekly/obs_vec.jld2"),
        observation_vector,
    )

    @info "Saved observation vector with $(length(observation_vector)) samples"
end
