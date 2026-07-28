"""
    Era5InitialConditions

Fetch, preprocess, validate, and cache ERA5 initial conditions for
subseasonal (WeatherModel) simulations.

Given a simulation start date, this module downloads ERA5 reanalysis data
from the Copernicus Climate Data Store (CDS) with CDSAPI.jl and produces the
NetCDF files that the coupler, ClimaLand, and ClimaAtmos read:

- `era5_raw_YYYYMMDD_0000.nc`: pressure-level atmosphere state for ClimaAtmos
- `sst_processed_YYYYMMDD_0000.nc`: sea surface temperature (Celsius)
- `sic_processed_YYYYMMDD_0000.nc`: sea ice concentration (percent)
- `era5_land_processed_YYYYMMDD_0000.nc`: integrated land initial conditions
- `era5_bucket_processed_YYYYMMDD_0000.nc`: bucket land initial conditions
- `albedo_processed_YYYYMMDD_0000.nc`: surface albedo

Files are cached in a Scratch.jl directory (override with the
`CLIMACOUPLER_ERA5_CACHE_DIR` environment variable), so repeat runs for the
same date make no network requests.

The CDS requests and the preprocessing are ported from the private
CliMA/WeatherQuest repository, which
builds the `wxquest_initial_conditions` artifact. Ported functions carry a
provenance note in their docstring. The fetch, cache, and validation
machinery is new. Unlike WeatherQuest, the atmosphere state uses
pressure-level data (ClimaAtmos interpolates it vertically); the
model-level path is not ported.

CDS credentials come from the `CDSAPI_URL`/`CDSAPI_KEY` environment variables
or a `~/.cdsapirc` file. See
https://cds.climate.copernicus.eu/how-to-api for setup instructions.
"""
module Era5InitialConditions

import CDSAPI
import ClimaComms
import Dates
import NCDatasets
import Scratch
import p7zip_jll

export fetch_era5_initial_conditions

# ============================================================================
# Cache directory
# ============================================================================

const ERA5_CACHE_DIR = Ref{String}("")

function __init__()
    ERA5_CACHE_DIR[] = get(ENV, "CLIMACOUPLER_ERA5_CACHE_DIR") do
        Scratch.@get_scratch!("era5_initial_conditions")
    end
    return nothing
end

"""
    default_cache_dir()

The directory where fetched ERA5 initial conditions are cached.
"""
function default_cache_dir()
    isempty(ERA5_CACHE_DIR[]) &&
        (ERA5_CACHE_DIR[] = Scratch.@get_scratch!("era5_initial_conditions"))
    return ERA5_CACHE_DIR[]
end

# ============================================================================
# File names
# ============================================================================

datestamp(date::Dates.DateTime) = Dates.format(date, Dates.dateformat"yyyymmdd")

raw_filename(date) = "era5_raw_$(datestamp(date))_0000.nc"
sst_filename(date) = "sst_processed_$(datestamp(date))_0000.nc"
sic_filename(date) = "sic_processed_$(datestamp(date))_0000.nc"
land_filename(date) = "era5_land_processed_$(datestamp(date))_0000.nc"
bucket_filename(date) = "era5_bucket_processed_$(datestamp(date))_0000.nc"
albedo_filename(date) = "albedo_processed_$(datestamp(date))_0000.nc"
sentinel_filename(date) = "era5_ic_$(datestamp(date))_done.toml"

output_filenames(date) = [
    raw_filename(date),
    sst_filename(date),
    sic_filename(date),
    land_filename(date),
    bucket_filename(date),
    albedo_filename(date),
]

# ============================================================================
# CDS requests
# ============================================================================

"""
All 37 pressure levels of the ERA5 pressure-level dataset, in hPa.
"""
const PRESSURE_LEVELS = [
    "1", "2", "3", "5", "7", "10", "20", "30", "50", "70",
    "100", "125", "150", "175", "200", "225", "250", "300", "350", "400",
    "450", "500", "550", "600", "650", "700", "750", "775", "800", "825",
    "850", "875", "900", "925", "950", "975", "1000",
]

"""
Pressure-level variables for the atmosphere initial condition, as
(CDS request name => NetCDF short name). From WeatherQuest
`era5_variables.py` (basic + full atmosphere sets).
"""
const PRESSURE_LEVEL_VARIABLES = [
    "geopotential" => "z",
    "temperature" => "t",
    "specific_humidity" => "q",
    "u_component_of_wind" => "u",
    "v_component_of_wind" => "v",
    "vertical_velocity" => "w",
    "specific_cloud_liquid_water_content" => "clwc",
    "specific_cloud_ice_water_content" => "ciwc",
    "specific_rain_water_content" => "crwc",
    "specific_snow_water_content" => "cswc",
]

"""
Instantaneous single-level surface variables, as
(CDS request name => NetCDF short name). From WeatherQuest
`era5_variables.py` (surface set, without accumulated variables).
"""
const SURFACE_VARIABLES = [
    "surface_pressure" => "sp",
    "skin_temperature" => "skt",
    "geopotential" => "z",
    "sea_surface_temperature" => "sst",
    "sea_ice_cover" => "siconc",
    "2m_temperature" => "t2m",
    "mean_sea_level_pressure" => "msl",
    "ice_temperature_layer_1" => "istl1",
    "ice_temperature_layer_2" => "istl2",
    "ice_temperature_layer_3" => "istl3",
    "ice_temperature_layer_4" => "istl4",
]

"""
Single-level land variables, as (CDS request name => NetCDF short name).
These live in the `reanalysis-era5-single-levels` dataset, not in ERA5-Land.
From WeatherQuest `era5_variables.py` (land set).
"""
const LAND_VARIABLES = [
    "volumetric_soil_water_layer_1" => "swvl1",
    "volumetric_soil_water_layer_2" => "swvl2",
    "volumetric_soil_water_layer_3" => "swvl3",
    "volumetric_soil_water_layer_4" => "swvl4",
    "soil_temperature_level_1" => "stl1",
    "soil_temperature_level_2" => "stl2",
    "soil_temperature_level_3" => "stl3",
    "soil_temperature_level_4" => "stl4",
    "snow_depth" => "sd",
    "temperature_of_snow_layer" => "tsn",
    "skin_temperature" => "skt",
    "skin_reservoir_content" => "src",
    "forecast_albedo" => "fal",
    "forecast_surface_roughness" => "fsr",
    "forecast_logarithm_of_surface_roughness_for_heat" => "flsr",
    "leaf_area_index_low_vegetation" => "lai_lv",
    "leaf_area_index_high_vegetation" => "lai_hv",
]

function base_request(date)
    return Dict{String, Any}(
        "product_type" => ["reanalysis"],
        "year" => Dates.format(date, "yyyy"),
        "month" => Dates.format(date, "mm"),
        "day" => Dates.format(date, "dd"),
        "time" => "00:00",
        "area" => [90, -180, -90, 180],
        "data_format" => "netcdf",
    )
end

function pressure_levels_request(date)
    request = base_request(date)
    request["variable"] = first.(PRESSURE_LEVEL_VARIABLES)
    request["pressure_level"] = PRESSURE_LEVELS
    return request
end

function surface_request(date)
    request = base_request(date)
    request["variable"] = first.(SURFACE_VARIABLES)
    return request
end

function land_request(date)
    request = base_request(date)
    request["variable"] = first.(LAND_VARIABLES)
    return request
end

# ============================================================================
# Credentials
# ============================================================================

"""
    cds_credentials_available()

Return `true` when CDS credentials are available through the
`CDSAPI_URL`/`CDSAPI_KEY` environment variables or a `~/.cdsapirc` file.
"""
function cds_credentials_available()
    haskey(ENV, "CDSAPI_URL") && haskey(ENV, "CDSAPI_KEY") && return true
    return isfile(joinpath(homedir(), ".cdsapirc"))
end

function assert_cds_credentials()
    cds_credentials_available() && return nothing
    error(
        "No CDS credentials found. Register at " *
        "https://cds.climate.copernicus.eu, accept the ERA5 licence, and " *
        "either create ~/.cdsapirc or set the CDSAPI_URL and CDSAPI_KEY " *
        "environment variables. See " *
        "https://cds.climate.copernicus.eu/how-to-api for instructions.",
    )
end

# ============================================================================
# Download
# ============================================================================

"""
    unpack_download(path)

Return the path of the NetCDF file for a completed CDS download. When CDS
delivers a zip archive (it can split requests into one file per data
stream), extract it and merge the contained NetCDF files.
"""
function unpack_download(path)
    magic = open(io -> read(io, 2), path)
    magic == UInt8['P', 'K'] || return path
    extract_dir = path * ".extracted"
    mkpath(extract_dir)
    run(pipeline(`$(p7zip_jll.p7zip()) x $path -o$extract_dir -y`; stdout = devnull))
    nc_files = filter(endswith(".nc"), readdir(extract_dir; join = true))
    isempty(nc_files) && error("CDS zip archive $path contains no NetCDF files")
    length(nc_files) == 1 && return only(nc_files)
    merged = path * ".merged.nc"
    merge_netcdf_variables(nc_files, merged)
    return merged
end

"""
    merge_netcdf_variables(paths, output_path)

Merge the variables of several NetCDF files on the same grid into one file.
The first file defines the dimensions; later files add their variables.
"""
function merge_netcdf_variables(paths, output_path)
    NCDatasets.NCDataset(output_path, "c") do ncout
        for (i, path) in enumerate(paths)
            NCDatasets.NCDataset(path) do ncin
                if i == 1
                    for (name, len) in ncin.dim
                        NCDatasets.defDim(ncout, name, len)
                    end
                end
                for (name, var) in ncin
                    haskey(ncout, name) && continue
                    all(d -> haskey(ncout.dim, d), NCDatasets.dimnames(var)) ||
                        continue
                    data = Array(var)
                    if Missing <: eltype(data) &&
                       nonmissingtype(eltype(data)) <: AbstractFloat
                        data = replace(data, missing => NaN)
                    end
                    NCDatasets.defVar(
                        ncout,
                        name,
                        data,
                        NCDatasets.dimnames(var);
                        attrib = clean_attributes(var),
                    )
                end
            end
        end
    end
    return output_path
end

"""
    download_era5(date, download_dir; retrieve_fn, wait)

Submit the three CDS requests for `date` and download the results into
`download_dir`. Return a NamedTuple with the paths of the pressure-level,
surface, and land NetCDF files. `retrieve_fn` has the signature of
`CDSAPI.retrieve` and exists so tests can inject a fake.
"""
function download_era5(
    date,
    download_dir;
    retrieve_fn = CDSAPI.retrieve,
    wait = 30.0,
)
    specs = (
        pressure = ("reanalysis-era5-pressure-levels", pressure_levels_request(date)),
        surface = ("reanalysis-era5-single-levels", surface_request(date)),
        land = ("reanalysis-era5-single-levels", land_request(date)),
    )
    paths = map(keys(specs)) do key
        (dataset, request) = specs[key]
        path = joinpath(download_dir, "cds_$(key)_$(datestamp(date)).nc")
        @info "Requesting ERA5 data from CDS (this can queue for a while)" dataset date group = key
        retrieve_fn(dataset, request, path; wait)
        return unpack_download(path)
    end
    return NamedTuple{keys(specs)}(paths)
end

# ============================================================================
# NetCDF reading helpers
# ============================================================================

const LON_DIM_NAMES = ("longitude", "lon")
const LAT_DIM_NAMES = ("latitude", "lat")
const TIME_DIM_NAMES = ("valid_time", "time")
const LEVEL_DIM_NAMES = ("pressure_level", "level", "plev")

"""
Attributes that must not be copied when re-writing decoded data.
"""
const ENCODING_ATTRIBUTES =
    ("scale_factor", "add_offset", "_FillValue", "missing_value")

function clean_attributes(var)
    return Dict(
        String(k) => v for
        (k, v) in var.attrib if !(String(k) in ENCODING_ATTRIBUTES)
    )
end

function find_dim(dims, candidates, what)
    idx = findfirst(d -> d in candidates, dims)
    isnothing(idx) && error("No $what dimension found among $(dims)")
    return idx
end

"""
    read_surface_field(ds, name)

Read the 2D field `name` and return it as a `(lon, lat)` matrix of
`Union{Missing, Float64}`. A singleton time dimension is dropped.
"""
function read_surface_field(ds, name)
    haskey(ds, name) || error("Variable $name not found in $(NCDatasets.path(ds))")
    var = ds[name]
    dims = NCDatasets.dimnames(var)
    lon_idx = find_dim(dims, LON_DIM_NAMES, "longitude")
    lat_idx = find_dim(dims, LAT_DIM_NAMES, "latitude")
    data = Array(var)
    slices = map(enumerate(dims)) do (i, _)
        (i == lon_idx || i == lat_idx) ? Colon() : 1
    end
    field = data[slices...]
    lon_idx_2d = lon_idx < lat_idx ? 1 : 2
    lon_idx_2d == 1 || (field = permutedims(field))
    return Array{Union{Missing, Float64}}(field)
end

function read_lonlat(ds)
    lon_name = findfirst(n -> haskey(ds, n), collect(LON_DIM_NAMES))
    lat_name = findfirst(n -> haskey(ds, n), collect(LAT_DIM_NAMES))
    (isnothing(lon_name) || isnothing(lat_name)) &&
        error("No longitude/latitude coordinates found in $(NCDatasets.path(ds))")
    lon = Float64.(Array(ds[collect(LON_DIM_NAMES)[lon_name]]))
    lat = Float64.(Array(ds[collect(LAT_DIM_NAMES)[lat_name]]))
    return lon, lat
end

"""
    reference_date(ds, fallback)

The analysis date of a downloaded file, from its time coordinate. Falls back
to `fallback` when no time coordinate is present.
"""
function reference_date(ds, fallback)
    for name in TIME_DIM_NAMES
        haskey(ds, name) || continue
        times = Array(ds[name])
        isempty(times) && continue
        return Dates.DateTime(first(times))
    end
    return fallback
end

# ============================================================================
# Field manipulation helpers
# ============================================================================

"""
    zero_fill(field)

Replace `missing` and `NaN` with 0. Used for masked fields whose consumers
require 0 outside the mask (soil moisture, soil temperature, snow).
"""
zero_fill(field) = replace(x -> (ismissing(x) || isnan(x)) ? 0.0 : Float64(x), field)

"""
    nearest_neighbor_fill(field)

Replace `missing`/`NaN` cells with the value of the nearest valid cell, with
distance measured in index space (a multi-source breadth-first search over
the 8-connected grid). Used for fields that must be globally defined but
have no meaningful fill value (SST over land, snow temperature over ocean).
Port of `fill_nans_nearest_neighbor` (WeatherQuest `interpolate.jl`), with
BFS in place of ScatteredInterpolation.
"""
function nearest_neighbor_fill(field::AbstractMatrix)
    filled = Float64[ismissing(x) ? NaN : Float64(x) for x in field]
    any(isnan, filled) || return filled
    all(isnan, filled) && error("Cannot fill a field with no valid values")
    (nx, ny) = size(filled)
    queue = Vector{Tuple{Int, Int}}()
    visited = falses(nx, ny)
    for j in 1:ny, i in 1:nx
        if !isnan(filled[i, j])
            push!(queue, (i, j))
            visited[i, j] = true
        end
    end
    head = 1
    while head <= length(queue)
        (i, j) = queue[head]
        head += 1
        for dj in -1:1, di in -1:1
            (di == 0 && dj == 0) && continue
            (ii, jj) = (i + di, j + dj)
            (1 <= ii <= nx && 1 <= jj <= ny) || continue
            visited[ii, jj] && continue
            visited[ii, jj] = true
            filled[ii, jj] = filled[i, j]
            push!(queue, (ii, jj))
        end
    end
    return filled
end

"""
    roll_longitudes(lon)

The permutation that sorts longitudes converted to [0, 360), and the sorted
longitudes. SST, SIC, and albedo outputs use a [0, 360) longitude axis.
Port of the longitude roll in WeatherQuest `interpolate.jl`.
"""
function roll_longitudes(lon)
    lon360 = mod.(lon, 360)
    perm = sortperm(lon360)
    return lon360[perm], perm
end

"""
    monthly_time_points(date)

Four monthly time points [date - 1 month, date, date + 1 month,
date + 2 months] as seconds since 1970-01-01. The SST, SIC, and albedo
files replicate one field across these points so time interpolation is
constant over a subseasonal run. Port of the time axis construction in
WeatherQuest `interpolate.jl` (`preprocess_sst`).
"""
function monthly_time_points(date)
    points = [
        date - Dates.Month(1),
        date,
        date + Dates.Month(1),
        date + Dates.Month(2),
    ]
    epoch = Dates.DateTime(1970, 1, 1)
    return [Dates.value(Dates.Second(t - epoch)) for t in points]
end

const TIME_ATTRIB = Dict(
    "standard_name" => "time",
    "long_name" => "time",
    "units" => "seconds since 1970-01-01",
    "calendar" => "proleptic_gregorian",
)

const LON_ATTRIB = Dict(
    "standard_name" => "longitude",
    "long_name" => "Longitude",
    "units" => "degrees_east",
    "axis" => "X",
)

const LAT_ATTRIB = Dict(
    "standard_name" => "latitude",
    "long_name" => "Latitude",
    "units" => "degrees_north",
    "axis" => "Y",
)

"""
    define_lonlat_time!(ncout, lon, lat, time_points)

Define the lon/lat dimensions, and a time dimension when `time_points` is
not `nothing`, with their coordinate variables.
"""
function define_lonlat_time!(ncout, lon, lat, time_points)
    NCDatasets.defDim(ncout, "lon", length(lon))
    NCDatasets.defDim(ncout, "lat", length(lat))
    lon_var =
        NCDatasets.defVar(ncout, "lon", Float32, ("lon",), attrib = LON_ATTRIB)
    lat_var =
        NCDatasets.defVar(ncout, "lat", Float32, ("lat",), attrib = LAT_ATTRIB)
    lon_var[:] = lon
    lat_var[:] = lat
    if !isnothing(time_points)
        NCDatasets.defDim(ncout, "time", length(time_points))
        time_var = NCDatasets.defVar(
            ncout,
            "time",
            Int64,
            ("time",),
            attrib = TIME_ATTRIB,
        )
        time_var[:] = time_points
    end
    return nothing
end

"""
    write_replicated_time_var!(ncout, name, field, ntimes; attrib)

Define the `(lon, lat, time)` variable `name` and write `field` into every
time slice.
"""
function write_replicated_time_var!(ncout, name, field, ntimes; attrib)
    var = NCDatasets.defVar(
        ncout,
        name,
        Float32,
        ("lon", "lat", "time"),
        attrib = attrib,
    )
    for t in 1:ntimes
        var[:, :, t] = field
    end
    return nothing
end

# ============================================================================
# Soil layer constants (ERA5 soil discretization, values as in WeatherQuest
# interpolate.jl: interpolate_land and interpolate_bucket)
# ============================================================================

"""
ERA5 soil layer midpoint depths [m].
"""
const SOIL_LAYER_MIDPOINTS = [0.035, 0.175, 0.64, 1.945]

"""
ERA5 soil layer top depths [m].
"""
const SOIL_LAYER_TOPS = [0.0, 0.07, 0.28, 1.00]

"""
ERA5 soil layer bottom depths [m].
"""
const SOIL_LAYER_BOTTOMS = [0.07, 0.28, 1.00, 2.89]

"""
The `z` coordinate of the 3D soil outputs: negative depths sorted
increasing, so index 1 is the deepest layer.
"""
const SOIL_Z = -reverse(SOIL_LAYER_MIDPOINTS)

"""
Reverse the layer dimension of a `(lon, lat, layer)` array so it matches the
`SOIL_Z` ordering (deepest first).
"""
soil_layers_to_z(field) = reverse(field; dims = 3)

const SOIL_Z_ATTRIB = Dict(
    "standard_name" => "depth",
    "long_name" => "Soil depth (negative below surface)",
    "units" => "m",
    "axis" => "Z",
)

# ============================================================================
# Preprocessing: SST and SIC
# ============================================================================

"""
    process_sst(surface_path, output_path, date)

Write the SST file: variable `SST` in Celsius on `(lon, lat, time)`, with
land filled by nearest neighbor, a [0, 360) longitude axis, and the field
replicated over four monthly time points.
Port of `preprocess_sst` (WeatherQuest `interpolate.jl`).
"""
function process_sst(surface_path, output_path, date)
    NCDatasets.NCDataset(surface_path) do ncin
        lon, lat = read_lonlat(ncin)
        ref_date = reference_date(ncin, date)
        sst = nearest_neighbor_fill(read_surface_field(ncin, "sst")) .- 273.15
        lon360, perm = roll_longitudes(lon)
        sst = sst[perm, :]
        time_points = monthly_time_points(ref_date)
        NCDatasets.NCDataset(output_path, "c") do ncout
            define_lonlat_time!(ncout, lon360, lat, time_points)
            write_replicated_time_var!(
                ncout,
                "SST",
                sst,
                length(time_points);
                attrib = Dict(
                    "standard_name" => "sea_surface_temperature",
                    "long_name" => "Sea Surface Temperature",
                    "units" => "celsius",
                    "varname" => "SST",
                ),
            )
        end
    end
    return output_path
end

"""
    process_sic(surface_path, output_path, date)

Write the sea ice file: variable `SEAICE` in percent on `(lon, lat, time)`
(missing treated as 0), plus nearest-neighbor-filled ice temperature layers
`ISTL1..4` in Kelvin.
Port of `preprocess_sic` (WeatherQuest `interpolate.jl`).
"""
function process_sic(surface_path, output_path, date)
    NCDatasets.NCDataset(surface_path) do ncin
        lon, lat = read_lonlat(ncin)
        ref_date = reference_date(ncin, date)
        sic = clamp.(zero_fill(read_surface_field(ncin, "siconc")) .* 100, 0, 100)
        lon360, perm = roll_longitudes(lon)
        sic = sic[perm, :]
        istl = map(1:4) do k
            name = "istl$k"
            haskey(ncin, name) || return nothing
            nearest_neighbor_fill(read_surface_field(ncin, name))[perm, :]
        end
        time_points = monthly_time_points(ref_date)
        NCDatasets.NCDataset(output_path, "c") do ncout
            define_lonlat_time!(ncout, lon360, lat, time_points)
            write_replicated_time_var!(
                ncout,
                "SEAICE",
                sic,
                length(time_points);
                attrib = Dict(
                    "standard_name" => "sea_ice_cover",
                    "long_name" => "Sea Ice Concentration",
                    "units" => "%",
                    "varname" => "SEAICE",
                ),
            )
            for (k, field) in enumerate(istl)
                isnothing(field) && continue
                write_replicated_time_var!(
                    ncout,
                    "ISTL$k",
                    field,
                    length(time_points);
                    attrib = Dict(
                        "long_name" => "Ice temperature layer $k",
                        "units" => "K",
                        "varname" => "ISTL$k",
                    ),
                )
            end
        end
    end
    return output_path
end

# ============================================================================
# Preprocessing: integrated land
# ============================================================================

"""
    process_land(land_path, output_path, date)

Write the integrated-land initial condition file. 2D variables `skt`, `tsn`
(Kelvin, nearest-neighbor filled), `swe` (m of water equivalent), and `lai`;
3D variables `swvl` (total volumetric water), `stl` (Kelvin), `si`
(volumetric ice), and `sie` (volumetric internal energy) on the negative
depth coordinate `z`. Ocean points are 0 (the ClimaLand reader masks with
`> 0`), and the latitude axis is sorted increasing.
Port of `interpolate_land` (WeatherQuest `interpolate.jl`), including the
binary freeze partition and internal energy formulas.
"""
function process_land(land_path, output_path, date)
    NCDatasets.NCDataset(land_path) do ncin
        lon, lat = read_lonlat(ncin)
        lat_perm = sortperm(lat)
        nlayers = length(SOIL_LAYER_MIDPOINTS)

        read_soil = prefix -> begin
            layers = map(1:nlayers) do k
                zero_fill(read_surface_field(ncin, "$prefix$k"))[:, lat_perm]
            end
            cat(layers...; dims = 3)
        end
        swvl = read_soil("swvl")
        stl = read_soil("stl")

        swe = zero_fill(read_surface_field(ncin, "sd"))[:, lat_perm]
        tsn = nearest_neighbor_fill(read_surface_field(ncin, "tsn"))[:, lat_perm]
        skt = nearest_neighbor_fill(read_surface_field(ncin, "skt"))[:, lat_perm]

        lai_lv = read_surface_field(ncin, "lai_lv")
        lai_hv = read_surface_field(ncin, "lai_hv")
        lai = max.(nearest_neighbor_fill(lai_lv .+ lai_hv), 0.0)[:, lat_perm]

        # Frozen soil water partition: below freezing, all water is ice.
        # The internal energy follows the WeatherQuest/ClimaLand convention
        # with soil porosity 0.4 and mineral soil density 2650 kg/m^3.
        T_0 = 273.15
        rho_i = 917.0
        L_f_0 = 3.34e5
        si = ifelse.(stl .< T_0, swvl, 0.0)
        theta_l = swvl .- si
        rho_s = (1 - 0.4) * 2650 .+ theta_l .* 1000 .+ si .* rho_i
        c_s =
            (
                (1 - 0.4) * 2650 * 800 .+ theta_l .* 1000 .* 4186 .+
                si .* rho_i .* 2100
            ) ./ rho_s
        sie = rho_s .* c_s .* (stl .- T_0) .- si .* rho_i .* L_f_0

        NCDatasets.NCDataset(output_path, "c") do ncout
            define_lonlat_time!(ncout, lon, lat[lat_perm], nothing)
            NCDatasets.defDim(ncout, "z", nlayers)
            z_var = NCDatasets.defVar(
                ncout,
                "z",
                Float32,
                ("z",),
                attrib = SOIL_Z_ATTRIB,
            )
            z_var[:] = SOIL_Z

            for (name, field, units, long_name) in (
                ("swvl", swvl, "m^3/m^3", "Volumetric fraction of water"),
                ("stl", stl, "K", "Soil temperature"),
                ("si", si, "m^3/m^3", "Volumetric fraction of ice"),
                ("sie", sie, "J/m^3", "Soil volumetric internal energy"),
            )
                var = NCDatasets.defVar(
                    ncout,
                    name,
                    Float32,
                    ("lon", "lat", "z"),
                    attrib = Dict(
                        "units" => units,
                        "longname" => long_name,
                        "varname" => name,
                    ),
                )
                var[:, :, :] = soil_layers_to_z(field)
            end

            for (name, field, units, long_name) in (
                ("swe", swe, "m", "Snow water equivalent"),
                ("tsn", tsn, "K", "Temperature of snow layer"),
                ("skt", skt, "K", "Skin temperature"),
                ("lai", lai, "m^2 m^-2", "Leaf area index"),
            )
                var = NCDatasets.defVar(
                    ncout,
                    name,
                    Float32,
                    ("lon", "lat"),
                    attrib = Dict(
                        "units" => units,
                        "longname" => long_name,
                        "varname" => name,
                    ),
                )
                var[:, :] = field
            end
        end
    end
    return output_path
end

# ============================================================================
# Preprocessing: bucket land
# ============================================================================

"""
    process_bucket(land_path, output_path, date; subsurface_water_z_max = 0.5)

Write the bucket initial condition file. 2D variables `W` (soil water
column integrated down to `subsurface_water_z_max`), `Ws` (skin reservoir),
`S` (snow water equivalent), `tsn`, and `skt`; 3D variable `T` (soil
temperature, nearest-neighbor filled over ocean) on the negative depth
coordinate `z`.
Port of `interpolate_bucket` (WeatherQuest `interpolate.jl`), including the
layer-thickness weights of the `W` column integral.
"""
function process_bucket(
    land_path,
    output_path,
    date;
    subsurface_water_z_max = 0.5,
)
    NCDatasets.NCDataset(land_path) do ncin
        lon, lat = read_lonlat(ncin)
        nlayers = length(SOIL_LAYER_MIDPOINTS)

        W = zeros(length(lon), length(lat))
        for k in 1:nlayers
            thickness = max(
                0.0,
                min(subsurface_water_z_max, SOIL_LAYER_BOTTOMS[k]) -
                SOIL_LAYER_TOPS[k],
            )
            thickness > 0 || continue
            W .+= zero_fill(read_surface_field(ncin, "swvl$k")) .* thickness
        end

        Ws = zero_fill(read_surface_field(ncin, "src"))
        S = zero_fill(read_surface_field(ncin, "sd"))

        # Soil temperature is filled over ocean so the bucket has values
        # everywhere. Zeros mark masked points in the source data.
        T_layers = map(1:nlayers) do k
            field = read_surface_field(ncin, "stl$k")
            with_nan = Union{Missing, Float64}[
                (ismissing(x) || x == 0) ? missing : x for x in field
            ]
            nearest_neighbor_fill(with_nan)
        end
        T = cat(T_layers...; dims = 3)

        tsn = nearest_neighbor_fill(read_surface_field(ncin, "tsn"))
        skt = nearest_neighbor_fill(read_surface_field(ncin, "skt"))

        NCDatasets.NCDataset(output_path, "c") do ncout
            define_lonlat_time!(ncout, lon, lat, nothing)
            NCDatasets.defDim(ncout, "z", nlayers)
            z_var = NCDatasets.defVar(
                ncout,
                "z",
                Float32,
                ("z",),
                attrib = SOIL_Z_ATTRIB,
            )
            z_var[:] = SOIL_Z

            T_var = NCDatasets.defVar(
                ncout,
                "T",
                Float32,
                ("lon", "lat", "z"),
                attrib = Dict(
                    "units" => "K",
                    "longname" => "Soil temperature profile",
                    "varname" => "T",
                ),
            )
            T_var[:, :, :] = soil_layers_to_z(T)

            for (name, field, units, long_name) in (
                ("W", W, "m", "Subsurface water content"),
                ("Ws", Ws, "m", "Surface water content"),
                ("S", S, "m", "Snow water equivalent"),
                ("tsn", tsn, "K", "Temperature of snow layer"),
                ("skt", skt, "K", "Skin temperature"),
            )
                var = NCDatasets.defVar(
                    ncout,
                    name,
                    Float32,
                    ("lon", "lat"),
                    attrib = Dict(
                        "units" => units,
                        "longname" => long_name,
                        "varname" => name,
                    ),
                )
                var[:, :] = field
            end
        end
    end
    return output_path
end

# ============================================================================
# Preprocessing: albedo
# ============================================================================

"""
    process_albedo(land_path, output_path, date)

Write the albedo file: variable `sw_alb_clr` (ERA5 forecast albedo) on
`(lon, lat, time)`, plus the roughness fields `fsr` and `flsr` when present.
Port of `preprocess_albedo` (WeatherQuest `interpolate.jl`).
"""
function process_albedo(land_path, output_path, date)
    NCDatasets.NCDataset(land_path) do ncin
        lon, lat = read_lonlat(ncin)
        ref_date = reference_date(ncin, date)
        lon360, perm = roll_longitudes(lon)
        albedo = nearest_neighbor_fill(read_surface_field(ncin, "fal"))[perm, :]
        time_points = monthly_time_points(ref_date)
        NCDatasets.NCDataset(output_path, "c") do ncout
            define_lonlat_time!(ncout, lon360, lat, time_points)
            write_replicated_time_var!(
                ncout,
                "sw_alb_clr",
                albedo,
                length(time_points);
                attrib = Dict(
                    "standard_name" => "surface_albedo",
                    "long_name" => "Clear-sky shortwave albedo",
                    "units" => "1",
                    "varname" => "sw_alb_clr",
                ),
            )
            for (name, units, long_name) in (
                ("fsr", "m", "Forecast surface roughness"),
                (
                    "flsr",
                    "1",
                    "Forecast logarithm of surface roughness for heat",
                ),
            )
                haskey(ncin, name) || continue
                field = nearest_neighbor_fill(read_surface_field(ncin, name))[perm, :]
                write_replicated_time_var!(
                    ncout,
                    name,
                    field,
                    length(time_points);
                    attrib = Dict(
                        "long_name" => long_name,
                        "units" => units,
                        "varname" => name,
                    ),
                )
            end
        end
    end
    return output_path
end

# ============================================================================
# Preprocessing: raw atmosphere file
# ============================================================================

"""
    build_raw(pressure_path, surface_path, output_path)

Build the raw atmosphere file that ClimaAtmos's WeatherModel fallback path
reads: the pressure-level variables merged with the single-level `skt` and
`sp`, and the single-level geopotential renamed to `surface_geopotential`.
The dimension names `longitude`, `latitude`, `pressure_level`, and
`valid_time` are preserved because ClimaAtmos asserts them.
Port of `combine_era5_datasets` (WeatherQuest `get_initial_conditions.py`),
including the geopotential rename.
"""
function build_raw(pressure_path, surface_path, output_path)
    NCDatasets.NCDataset(output_path, "c") do ncout
        NCDatasets.NCDataset(pressure_path) do ncp
            for name in ("longitude", "latitude", "pressure_level", "valid_time")
                haskey(ncp.dim, name) || error(
                    "Expected dimension $name in the CDS pressure-level " *
                    "file, found $(collect(keys(ncp.dim))). The CDS output " *
                    "format may have changed.",
                )
            end
            for (name, len) in ncp.dim
                NCDatasets.defDim(ncout, name, len)
            end
            for (name, var) in ncp
                data = Array(var)
                attrib = clean_attributes(var)
                if eltype(data) <: Union{Missing, Number} &&
                   !(eltype(data) <: Union{Missing, Integer})
                    data = Float32.(replace(data, missing => NaN))
                end
                NCDatasets.defVar(
                    ncout,
                    name,
                    data,
                    NCDatasets.dimnames(var);
                    attrib = attrib,
                )
            end
        end
        NCDatasets.NCDataset(surface_path) do ncs
            for (src_name, dst_name) in
                (("skt", "skt"), ("sp", "sp"), ("z", "surface_geopotential"))
                haskey(ncs, src_name) || error(
                    "Variable $src_name not found in the CDS surface file",
                )
                var = ncs[src_name]
                dims = NCDatasets.dimnames(var)
                all(d -> haskey(ncout.dim, d), dims) || error(
                    "The CDS surface file grid does not match the " *
                    "pressure-level file grid (dims $(dims))",
                )
                data = Float32.(replace(Array(var), missing => NaN))
                attrib = clean_attributes(var)
                attrib["varname"] = dst_name
                NCDatasets.defVar(ncout, dst_name, data, dims; attrib = attrib)
            end
        end
    end
    return output_path
end

# ============================================================================
# Validation
# ============================================================================

struct OutputSpec
    filename_fn::Function
    required_vars::Vector{String}
    check::Function
end

function check_no_nan(ds, vars, filename)
    for name in vars
        data = Array(ds[name])
        if any(x -> ismissing(x) || (x isa AbstractFloat && isnan(x)), data)
            error("Variable $name in $filename contains missing or NaN values")
        end
    end
    return nothing
end

const OUTPUT_SPECS = [
    OutputSpec(
        raw_filename,
        ["u", "v", "w", "t", "q", "z", "skt", "sp", "surface_geopotential"],
        (ds, filename) -> begin
            for dim in ("longitude", "latitude", "pressure_level", "valid_time")
                haskey(ds.dim, dim) ||
                    error("Missing dimension $dim in $filename")
            end
            check_no_nan(
                ds,
                ["u", "v", "t", "q", "z", "skt", "sp", "surface_geopotential"],
                filename,
            )
        end,
    ),
    OutputSpec(
        sst_filename,
        ["SST"],
        (ds, filename) -> begin
            check_no_nan(ds, ["SST"], filename)
            sst = Array(ds["SST"])
            (minimum(sst) > -60 && maximum(sst) < 60) ||
                error("SST in $filename is outside the plausible range in Celsius")
        end,
    ),
    OutputSpec(
        sic_filename,
        ["SEAICE"],
        (ds, filename) -> begin
            check_no_nan(ds, ["SEAICE"], filename)
            sic = Array(ds["SEAICE"])
            (minimum(sic) >= 0 && maximum(sic) <= 100) ||
                error("SEAICE in $filename is outside [0, 100] percent")
        end,
    ),
    OutputSpec(
        land_filename,
        ["skt", "tsn", "swe", "swvl", "stl", "si", "sie"],
        (ds, filename) -> begin
            check_no_nan(ds, ["skt", "tsn", "swe", "swvl", "stl"], filename)
            stl = Array(ds["stl"])
            all(x -> x == 0 || x > 100, stl) || error(
                "stl in $filename has values that are neither 0 (ocean) " *
                "nor a plausible temperature in Kelvin",
            )
        end,
    ),
    OutputSpec(
        bucket_filename,
        ["W", "Ws", "S", "T", "tsn", "skt"],
        (ds, filename) -> begin
            check_no_nan(ds, ["W", "Ws", "S", "T"], filename)
            T = Array(ds["T"])
            minimum(T) > 100 ||
                error("Bucket T in $filename has implausibly cold values")
        end,
    ),
    OutputSpec(
        albedo_filename,
        ["sw_alb_clr"],
        (ds, filename) -> begin
            check_no_nan(ds, ["sw_alb_clr"], filename)
            albedo = Array(ds["sw_alb_clr"])
            (minimum(albedo) >= 0 && maximum(albedo) <= 1) ||
                error("sw_alb_clr in $filename is outside [0, 1]")
        end,
    ),
]

"""
    validate_era5_dir(dir, date)

Check that `dir` holds a complete, well-formed set of ERA5 initial
condition files for `date`. Throws a descriptive error on the first
problem found.
"""
function validate_era5_dir(dir, date)
    for spec in OUTPUT_SPECS
        filename = spec.filename_fn(date)
        path = joinpath(dir, filename)
        isfile(path) || error("Missing ERA5 initial condition file $path")
        NCDatasets.NCDataset(path) do ds
            for name in spec.required_vars
                haskey(ds, name) ||
                    error("Missing variable $name in $filename")
            end
            spec.check(ds, filename)
        end
    end
    return nothing
end

# ============================================================================
# Cache management and orchestration
# ============================================================================

"""
    era5_files_complete(dir, date)

Whether the cache at `dir` holds all output files and the completion
sentinel for `date`.
"""
function era5_files_complete(dir, date)
    isfile(joinpath(dir, sentinel_filename(date))) || return false
    return all(name -> isfile(joinpath(dir, name)), output_filenames(date))
end

function write_sentinel(dir, date)
    version = string(pkgversion(parentmodule(@__MODULE__)))
    open(joinpath(dir, sentinel_filename(date)), "w") do io
        println(io, "date = \"$(date)\"")
        println(io, "created = \"$(Dates.now())\"")
        println(io, "climacoupler_version = \"$version\"")
    end
    return nothing
end

function remove_cached_files(dir, date)
    for name in [output_filenames(date); sentinel_filename(date)]
        path = joinpath(dir, name)
        isfile(path) && rm(path)
    end
    return nothing
end

"""
    cleanup_stale_tmpdirs(dir; max_age_seconds = 86400)

Remove leftover download directories from fetches that died without
cleanup. Only directories older than `max_age_seconds` are removed, so a
concurrent fetch into the same cache keeps its live download directory.
"""
function cleanup_stale_tmpdirs(dir; max_age_seconds = 86400)
    for name in readdir(dir)
        path = joinpath(dir, name)
        (startswith(name, "tmp_era5_") && isdir(path)) || continue
        (time() - mtime(path)) > max_age_seconds || continue
        rm(path; recursive = true)
    end
    return nothing
end

"""
    fetch_era5_initial_conditions(start_date; dir, force, comms_ctx, retrieve_fn, wait)

Return a directory that holds the complete set of ERA5 initial condition
files for `start_date`, downloading and preprocessing them from the CDS API
when they are not already cached.

# Arguments
- `start_date`: the simulation start date (files are named from its day at
  00:00 UTC).
- `dir`: the cache directory (default: a Scratch.jl directory, override
  with the `CLIMACOUPLER_ERA5_CACHE_DIR` environment variable).
- `force`: delete any cached files for this date and fetch again.
- `comms_ctx`: an optional `ClimaComms` context. Only the root process
  downloads; the others wait on a barrier.
- `retrieve_fn`: the retrieval function (`CDSAPI.retrieve` by default);
  tests inject a fake here.
- `wait`: seconds between CDS job status checks.

Downloads land in a temporary directory inside `dir` and files move into
place only after validation, so an interrupted fetch never leaves a
partial cache.
"""
function fetch_era5_initial_conditions(
    start_date;
    dir = default_cache_dir(),
    force = false,
    comms_ctx = nothing,
    retrieve_fn = CDSAPI.retrieve,
    wait = 30.0,
)
    date = Dates.DateTime(Dates.Date(start_date))
    mkpath(dir)

    if !isnothing(comms_ctx) && !ClimaComms.iamroot(comms_ctx)
        ClimaComms.barrier(comms_ctx)
        era5_files_complete(dir, date) ||
            error("ERA5 initial condition fetch failed on the root process")
        return dir
    end

    try
        force && remove_cached_files(dir, date)
        if era5_files_complete(dir, date)
            @info "Using cached ERA5 initial conditions" dir date
            return dir
        end
        # Tests inject a fake retrieve_fn and need no credentials
        retrieve_fn === CDSAPI.retrieve && assert_cds_credentials()
        cleanup_stale_tmpdirs(dir)
        mktempdir(dir; prefix = "tmp_era5_") do tmpdir
            files = download_era5(date, tmpdir; retrieve_fn, wait)
            @info "Preprocessing ERA5 initial conditions" date
            build_raw(files.pressure, files.surface, joinpath(tmpdir, raw_filename(date)))
            process_sst(files.surface, joinpath(tmpdir, sst_filename(date)), date)
            process_sic(files.surface, joinpath(tmpdir, sic_filename(date)), date)
            process_land(files.land, joinpath(tmpdir, land_filename(date)), date)
            process_bucket(files.land, joinpath(tmpdir, bucket_filename(date)), date)
            process_albedo(files.land, joinpath(tmpdir, albedo_filename(date)), date)
            validate_era5_dir(tmpdir, date)
            for name in output_filenames(date)
                mv(joinpath(tmpdir, name), joinpath(dir, name); force = true)
            end
            write_sentinel(dir, date)
        end
        @info "ERA5 initial conditions ready" dir date
        return dir
    finally
        isnothing(comms_ctx) || ClimaComms.barrier(comms_ctx)
    end
end

end # module Era5InitialConditions
