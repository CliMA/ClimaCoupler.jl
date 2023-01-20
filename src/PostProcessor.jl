"""
    PostProcessor

This module contains functions for postprocessing model data (saved during the simulation by `Diagnostics.jl`) after the simulation is complete.
"""
module PostProcessor

export PostProcessedData, ZLatLonData, ZLatData, LatLonData, LatData, RawData, DataPackage, postprocess

using Statistics
using NCDatasets: NCDataset

using ClimaCoupler.Regridder: remap_field_cgll_to_rll
using ClimaCore: Fields

# data types for postprocessing
abstract type PostProcessedData end

struct ZLatLonData <: PostProcessedData end
struct ZLatData <: PostProcessedData end
struct LatLonData <: PostProcessedData end
struct LatData <: PostProcessedData end
struct RawData <: PostProcessedData end

"""
    DataPackage(tag::PostProcessedData, name::Symbol, data::Union{Array, Field}; coords = coords)

A container for storing the tyoe, name, data and coordinates of a variable.
"""
struct DataPackage{PPD <: PostProcessedData, NTC <: NamedTuple} # TODO: add long name (ppp info) and units
    tag::PPD
    name::String
    data::Union{Fields.Field, Array}
    coords::NTC
end
function DataPackage(tag::ZLatLonData, name::Symbol, data::Array; coords = coords)
    DataPackage(tag, string(name), data, (; lev = coords.lev, lat = coords.lat, lon = coords.lon))
end
function DataPackage(tag::ZLatData, name::Symbol, data::Array; coords = coords)
    DataPackage(tag, string(name), data, (; lev = coords.lev, lat = coords.lat))
end
function DataPackage(tag::LatLonData, name::Symbol, data::Array; coords = coords)
    DataPackage(tag, string(name), data, (; lat = coords.lat, lon = coords.lon))
end
function DataPackage(tag::LatData, name::Symbol, data::Array; coords = coords)
    DataPackage(tag, string(name), data, (; lat = coords.lat))
end
function DataPackage(tag::RawData, name::Symbol, data::Union{Array, Fields.Field}; coords = nothing)
    DataPackage(tag, string(name), data, (;))
end

"""
    postprocess(
        name::Symbol,
        raw_data::Union{Fields.Field, Array},
        p_methods::Tuple;
        lev_slice = 1,
        datafile_latlon = nothing,
        nlat = 90,
        nlon = 180,
    )

Coordinates regridding, averaging or slicing of variable `name` corresponding
to `raw_data`. Postprocessing methods are specified in `p_methods`. `raw_data` is
assumed to be a `Field` (dimensions corresponding to the model's CGLL grid), a 2D
Array with [longitude, latitude] or a 3D Array [longitude, latitude, level].

# Arguments:
- `name`: [Symbol] variable name
- `raw_data`: [Union{Fields.Field, Array}] variable data
- `p_methods`: [Tuple] postproessing methods (`:regrid`, `:horizontal_slice`, `:zonal_mean`)
- `lev_slice`: [Int] level index along which the `:horizontal_slice` is applied
- `datafile_latlon`: [String] name of the regrid file
- `nlat`: [Int] number of latitudes of the regridded array
- `nlon`: [Symbol] number of longitudes of the regridded array

"""
function postprocess(
    name::Symbol,
    raw_data::Union{Fields.Field, Array},
    p_methods::Tuple;
    lev_slice = 1,
    datafile_latlon = nothing,
    nlat = 90,
    nlon = 180,
    REGRID_DIR = "posrprocess_regrid_tmp/",
    coords = (;),
    raw_tag = RawData(),
)
    # regridding
    if :regrid in p_methods
        DIR = joinpath(REGRID_DIR, "cgll2rll")
        isdir(DIR) ? nothing : mkpath(DIR)
        datafile_latlon =
            (datafile_latlon == nothing) ? datafile_latlon = DIR * "/remapped_" * string(name) * ".nc" : datafile_latlon
        remap_field_cgll_to_rll(name, raw_data, DIR, datafile_latlon, nlat = nlat, nlon = nlon)
        new_data, coords = read_remapped_field(name, datafile_latlon)
        raw_tag = length(size(new_data)) == 3 ? ZLatLonData() : LatLonData()
        package = DataPackage(raw_tag, name, new_data, coords = coords)
    else
        @info "No regridding required."
        package = DataPackage(raw_tag, name, raw_data, coords = coords)
    end

    # spatial slicing and averaging
    if :horizontal_slice in p_methods
        package.tag == RawData() ? @error("Cannot perform horizontal slicing on raw model data. Specify :regrid") :
        nothing
        if package.tag == ZLatLonData()
            package = DataPackage(ZLatData(), name, package.data[:, :, lev_slice], coords = package.coords)
        end
    end

    if :zonal_mean in p_methods
        package.tag == RawData() ? @error("Cannot perform zonal mean on raw model data. Specify :regrid") : nothing
        if package.tag == ZLatLonData()
            package = DataPackage(
                ZLatData(),
                name,
                dropdims(Statistics.mean(package.data, dims = 1), dims = 1),
                coords = package.coords,
            )
        else
            package.tag == LatLonData()
            package = DataPackage(
                LatData(),
                name,
                dropdims(Statistics.mean(package.data, dims = 1), dims = 1),
                coords = (; lat = package.coords.lat),
            )
        end
    end

    return package
end


"""
    read_remapped_field(name::Symbol, datafile_latlon::String, lev_name = "z")

Extract data and coordinates from `datafile_latlon`.
"""
function read_remapped_field(name::Symbol, datafile_latlon::String, lev_name = "z")
    out = NCDataset(datafile_latlon, "r") do nc
        lon = nc["lon"][:]
        lat = nc["lat"][:]
        lev = lev_name in keys(nc) ? nc[lev_name][:] : Float64(-999)
        var = nc[name][:]
        coords = (; lon = lon, lat = lat, lev = lev)

        (var, coords)
    end

    return out
end

end # module
