"""
    PostProcessor

This module contains functions for postprocessing model data (saved during the simulation by `Diagnostics.jl`) after the simulation is complete.
"""
module PostProcessor

export PostProcessedData, ZLatLonData, ZLatData, LatLonData, LatData, RawData, DataPackage, postprocess

using Statistics
using NCDatasets: NCDataset

using ClimaCoupler: Regridder
using ClimaCore: Fields

# data types for postprocessing
"""
    PostProcessedData

Abstract type for postprocessed data.
"""
abstract type PostProcessedData end

"""
    ZLatLonData <: PostProcessedData
Concrete type for 3D data.
"""
struct ZLatLonData <: PostProcessedData end

"""
    ZLatData <: PostProcessedData
Concrete type for 2D data with latitude and level.
"""
struct ZLatData <: PostProcessedData end

"""
    LatLonData <: PostProcessedData
Concrete type for 2D data with latitude and longitude.
"""
struct LatLonData <: PostProcessedData end

"""
    LatData <: PostProcessedData
Concrete type for 1D data with latitude.
"""
struct LatData <: PostProcessedData end

"""
    RawData <: PostProcessedData
Concrete type for raw model data.
"""
struct RawData <: PostProcessedData end

"""
    DataPackage(tag::PostProcessedData, name::Symbol, data::Union{Array, Field}; coords = coords)

A container for storing the tyoe, name, data and coordinates of a variable.
"""
struct DataPackage{PPD <: PostProcessedData, NTC <: NamedTuple, A <: Union{Fields.Field, AbstractArray}} # TODO: add long name (ppp info) and units
    tag::PPD
    name::String
    data::A
    coords::NTC
end
function DataPackage(tag::ZLatLonData, name::Symbol, data::AbstractArray; coords = coords)
    DataPackage(tag, string(name), data, (; lev = coords.lev, lat = coords.lat, lon = coords.lon))
end
function DataPackage(tag::ZLatData, name::Symbol, data::AbstractArray; coords = coords)
    DataPackage(tag, string(name), data, (; lev = coords.lev, lat = coords.lat))
end
function DataPackage(tag::LatLonData, name::Symbol, data::AbstractArray; coords = coords)
    DataPackage(tag, string(name), data, (; lat = coords.lat, lon = coords.lon))
end
function DataPackage(tag::LatData, name::Symbol, data::AbstractArray; coords = coords)
    DataPackage(tag, string(name), data, (; lat = coords.lat))
end
function DataPackage(tag::RawData, name::Symbol, data::Union{AbstractArray, Fields.Field}; coords = nothing)
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
    REGRID_DIR = "postprocess_regrid_tmp/",
    coords = (;),
    raw_tag = RawData(),
)
    # regridding
    if :regrid in p_methods
        DIR = joinpath(REGRID_DIR, "cgll2rll")
        isdir(DIR) ? nothing : mkpath(DIR)
        datafile_latlon =
            (datafile_latlon == nothing) ? datafile_latlon = DIR * "/remapped_" * string(name) * ".nc" : datafile_latlon
        Regridder.remap_field_cgll_to_rll(name, raw_data, DIR, datafile_latlon, nlat = nlat, nlon = nlon)
        new_data, coords = Regridder.read_remapped_field(name, datafile_latlon)
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


end # module
