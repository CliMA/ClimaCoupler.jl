
using Statistics
using ClimaCoupler.Regridder: remap_field_cgll_to_rll

# data types for postprocessing
abstract type PostProcessedData end

struct ZonalMeanData <: PostProcessedData
    data_2d::Array
    name::String
    pressure_interpol::Bool
    lev::Array
    lat::Array
end
function ZonalMeanData(post_data, name; pressure_interpol = false)
    ZonalMeanData(
        getproperty(post_data, :data),
        name,
        pressure_interpol,
        getproperty(post_data, :lev),
        getproperty(post_data, :lat),
    )
end

struct HorizontalSliceData <: PostProcessedData
    data_2d::Array
    name::String
    lat::Array
    lon::Array
end
function HorizontalSliceData(post_data, name)
    HorizontalSliceData(getproperty(post_data, :data), name, getproperty(post_data, :lat), getproperty(post_data, :lon))
end

"""
    postprocess(
        name,
        raw_data::Union{Fields.Field, NamedTuple},
        p_methods;
        lev_slice = 1,
        datafile_latlon = nothing,
    ) 

Coordinates regridding, averaging or slicing of variable `name` corresponding 
to `raw_data`. Postprocessing methods are specified in `p_methods`.
"""

function postprocess(
    name,
    raw_data::Union{Fields.Field, NamedTuple},
    p_methods;
    lev_slice = 1,
    datafile_latlon = nothing,
) # TODO: this should be extended to acount for different vartical spaces
    if (:regridded_3d in p_methods) || (:regridded_2d in p_methods)
        DIR = joinpath(REGRID_DIR, "cgll2rll")
        isdir(DIR) ? nothing : mkpath(DIR)
        datafile_latlon =
            (datafile_latlon == nothing) ? datafile_latlon = DIR * "/remapped_" * string(name) * ".nc" : datafile_latlon
        remap_field_cgll_to_rll(name, raw_data, DIR, datafile_latlon)
        grid_data = read_remapped_field(name, datafile_latlon)
    else
        @info "No regridding required."
        grid_data = raw_data
    end

    if :zonal_mean in p_methods
        post_data = ZonalMeanData(
            (;
                data = Statistics.mean(Array(FT.(getproperty(grid_data, :data))), dims = 1)[1, :, :],
                lat = getproperty(grid_data, :lat),
                lev = getproperty(grid_data, :lev),
            ),
            string(name),
        )
    end
    if :horizontal_2d in p_methods
        if length(size(getproperty(grid_data, :data))) == 2
            post_data = HorizontalSliceData(
                (;
                    data = getproperty(grid_data, :data),
                    lat = getproperty(grid_data, :lat),
                    lon = getproperty(grid_data, :lon),
                ),
                string(name),
            )

        else
            post_data = HorizontalSliceData(
                (;
                    data = getproperty(grid_data, :data)[:, :, lev_slice],
                    lat = getproperty(grid_data, :lat),
                    lon = getproperty(grid_data, :lon),
                    lev = getproperty(grid_data, :lev)[lev_slice],
                ),
                string(name),
            )
        end
    end
    return post_data
end

function read_remapped_field(name, datafile_latlon, lev_name = "z")
    nc = NCDataset(datafile_latlon, "r")
    lon = nc["lon"][:]
    lat = nc["lat"][:]
    lev = lev_name in keys(nc) ? nc[lev_name][:] : Float64(-999)
    var = nc[name][:]
    close(nc)

    (; lon = lon, lat = lat, lev = lev, data = var)
end
