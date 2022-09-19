
using Statistics

abstract type OfflineDataType end

struct PostProcessedData <: OfflineDataType
    info::NamedTuple
end

struct ZonalMeanData <: OfflineDataType
    data_2d::Array
    clims::Tuple
    title::String
    pressure_interpol::Bool
    lev::Array
    lat::Array
end
function ZonalMeanData(post_data, name; clims = nothing, pressure_interpol = false)
    info = post_data.info
    clims = clims == nothing ? extrema(getproperty(info, :data)) : clims
    ZonalMeanData(
        getproperty(info, :data),
        clims,
        name,
        pressure_interpol,
        getproperty(info, :lev),
        getproperty(info, :lat),
    )
end

struct HorizonalData <: OfflineDataType
    data_2d::Array
    clims::Tuple
    title::String
    lat::Array
    lon::Array
end
function HorizonalData(post_data, name; clims = nothing)
    info = post_data.info
    clims = clims == nothing ? extrema(getproperty(info, :data)) : clims
    HorizonalData(getproperty(info, :data), clims, name, getproperty(info, :lat), getproperty(info, :lon))
end

function postprocess(name, field, p_methods; lev_slice = 5) # TODO: this should be extended to acount for different vartical spaces
    if (:regridded_3d in p_methods) || (:regridded_2d in p_methods)
        # DIR = :regridded_3d in propertynames(p_methods) ? joinpaths(REGRID_DIR,"cgll2rll_3d") : joinpaths(REGRID_DIR,"cgll2rll_2d")
        DIR = joinpath(REGRID_DIR, "cgll2rll")
        isdir(DIR) ? nothing : mkpath(DIR)
        datafile_latlon = DIR * "/remapped_" * string(name) * ".nc"
        remap_field_cgll2rll(name, field, DIR, datafile_latlon)
        grid_data = read_remapped_field(name, datafile_latlon)
    else
        @error "No regridding required. Use ClimaCorePlots for plots on the model grid."
    end


    if :zonal_mean in p_methods
        @assert (:regridded_3d in p_methods) || (:regridded_2d in p_methods)
        post_data = (;
            data = Statistics.mean(Array(FT.(getproperty(grid_data, :data))), dims = 1)[1, :, :],
            lat = getproperty(grid_data, :lat),
            lev = getproperty(grid_data, :lev),
            slice_type = :zonal_mean,
        )
    end
    if :horizontal_2d in p_methods
        @assert (:regridded_3d in p_methods) || (:regridded_2d in p_methods)
        if length(size(getproperty(grid_data, :data))) == 2
            post_data = (;
                data = getproperty(grid_data, :data),
                lat = getproperty(grid_data, :lat),
                lon = getproperty(grid_data, :lon),
                slice_type = :horizontal_2d,
            )

        else
            post_data = (;
                data = getproperty(grid_data, :data)[:, :, lev],
                lat = getproperty(grid_data, :lat),
                lon = getproperty(grid_data, :lon),
                lev = getproperty(grid_data, :lev)[lev_slice],
                slice_type = :horizontal_3d,
            )

        end
    end
    return PostProcessedData(post_data)
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
