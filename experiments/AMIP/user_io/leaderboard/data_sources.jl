"""
    A struct to describe some observation data that we want to compare against.
"""
struct ObsDataSource

    # NOTE: This struct is not concretely typed, but we don't care about performance here.
    # We only care about beauty.

    """Path of the NetCDF file"""
    path::AbstractString

    """Name of the variable of interest in the NetCDF file"""
    var_name::AbstractString

    """Name of the time dimension in the NetCDF file"""
    time_name::AbstractString
    """Name of the longitude dimension in the NetCDF file"""
    lon_name::AbstractString
    """Name of the latitude dimension in the NetCDF file"""
    lat_name::AbstractString

    """Function that has to be applied to the data to convert it to the same conventions
    as CliMA"""
    preprocess_data_fn::Function

    """The NCDataset associated to the file"""
    ncdataset::NCDatasets.NCDataset
end

function ObsDataSource(;
    path,
    var_name,
    time_name = "time",
    lon_name = "lon",
    lat_name = "lat",
    preprocess_data_fn = identity,
)

    ncdataset = NCDatasets.NCDataset(path)

    return ObsDataSource(path, var_name, time_name, lon_name, lat_name, preprocess_data_fn, ncdataset)
end

"""
    Base.close(ds::ObsDataSource)

Close the file associated to `ds`.
"""
function Base.close(ds::ObsDataSource)
    NCDatasets.close(ds.ncdataset)
end

"""
    A struct to describe some simulation data.
"""
struct SimDataSource

    # This struct is not concretely typed, but we don't care about performance here. We only
    # care about beauty.

    """Path of the simulation output"""
    path::AbstractString

    """Short name of the variable of interest"""
    short_name::AbstractString

    """Reduction to consider"""
    reduction::AbstractString

    """Period of the reduction (e.g., 1d, 30d)"""
    period::AbstractString

    """ClimaAnalysis OutputVar"""
    var::ClimaAnalysis.OutputVar

    """Simulation longitudes and latitudes"""
    lonlat::Tuple{AbstractArray, AbstractArray}
end

function SimDataSource(; path, short_name, reduction = "average", period = "10d")

    sim = ClimaAnalysis.SimDir(path)
    # TODO: Add period, for the time-being, we just pick up what's there
    var = get(sim; short_name, reduction)

    lonlat = (var.dims["lon"], var.dims["lat"])

    return SimDataSource(path, short_name, reduction, period, var, lonlat)
end

"""
    data_at_date(sim_ds::SimDataSource, date::Dates.DateTime)

Return the simulation data at the given date.
"""
function data_at_date(sim_ds::SimDataSource, date::Dates.DateTime)
    start_date = Dates.DateTime(sim_ds.var.attributes["start_date"])
    time_diff_seconds = (date - start_date) / Dates.Second(1)
    return ClimaAnalysis.slice(sim_ds.var, time = time_diff_seconds).data
end
