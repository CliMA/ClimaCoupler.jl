import OrderedCollections: OrderedDict

"""
    isequispaced(arr::Vector)

Return whether the array is equispaced or not.
"""
function isequispaced(arr::Vector)
    return all(diff(arr) .≈ arr[begin + 1] - arr[begin])
end

"""
    resample(
        data::AbstractArray,
        src_lonlat::Tuple{<:AbstractArray, <:AbstractArray},
        dest_lonlat::Tuple{<:AbstractArray, <:AbstractArray},
        )

Resample 2D `data` from `src_lonlat` to `dest_lonlat`.

Note, "flat" boundary conditions are imposed for the outer edges. This should make sense for
most data sources (that are defined all over the globe but do not necessarily have points at
the poles). Be sure to check that it makes sense for your data.

The resampling performed here is a 0th-order constant resampling.
"""

function resample(
    data::AbstractArray,
    src_lonlat::Tuple{<:AbstractArray, <:AbstractArray},
    dest_lonlat::Tuple{<:AbstractArray, <:AbstractArray},
)

    # Interpolations.jl wants ranges, so we have to make ranges out of the given src_lonlat
    #
    # NOTE: We are assuming that lonlat are equispaced. Check this!
    vec_to_range(v) = range(v[begin], v[end], length = length(v))
    src_lonlat_ranges = vec_to_range.(src_lonlat)

    itp = Interpolations.constant_interpolation(
        src_lonlat_ranges,
        data,
        extrapolation_bc = (Interpolations.Periodic(), Interpolations.Flat()),
    )

    dest_lon, dest_lat = dest_lonlat

    return [itp(lon, lat) for lon in dest_lon, lat in dest_lat]
end

"""
    integration_weights(lonlat::Tuple{<:AbstractArray, <:AbstractArray})

Compute the integration weights for a first-order spherical integration.
"""
function integration_weights(lonlat::Tuple{<:AbstractArray, <:AbstractArray})
    lon, lat = lonlat

    # We are also assuming that they are in degrees
    abs(maximum(lon)) <= π && error("longitude should be in degrees")
    abs(maximum(lat)) <= π / 2 && error("latitude should be in degrees")

    isequispaced(lon) || error("Longitude is not equispaced")
    isequispaced(lat) || error("Latitude is not equispaced")

    dlon_rad = deg2rad(abs(lon[begin + 1] - lon[begin]))
    dlat_rad = deg2rad(abs(lat[begin + 1] - lat[begin]))

    return [cosd(lat1) * dlon_rad * dlat_rad for _ in lon, lat1 in lat]
end

"""
    integrate_on_sphere(
        data::AbstractArray,
        lonlat::Tuple{<:AbstractArray, <:AbstractArray},
        )

Integrate `data` on a sphere with a first-order scheme. `data` has to be discretized on
`lonlat`.
"""
function integrate_on_sphere(data::AbstractArray, lonlat::Tuple{<:AbstractArray, <:AbstractArray})
    lon, lat = lonlat
    size_data = size(data)
    exp_size = (length(lon), length(lat))
    size_data == exp_size || error("Inconsistent dimensions $size_data != $exp_size")
    return sum(data .* integration_weights(lonlat)) ./ sum(integration_weights(lonlat))
end

"""
    mse(
        data1::AbstractArray,
        data2::AbstractArray,
        lonlat::Tuple{<:AbstractArray, <:AbstractArray},
        )

Compute the 2D map of the mean-square error.
"""
function mse(data1::AbstractArray, data2::AbstractArray, lonlat::Tuple{<:AbstractArray, <:AbstractArray})
    lon, lat = lonlat
    size_data1 = size(data1)
    size_data2 = size(data2)
    exp_size = (length(lon), length(lat))
    size_data1 == exp_size || error("Inconsistent dimensions $size_data1 != $exp_size")
    size_data1 == size_data2 || error("Inconsistent dimensions $size_data1 != $size_data2")

    return (data1 .- data2) .^ 2
end

"""
    bias(
        sim::AbstractArray,
        obs::AbstractArray,
        lonlat::Tuple{<:AbstractArray, <:AbstractArray},
        )

Compute the 2D map of the bias (simulated - observed).
"""
function bias(sim::AbstractArray, obs::AbstractArray, lonlat::Tuple{<:AbstractArray, <:AbstractArray})
    lon, lat = lonlat
    size_sim = size(sim)
    size_obs = size(obs)
    exp_size = (length(lon), length(lat))
    size_sim == exp_size || error("Inconsistent dimensions $size_sim != $exp_size")
    size_sim == size_obs || error("Inconsistent dimensions $size_obs != $size_data2")

    return sim .- obs
end

"""
    bias(
        obs_ds::ObsDataSource,
        sim_ds::SimDataSource,
        target_dates::AbstractArray{<: Dates.DateTime},
        )

Compute the 2D map of the bias (simulated - observed) for data averaged over the given dates.

The return value is a `ClimaAnalysis.OutputVar` with the computed `rmse` and `bias` among
its attributes.
"""
function bias(obs_ds::ObsDataSource, sim_ds::SimDataSource, target_dates::AbstractArray{<:Dates.DateTime})
    lonlat = sim_ds.lonlat
    simulated_data = map(d -> data_at_date(sim_ds, d), target_dates) |> mean
    observational_data = map(d -> find_and_resample(obs_ds, d, lonlat), target_dates) |> mean

    bias_arr = bias(simulated_data, observational_data, lonlat)
    mse_arr = mse(simulated_data, observational_data, lonlat)

    short_name = ClimaAnalysis.short_name(sim_ds.var)

    bias_dims = OrderedDict("lon" => lonlat[1], "lat" => lonlat[2])
    bias_dim_attribs = Dict{String, Any}()

    rmse = round(sqrt(integrate_on_sphere(mse_arr, lonlat)); sigdigits = 3)
    global_bias = round(integrate_on_sphere(bias_arr, lonlat); sigdigits = 3)

    units = isnothing(sim_ds.new_units) ? sim_ds.var.attributes["units"] : sim_ds.new_units

    bias_attribs = Dict{String, Any}(
        "short_name" => "sim-obs_$short_name",
        "var_short_name" => "$short_name",
        "long_name" => "SIM - OBS mean $short_name\n(RMSE: $rmse $units, Global bias: $global_bias $units)",
        "rmse" => rmse,
        "bias" => global_bias,
        "units" => units,
    )

    return ClimaAnalysis.OutputVar(bias_attribs, bias_dims, bias_dim_attribs, bias_arr)
end

"""
    integrate_on_sphere(var::ClimaAnalysis.OutputVar)

Integrate the given `var` onto the sphere with a first order integration scheme.
"""
function integrate_on_sphere(var::ClimaAnalysis.OutputVar)
    lonlat = (var.dims["lon"], var.dims["lat"])
    return integrate_on_sphere(var.data, lonlat)
end

"""
    find_and_resample(
        observed_data::ObsDataSource,
        date::Dates.DateTime,
        dest_lonlat::Tuple{<:AbstractArray, <:AbstractArray},
        )

Find the data corresponding to the given `date` in `observed_data` and resample it to
`dest_lonlat`.
"""
function find_and_resample(
    observed_data::ObsDataSource,
    date::Dates.DateTime,
    dest_lonlat::Tuple{<:AbstractArray, <:AbstractArray},
)
    obs = observed_data

    available_times = obs.ncdataset[observed_data.time_name]

    if observed_data.shift_to_end_of_month
        available_times = Dates.DateTime.(Dates.lastdayofmonth.(available_times))
    end

    time_index = ClimaAnalysis.Utils.nearest_index(available_times, date)

    # NOTE: We are hardcoding that the time index is the last one and that there are three
    # dimensions! We should generalize this (the NetCDF file contains all the information to
    # deduce this).

    data_arr = obs.preprocess_data_fn(obs.ncdataset[obs.var_name][:, :, time_index])

    lon_arr = obs.ncdataset[obs.lon_name][:]
    lat_arr = obs.ncdataset[obs.lat_name][:]

    return resample(data_arr, (lon_arr, lat_arr), dest_lonlat)
end

"""
    split_by_season(dates::AbstractArray{<: Dates.DateTime})

Take an array of dates and return 4 split into seasons.
"""
function split_by_season(dates::AbstractArray{<:Dates.DateTime})
    MAM, JJA, SON, DJF =
        Vector{Dates.DateTime}(), Vector{Dates.DateTime}(), Vector{Dates.DateTime}(), Vector{Dates.DateTime}()

    for date in dates
        if Dates.Month(3) <= Dates.Month(date) <= Dates.Month(5)
            push!(MAM, date)
        elseif Dates.Month(6) <= Dates.Month(date) <= Dates.Month(8)
            push!(JJA, date)
        elseif Dates.Month(9) <= Dates.Month(date) <= Dates.Month(11)
            push!(SON, date)
        else
            push!(DJF, date)
        end
    end

    return (MAM, JJA, SON, DJF)
end
