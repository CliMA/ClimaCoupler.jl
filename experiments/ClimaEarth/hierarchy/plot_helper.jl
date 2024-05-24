# plot_helper

"""
    get_nc_data_all(var, red, DATA_DIR)

Reads the netcdf file for the variable `var` and reduction `red` from the directory `DATA_DIR` and returns the variable, lat, lon, z, and time.
"""
get_nc_data_all = (var, red, DATA_DIR) -> begin
    ds = NCDataset("$DATA_DIR/$(var)_$red.nc")
    var = ds["$var"][:, :, :, :]
    lat = ds["lat"][:]
    lon = ds["lon"][:]
    z = ds["z"][:]
    time = ds["time"][:]
    close(ds)
    return var, lat, lon, z, time
end

"""
    mean_climate_data(var, red, DATA_DIR; lev_i = 1, spinup=1)

Postprocesses the climate data for the variable `var` and reduction `red` from the directory `DATA_DIR`. Returns the zonal mean and horizontal surface slice mean of the variable.
"""
mean_climate_data =
    (varname, reduction, DATA_DIR; lev_i = 2, spinup = 1) -> begin

        var, lat, lon, z, time = get_nc_data_all(varname, reduction, DATA_DIR)
        @assert spinup < size(var, 1)

        var_time_zonal_mean = mean(var[spinup:end, :, :, :], dims = (1, 2))[1, 1, :, :]
        var_time_mean_sfc = mean(var[spinup:end, :, :, :], dims = (1))[1, :, :, lev_i]

        return var_time_zonal_mean, var_time_mean_sfc, lat, lon, z
    end

"""
    point_timeseries_data(variable, lon_i, lat_i, lev_i)

Returns the time series data for the variable `variable` at the indices `lon_i`, `lat_i`, and `lev_i`.
"""
point_timeseries_data =
    (variable, lon_i, lat_i, lev_i) -> begin

        variable_time_mean = mean(variable[:, lon_i[1]:lon_i[2], lat_i[1]:lat_i[2], lev_i], dims = (2, 3))[:, 1, 1]

        return variable_time_mean
    end

"""
    plot_climate(var, DATA_DIR, PLOT_DIR; reduction = "inst")

Plots the zonal mean and horizontal surface slice mean of the variable `var` from the directory `DATA_DIR` and saves the plots in the directory `PLOT_DIR`.
"""
function plot_climate(var, DATA_DIR, PLOT_DIR, job_id; reduction = "inst", interpolate_to_pressure = false)
    strf_zm, strf_sfc, lat, lon, z = mean_climate_data(var, reduction, DATA_DIR)
    strf_zm, strf_upper, lat, lon, z = mean_climate_data(var, reduction, DATA_DIR, lev_i = 10)

    # vertical-lat plot of zonal and time mean
    if interpolate_to_pressure
        pa_zm, ~, ~, ~, ~ = mean_climate_data("pfull", reduction, DATA_DIR)
        pa_zm = pa_zm ./ 100 # convert to hPa
        pa_grid = [950, 800, 700, 600, 500, 400, 300, 200, 50]
        strf_zm = interpolate_to_pressure_coord_2d(strf_zm, pa_zm, pa_grid)
        Plots.contourf(
            lat,
            -pa_grid,
            strf_zm',
            xlabel = "Latitude (deg N)",
            ylabel = "Pressure (hPa)",
            title = "$var",
            color = :viridis,
            ylims = (-pa_grid[1], -pa_grid[end]),
            yticks = (-pa_grid, pa_grid),
        )# , clims=(-1e10, 1e10))
        png(joinpath(PLOT_DIR, "$(job_id)_$(var)_pa.png"))
    else
        Plots.contourf(
            lat,
            z,
            strf_zm',
            xlabel = "Latitude",
            ylabel = "Height (km)",
            title = "$var",
            color = :viridis,
            ylims = (0, 3e4),
            yscale = :log10,
            yticks = ([1e3, 5e3, 10e3, 20e3, 30e3], ["1", "5", "10", "20", "30"]),
        )# , clims=(-1e10, 1e10)
        png(joinpath(PLOT_DIR, "$(job_id)_$var.png"))
    end

    # horizontal slices
    Plots.contourf(lon, lat, strf_sfc', xlabel = "Longitude", ylabel = "Latitude", title = "$var", color = :viridis)#, clims=(-1e10, 1e10))
    png(joinpath(PLOT_DIR, "$(job_id)_$(var)_sfc.png"))

    Plots.contourf(lon, lat, strf_upper', xlabel = "Longitude", ylabel = "Latitude", title = "$var", color = :viridis)#, clims=(-1e10, 1e10))
    png(joinpath(PLOT_DIR, "$(job_id)_$(var)_10km.png"))
end

"""
    interpolate_to_pressure_coord_2d(var_zm, pa, pa_grid)

Interpolates the 2D variable `var_zm` to the pressure grid `pa_grid` using the pressure values `pa`.
"""
function interpolate_to_pressure_coord_2d(var_zm, pa, pa_grid)
    var_on_pa = zeros(size(var_zm, 1), length(pa_grid))
    for lat_i in collect(1:size(var_zm, 1))
        # Extract ua and corresponding ta values
        var_values = var_zm[lat_i, :]
        pa_values = pa[lat_i, :]

        # Interpolate ua onto ta_grid
        for (pa_j, pa_val) in enumerate(pa_grid)
            itp_var = LinearInterpolation(-pa_values, var_values)
            var_on_pa[lat_i, pa_j] = itp_var(-pa_val)
        end
    end
    return var_on_pa
end
