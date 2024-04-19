# climate diagnostics: T, u, (q), rho v, N, MSE
# storm track diagnostics: [vT], [v][T], EGR = f/N dTdy

using NCDatasets
using Statistics
using Plots

# PLOT_DIR = "dry_held_suarez3/dry_held_suarez3/clima_atmos"
build = "624"
DATA_DIR = "/scratch/clima/slurm-buildkite/climacoupler-longruns/$build/climacoupler-longruns/dry_held_suarez/dry_held_suarez/clima_atmos/"

PLOT_DIR = build
mkpath(PLOT_DIR)

get_nc_data_all = (var, red, DATA_DIR) -> begin
    ds = NCDataset("$DATA_DIR/$(var)_$red.nc")
    var = ds["$var"][:,:,:,:]
    lat = ds["lat"][:]
    lon = ds["lon"][:]
    z = ds["z"][:]
    time = ds["time"][:]
    close(ds)
    return var, lat, lon, z, time
end

postprocess_climate_data = (var, red, DATA_DIR; lev_i = 1, spinup=15) -> begin

    var, lat, lon, z, time = get_nc_data_all(var, red, DATA_DIR)

    var_time_zonal_mean = mean(var[spinup:end,:,:,:], dims=(1,2))[1, 1, :, :]
    var_time_mean_sfc = mean(var[spinup:end,:,:,:], dims=(1))[1, :, :, lev_i]

    return var_time_zonal_mean, var_time_mean_sfc, lat, lon, z
end

get_timeseries_data = (variable, lon_i, lat_i, lev_i) -> begin

    variable_time_mean = mean(variable[:, lon_i[1]:lon_i[2], lat_i[1]:lat_i[2], lev_i], dims=(2,3))[:,1,1]

    return variable_time_mean
end




function plot_climate(var, DATA_DIR, PLOT_DIR; reduction = "inst")
    strf_zm, strf_sfc, lat, lon, z = postprocess_climate_data(var, reduction, DATA_DIR)
    strf_zm, strf_upper, lat, lon, z = postprocess_climate_data(var, reduction, DATA_DIR, lev_i = 39)

    # plot

    contourf(lat, z, strf_zm', xlabel="Latitude", ylabel="Height (km)", title="$var", color=:viridis, ylims = (0, 3e4), yscale=:log10, yticks = ([1e3,5e3,10e3,20e3,30e3], ["1", "5", "10", "20", "30"]))# , clims=(-1e10, 1e10))
    # yticks = ([1e3,5e3,10e3,20e3,30e3], ["1", "5", "10", "20", "30"])
    png(joinpath(PLOT_DIR, "$var.png"))

    # repeat the plot with exponential y axis
    contourf(lon, lat, strf_sfc', xlabel="Longitude", ylabel="Latitude", title="$var", color=:viridis)#, clims=(-1e10, 1e10))
    png(joinpath(PLOT_DIR, "$(var)_sfc.png"))

    contourf(lon, lat, strf_upper', xlabel="Longitude", ylabel="Latitude", title="$var", color=:viridis)#, clims=(-1e10, 1e10))
    png(joinpath(PLOT_DIR, "$(var)_10km.png"))
end

# loop plotting over variables
vars = ["mass_streamfunction", "va", "ua", "ta", "egr"]
for var in vars
    plot_climate(var, DATA_DIR, PLOT_DIR)
end

# animation of sfc lat-lon ta in time
ta_sfc, lat, lon, z, time = get_nc_data_all("ta", "inst", DATA_DIR)
ta_sfc = ta_sfc[:, :, :, 1]
anim = Plots.@animate for i in 1:size(ta_sfc, 1)
    Plots.contourf(lon, lat, ta_sfc[i, :, :, 1]', xlabel="Longitude", ylabel="Latitude", title="$var", color=:viridis, clims=(260, 315))
end
Plots.mp4(anim, joinpath(PLOT_DIR, "anim_ta_sfc.mp4"), fps = 10)

# storm tracks maps
ta_zm, ta_sfc,lat, lon, z = postprocess_climate_data("ta", "inst", DATA_DIR)
va_zm, va_sfc,lat, lon, z = postprocess_climate_data("va", "inst", DATA_DIR)
vT_zm, vT_sfc,lat, lon, z = postprocess_climate_data("vT", "inst", DATA_DIR)

heat_flux = vT_zm .- va_zm .* ta_zm
contourf(lat, z, heat_flux', xlabel="Latitude", ylabel="Height (m)", title="Heat flux", color=:viridis, ylims = (0, 1e4))# , clims=(-1e10, 1e10))
png(joinpath(PLOT_DIR, "heat_flux.png"))

vT_sfc = vT_sfc .- va_sfc .* ta_sfc
contourf(lon, lat, vT_sfc', xlabel="Longitude", ylabel="Latitude", title="Heat flux", color=:viridis)#, clims=(-1e10, 1e10))
png(joinpath(PLOT_DIR, "heat_flux)_sfc.png"))

# compute timeseries of storm track diagnostics
# egr
egr_all, lat, lon, z, time = get_nc_data_all("egr", "inst", DATA_DIR)
egr_t = get_timeseries_data(egr_all, [1, 180], [1, 90], 1)

# heat flux
vT_all, lat, lon, z, time = get_nc_data_all("vT", "inst", DATA_DIR)
va_all, lat, lon, z, time = get_nc_data_all("va", "inst", DATA_DIR)
ta_all, lat, lon, z, time = get_nc_data_all("ta", "inst", DATA_DIR)

heat_flux_all = vT_all .- va_all .* ta_all
heat_flux_t = get_timeseries_data(heat_flux_all, [1, 180], [1, 90], 1)

# plot timeseries together on different yaxis, same xaxis
nans2zero(x) = isnan(x) ? 0 : x
plot(time, egr_t, xlabel="Time", ylabel="EGR", label="EGR", color=:blue, ylims = extrema(nans2zero.(egr_t)))
# plot heatflux on the same plot with different yaxis and twinx
plt = twinx()
plot!(plt, time, heat_flux_t, ylabel="Heat flux", label="Heat flux", color=:red, ylims = extrema(nans2zero.(heat_flux_t)))

png(joinpath(PLOT_DIR, "timeseries.png"))


### interpolations to pressure


# ua, lat, lon, z, time = get_nc_data_all("ua", "inst", DATA_DIR)

# itp = LinearInterpolation((time, lon, lat, z), ua)


# pa, lat, lon, z, time = get_nc_data_all("pa", "inst", DATA_DIR)
