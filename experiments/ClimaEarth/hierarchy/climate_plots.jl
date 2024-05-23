# paper_figs
using NCDatasets
using Statistics
using Plots

import Interpolations: LinearInterpolation
import DelimitedFiles: writedlm, readdlm

include("plot_helper.jl")

for job_id in ["dry_held_suarez", "moist_held_suarez"]
    if isinteractive()
        DATA_DIR = "experiments/ClimaEarth/$job_id/$job_id/clima_atmos/output_active/"
    else
        build = ENV["BUILDKITE_BUILD_NUMBER"]
        DATA_DIR = "/scratch/clima/slurm-buildkite/climacoupler-hierarchies/$build/climacoupler-hierarchies/$job_id/$job_id/clima_atmos/output_active/"
    end

    reduction = "6h_inst"
    PLOT_DIR = "paper_figs"

    mkpath(PLOT_DIR)

    # SUPPLEMENTAL: animation of surface temperature
    ta_sfc, lat, lon, z, time = get_nc_data_all("ta", reduction, DATA_DIR)
    anim = Plots.@animate for i in 1:size(ta_sfc, 1)
        Plots.contourf(
            lon,
            lat,
            ta_sfc[i, :, :, 1]',
            xlabel = "Longitude",
            ylabel = "Latitude",
            title = "$var",
            color = :viridis,
            clims = (260, 315),
        )
    end
    Plots.mp4(anim, joinpath(PLOT_DIR, "anim_ta_sfc.mp4"), fps = 10)


    upper_level = 10

    # F2: climatology
    vars = ["mass_strf", "va", "ua", "ta", "egr", "wa"]
    for var in vars
        plot_climate(var, DATA_DIR, PLOT_DIR, job_id, reduction = reduction, interpolate_to_pressure = true)
    end

    # F4: storm track diagnostics: time mean maps
    lev_st = 6
    ta_zm, ta_sfc, lat, lon, z = mean_climate_data("ta", reduction, DATA_DIR, lev_i = lev_st)
    va_zm, va_sfc, lat, lon, z = mean_climate_data("va", reduction, DATA_DIR, lev_i = lev_st)
    vT_zm, vT_sfc, lat, lon, z = mean_climate_data("vt", reduction, DATA_DIR, lev_i = lev_st)
    heat_flux_zm = vT_zm .- va_zm .* ta_zm

    pa_zm, ~, ~, ~, ~ = mean_climate_data("pfull", reduction, DATA_DIR)
    pa_zm = pa_zm ./ 100 # convert to hPa
    pa_grid = [950, 800, 700, 600, 500, 400, 300, 200, 50]

    heat_flux_int_zm = interpolate_to_pressure_coord_2d(heat_flux_zm, pa_zm, pa_grid)
    Plots.contourf(
        lat,
        -pa_grid,
        heat_flux_int_zm',
        xlabel = "Latitude (deg N)",
        ylabel = "Pressure (hPa)",
        title = "Heat flux",
        color = :viridis,
        ylims = (-pa_grid[1], -pa_grid[end]),
        yticks = (-pa_grid, pa_grid),
    )# , clims=(-1e10, 1e10))
    png(joinpath(PLOT_DIR, "$(job_id)_heat_flux.png"))

    # F5: storm track diagnostics: timeseries
    lev_i, lat_s_i, lat_n_i, lon_w_i, lon_e_i = lev_st, 60, 75, 1, 30
    println(
        "Sectorial selevtion for timeseries: \n level: $(z[lev_i]), lat: $(lat[lat_s_i]) to $(lat[lat_n_i]), lon: $(lon[lon_w_i]) to $(lon[lon_e_i])",
    )

    egr_all, lat, lon, z, time = get_nc_data_all("egr", reduction, DATA_DIR)
    egr_t = point_timeseries_data(egr_all, [lon_w_i, lon_e_i], [lat_s_i, lat_n_i], lev_i)

    vT_all, lat, lon, z, time = get_nc_data_all("vt", reduction, DATA_DIR)
    va_all, lat, lon, z, time = get_nc_data_all("va", reduction, DATA_DIR)
    ta_all, lat, lon, z, time = get_nc_data_all("ta", reduction, DATA_DIR)
    heat_flux_all = vT_all .- va_all .* ta_all
    heat_flux_t = point_timeseries_data(heat_flux_all, [lon_w_i, lon_e_i], [lat_s_i, lat_n_i], lev_i)
end
