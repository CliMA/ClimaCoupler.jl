# paper_figs
using NCDatasets
using Statistics
import Makie
import CairoMakie

import Interpolations: LinearInterpolation
import DelimitedFiles: writedlm, readdlm

include("plot_helper.jl")

for job_id in ["dry_held_suarez", "moist_held_suarez"]
    if isinteractive()
        DATA_DIR = "experiments/ClimaEarth/$job_id/$job_id/clima_atmos/output_active/"
    else
        build = ENV["BUILDKITE_BUILD_NUMBER"]
        DATA_DIR = "/central/scratch/esm/slurm-buildkite/climacoupler-ci/$build/climacoupler-ci/$job_id/output_active/clima_atmos/"
    end

    reduction = "6h_inst"
    PLOT_DIR = "paper_figs"

    mkpath(PLOT_DIR)

    # SUPPLEMENTAL: animation of surface temperature
    ta_sfc, lat, lon, z, time = get_nc_data_all("ta", reduction, DATA_DIR)
    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1], xlabel = "Longitude", ylabel = "Latitude", title = "ta_sfc")
    # plot once before animation to set axis
    co_ta_sfc = Makie.contourf!(
        ax,
        lon,
        lat,
        ta_sfc[1, :, :, 1],
        colormap = :viridis,
        levels = [k for k in range(260.0, 315.0)],
    )
    Makie.Colorbar(f[1, 2], co_ta_sfc)

    Makie.record(f, joinpath(PLOT_DIR, "anim_ta_sfc.mp4"), 1:size(ta_sfc, 1); framerate = 10) do t
        Makie.contourf!(lon, lat, ta_sfc[t, :, :, 1], colormap = :viridis, levels = [k for k in range(260.0, 315.0)])
    end
    # Figure 2: climatology
    # this plots the time-mean (upper-level and surface) slices and zonal means of
    # the mass streamfunction, zonal wind, meridional wind, temperature, max. Eady growth rate, and vertical velocity
    vars = ["mass_strf", "va", "ua", "ta", "egr", "wa"]
    for var in vars
        plot_climate(var, DATA_DIR, PLOT_DIR, job_id, reduction = reduction, interpolate_to_pressure = true)
    end

    # Figure 4: storm track diagnostics
    # this extracts the eddy heat flux and plots its climatology
    lev_st = 6
    ta_zm, ta_sfc, lat, lon, z = mean_climate_data("ta", reduction, DATA_DIR, lev_i = lev_st)
    va_zm, va_sfc, lat, lon, z = mean_climate_data("va", reduction, DATA_DIR, lev_i = lev_st)
    vT_zm, vT_sfc, lat, lon, z = mean_climate_data("vt", reduction, DATA_DIR, lev_i = lev_st)
    heat_flux_zm = vT_zm .- va_zm .* ta_zm

    pa_zm, ~, ~, ~, ~ = mean_climate_data("pfull", reduction, DATA_DIR)
    pa_zm = pa_zm ./ 100 # convert to hPa
    pa_grid = [950, 800, 700, 600, 500, 400, 300, 200, 50]

    heat_flux_int_zm = interpolate_to_pressure_coord_2d(heat_flux_zm, pa_zm, pa_grid)
    co_heat_flux = Makie.contourf(lat, -pa_grid, heat_flux_int_zm, colormap = :viridis, levels = 16)
    co_heat_flux.axis.xlabel = "Latitude (deg N)"
    co_heat_flux.axis.ylabel = "Pressure (hPa)"
    co_heat_flux.axis.title = "Heat flux"
    Makie.ylims!(co_heat_flux.axis, -pa_grid[1], -pa_grid[end])
    co_heat_flux.axis.yticks = (-pa_grid, string.(pa_grid))
    Makie.save(joinpath(PLOT_DIR, "$(job_id)_heat_flux.png"), co_heat_flux)
end
