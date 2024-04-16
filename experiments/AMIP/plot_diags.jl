# climate plots
include("user_io/ci_plots.jl")
function make_plots(
    ::Union{Val{:general_plots}},
    output_paths::Vector{<:AbstractString},
    plot_path::AbstractString;
    reduction::String = "average",
)
    simdirs = CAN.SimDir.(output_paths)

    # Default output diagnostics
    short_names_3D = ["mse", "lr"]
    short_names_2D = ["ts"]

    available_periods = CAN.available_periods(simdirs[1]; short_name = short_names_3D[1], reduction)
    period = ""
    if "10d" in available_periods
        period = "10d"
    elseif "1d" in available_periods
        period = "1d"
    elseif "12h" in available_periods
        period = "12h"
    end

    # Creates diagnostics vector
    # 3D fields are zonally averaged platted onf the lat-z plane
    # 2D fields are plotted on the lon-lat plane
    vars_3D = map_comparison(simdirs, short_names_3D) do simdir, short_name
        get(simdir; short_name, reduction, period) |> CAN.average_lon
    end

    available_periods = CAN.available_periods(simdirs[1]; short_name = short_names_2D[1], reduction)

    vars_2D = map_comparison(simdirs, short_names_2D) do simdir, short_name
        get(simdir; short_name, reduction, period)
    end

    make_plots_generic(output_paths, plot_path, vars_3D, time = LAST_SNAP, more_kwargs = YLINEARSCALE)
    make_plots_generic(output_paths, plot_path, vars_2D, time = LAST_SNAP, output_name = "summary_2D")

end


COUPLER_OUTPUT_DIR = "/scratch/clima/slurm-buildkite/climacoupler-longruns/621/climacoupler-longruns/dry_held_suarez/dry_held_suarez/clima_atmos/"
make_plots(Val(:general_plots), [COUPLER_OUTPUT_DIR], COUPLER_OUTPUT_DIR)


# climate diagnostics: T, u, (q), rho v, N, MSE

# storm track diagnostics: [vT], [v][T], EGR = f/N dTdy


# explote data on sampo
using NCDatasets
using Statistics

DIAG_DIR = "data_623/"
var, red = ("mass_streamfunction", "inst")
# var, red = ("va", "30d_average")
var, red = ("ua", "30d_average")
# var, red = ("ta", "30d_average")
ds = NCDataset("$DIAG_DIR/$(var)_$red.nc")
lat = ds["lat"][:]
lon = ds["lon"][:]
z = ds["z"][:]
time = ds["time"][:]
strf = ds["$var"][:,:,:,:] #.* 2π * 6371e3 .* cosd.(reshape(lat, 1,1,length(lat),1)) # kg s^-1
close(ds)

strf_time_zonal_mean = mean(strf, dims=(1,2))[1, 1, :, :]
strf_time_mean_sfc = mean(strf, dims=(1))[1, :, :, 39]

# plot
using Plots
contourf(lat, z, strf_time_zonal_mean', xlabel="Latitude", ylabel="Height (m)", title="$var", color=:viridis, ylims = (0, 1e4))# , clims=(-1e10, 1e10))
png(joinpath(DIAG_DIR, "$var.png"))

contourf(lon, lat, strf_time_mean_sfc', xlabel="Longitude", ylabel="Latitude", title="$var", color=:viridis)#, clims=(-1e10, 1e10))
png(joinpath(DIAG_DIR, "$(var)_10k.png"))

using ClimaCorePlots


plot(CC.Geometry.UVVector.(atmos_sim.integrator.u.c.uₕ).components.data.:1)
png("cc_u")

sol_atm = atmos_sim.integrator.sol;

anim = Plots.@animate for u in sol_atm.u
    Plots.plot(CC.Fields.level(CC.Geometry.UVVector.(u.c.uₕ).components.data.:1, 5))
end
Plots.mp4(anim, joinpath(DIAG_DIR, "anim_u.mp4"), fps = 10)