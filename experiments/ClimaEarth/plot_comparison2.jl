using CairoMakie
using ClimaAnalysis

output_dir_ml = "/home/jschmitt/ClimaCoupler.jl/experiments/ClimaEarth/output/wxquest_progedmf/output_0000/clima_atmos"
output_dir_quad = "/home/jschmitt/ClimaCoupler.jl/experiments/ClimaEarth/output/wxquest_progedmf/output_0003/clima_atmos"

simdir_ml = SimDir(output_dir_ml)
simdir_quad = SimDir(output_dir_quad)


clw_dat = average_lon(get(simdir_ml; short_name = "clw", period = "1d"))

# dimensions
lat_dims = clw_dat.dims["lat"]
z_dims = clw_dat.dims["z"]
dat = clw_dat.data[1, :, :]
fig = Figure(size = (1000, 400))
ax1 = Axis(fig[1, 1], title = "ML Cloud Water", xlabel = "Latitude", ylabel = "z (m)")
CairoMakie.heatmap!(ax1, dat, colormap = :gist_ncar,
)

ax2 = Axis(fig[1, 2], title = "Quadrature Cloud Water", xlabel = "Latitude", ylabel = "z (m)")
simdir_quad = SimDir(output_dir_quad)
# average_lon(get(simdir_ml; short_name = "clw", period = "1d"))
clw_dat_quad = average_lon(get(simdir_quad; short_name = "clw", period = "1d"))
dat_quad = clw_dat_quad.data[1, :, :]
CairoMakie.heatmap!(ax2, dat_quad, colormap = :gist_ncar,
)

ax3 = Axis(fig[1, 3], title = "Difference (ML - Quad)", xlabel = "Latitude", ylabel = "z (m)")
CairoMakie.heatmap!(ax3, dat - dat_quad, colormap = Reverse(:RdBu),
)

save("clw_ml.png", fig)



simdir_ml = SimDir(output_dir_ml)
simdir_quad = SimDir(output_dir_quad)

clw_dat_ml = window(average_lon(get(simdir_ml; short_name = "clw", period = "1d")), "z", left = 0, right = 10000)
clw_dat_quad = window(average_lon(get(simdir_quad; short_name = "clw", period = "1d")), "z", left = 0, right = 10000)
# coordinates + data
lat = clw_dat_ml.dims["lat"]
z   = clw_dat_ml.dims["z_reference"]

ml   = clw_dat_ml.data[1, :, :]
quad = clw_dat_quad.data[1, :, :]
data_diff = ml .- quad

# color ranges
cr_main = extrema(vcat(ml, quad))
Δmax    = maximum(abs, data_diff)
cr_diff = (-Δmax, Δmax)

fig = Figure(size = (1300, 400))

ax1 = Axis(fig[1, 1]; title="ML Cloud Water", xlabel="Latitude", ylabel="z (m)")
ax2 = Axis(fig[1, 2]; title="Quadrature Cloud Water", xlabel="Latitude", ylabel="z (m)")
ax3 = Axis(fig[1, 4]; title="Difference (ML − Quad)", xlabel="Latitude", ylabel="z (m)")

hm1 = heatmap!(ax1, lat, z, ml;   colormap=:gist_ncar, colorrange=cr_main)
hm2 = heatmap!(ax2, lat, z, quad; colormap=:gist_ncar, colorrange=cr_main)
hm3 = heatmap!(ax3, lat, z, data_diff;
               colormap=Reverse(:RdBu), colorrange=cr_diff)

Colorbar(fig[1, 3], hm1; label="Cloud Liquid Water")
Colorbar(fig[1, 5],   hm3; label="Δ Cloud Liquid Water")

save("clw_ml.png", fig)




# estimate total cloud cover
# load cloud fraction
cl_dat = get(simdir_ml; short_name="cl", period="1d")

# coordinates
z = cl_dat.dims["z_reference"]          # layer centers [m]
cl = cl_dat.data[1, :, :, :]

# compute dz on staggered grid
dz = diff(z[:])
dz = vcat(dz[1], dz)                    # match size(z)

# vertical integral → cloud depth [m]
dz3 = reshape(dz, 1, 1, :)
cloud_depth = sum(cl .* dz3 .* .01, dims = 3)
cloud_depth = dropdims(cloud_depth; dims = 3)

fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1]; title="Estimated Total Cloud Depth", xlabel="Longitude", ylabel="Latitude")
hm = heatmap!(ax, cl_dat.dims["lon"], cl_dat.dims["lat"], cloud_depth; colormap=:gist_ncar)
Colorbar(fig[1, 2], hm; label="Cloud Depth (m)")

save("cloud_depth_ml.png", fig)




using CairoMakie

# --- helper: cloud depth from volumetric cloud fraction ---
cloud_depth_from_cl(cl_dat) = begin
    z  = cl_dat.dims["z_reference"]
    cl = cl_dat.data[1, :, :, :]              # (lon, lat, z)

    dz = diff(z[:])
    dz = vcat(dz[1], dz)                      # size(z)
    dz3 = reshape(dz, 1, 1, :)

    dropdims(sum(cl .* dz3 .* 0.01, dims=3); dims=3)
end

# --- load data ---
cl_ml   = get(simdir_ml;   short_name="cl", period="1d")
cl_quad = get(simdir_quad; short_name="cl", period="1d")

cd_ml   = cloud_depth_from_cl(cl_ml)
cd_quad = cloud_depth_from_cl(cl_quad)
cd_diff = cd_ml .- cd_quad

lon = cl_ml.dims["lon"]
lat = cl_ml.dims["lat"]

# color ranges
cr_main = extrema(vcat(cd_ml, cd_quad))
Δmax    = 1000# maximum(abs, cd_diff)
cr_diff = (-Δmax, Δmax)

# --- plotting ---
fig = Figure(size=(1300, 400))

ax1 = Axis(fig[1,1]; title="ML Cloud Depth",   xlabel="Longitude", ylabel="Latitude")
ax2 = Axis(fig[1,2]; title="Quad Cloud Depth", xlabel="Longitude", ylabel="Latitude")
ax3 = Axis(fig[1,4]; title="Difference (ML − Quad)", xlabel="Longitude", ylabel="Latitude")

hm1 = heatmap!(ax1, lon, lat, cd_ml;   colormap=:gist_ncar, colorrange=cr_main)
hm2 = heatmap!(ax2, lon, lat, cd_quad; colormap=:gist_ncar, colorrange=cr_main)
hm3 = heatmap!(ax3, lon, lat, cd_diff;
               colormap=Reverse(:RdBu), colorrange=cr_diff)

Colorbar(fig[1,3], hm1; label="Cloud Depth (m)")
Colorbar(fig[1,5],   hm3; label="Δ Cloud Depth (m)")

save("cloud_depth_ml_vs_quad.png", fig)
