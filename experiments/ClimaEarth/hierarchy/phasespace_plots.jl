using Dates
using DelimitedFiles
using PyPlot
using LinearAlgebra

# Parameters
σ_x = 0.08 # size of filter
σ_y = 0.08 # size of filter
grid_scale = 0.05 # grid scale; 1/grid_scale = number of grid points

dd = collect(range(0, 9, step=1))
dd_2d = repeat(reshape(dd, 1, length(dd)), length(dd), 1)

# Load data
vt_all = readdlm("624/vtLongTsrsPAC.txt")[:, :] # yr, t
egr_all = readdlm("624/bcLongTsrsPAC.txt")[:, :] # yr, t

# flatten
vt = vt_all'[:]
egr = egr_all'[:]

# Calculate tendencies

"""
    get_tendencies(vt)
Use centered differences to calculate the tendencies of a variable.
"""
function get_tendencies(v)
    Δv = zeros(size(v))
    for i in 2:size(v)[1]-1
        Δv[i] = (v[i+1] - v[i-1]) / 2
    end
    return Δv[2:end-1]
end

ddt_vt = get_tendencies(vt)
ddt_egr = get_tendencies(egr)

# Remove the first and last values
vt = vt[2:end-1]
egr = egr[2:end-1]

# Define some parameters
"""
    get_grid_axis(v)
Get the grid axis for a variable.
"""
function get_grid_axis(v, grid_scale = 0.1)
    v_min, v_max = extrema(v)
    Δv = (v_max - v_min) * grid_scale
    return range(v_min, v_max, step=Δv)
end

grid_vt = get_grid_axis(vt, grid_scale)
grid_egr = get_grid_axis(egr, grid_scale)

# Normalize variables
vt_n = vt ./ maximum(vt)
egr_n = egr ./ maximum(egr)
grid_vt_n = grid_vt ./ maximum(grid_vt)
grid_egr_n = grid_egr ./ maximum(grid_egr)

# Calculate smoothing weights
"""
    get_weights(vx, vy, gx, gy; σ_x = 0.05, σ_y = 0.05)
Calculate the smoothing weights for a variable.
# Arguments
- `vx::Array`: variable x
- `vy::Array`: variable y
- `gx::Array`: grid x
- `gy::Array`: grid y
- `σ_x::Float64`: standard deviation x
- `σ_y::Float64`: standard deviation y
"""
function get_weights(vx, vy, gx, gy; σ_x = 0.05, σ_y = 0.05)
    N = 1 / (σ_x * σ_y * 2π)
    TwoSigSqX = 2σ_x^2
    TwoSigSqY = 2σ_y^2
    weights = zeros(length(gx), length(gy), length(vx))
    for xx in 1:length(gx)
        for yy in 1:length(gy)
            for apt in 1:length(vx)
                wx = (gx[xx] - vx[apt])^2 / TwoSigSqX
                wy = (gy[yy] - vy[apt])^2 / TwoSigSqY
                weights[xx, yy, apt] = N * exp(-wx - wy)
            end
        end
    end
    return weights
end

weights = get_weights(vt_n, egr_n, grid_vt_n, grid_egr_n, σ_x = σ_x, σ_y = σ_y)

function kernel_smooth(v, gx, gy, weights)
    v_smooth = zeros(length(gx), length(gy))
    for xx in 1:length(gx)
        for yy in 1:length(gy)
            weights_n = weights[xx, yy, :] / sum(weights[xx, yy, :])
            weights_n[isnan.(weights_n)] .= 0.0
            v_smooth[xx, yy] = sum(v .* weights_n)
        end
    end
    return v_smooth
end

# Calculate smoothed tendencies
ddt_vt_smooth = kernel_smooth(ddt_vt, grid_vt_n, grid_egr_n, weights)
ddt_egr_smooth = kernel_smooth(ddt_egr, grid_vt_n, grid_egr_n, weights)

# make all grid points with 2 or fewer values NaN
function make_zero(v_x, v_y ,dv, grid_vt, grid_egr; threshold = 1)
    Δvt = grid_vt.step.hi
    Δegr = grid_egr.step.hi
    count = zeros(length(grid_vt), length(grid_egr))
    for xx in 1:length(grid_vt)
        for yy in 1:length(grid_egr)
            for apt in 1:length(v_x)
                if (grid_vt[xx] - Δvt / 2.0 < v_x[apt]) && (grid_vt[xx] + Δvt / 2.0 >= v_x[apt]) &&
                    (grid_egr[yy] - Δegr / 2.0 < v_y[apt]) && (grid_egr[yy] + Δegr / 2.0 >= v_y[apt])
                    count[xx, yy] += 1
                end
            end
        end
    end

    for xx in 1:length(grid_vt)
        for yy in 1:length(grid_egr)
            if count[xx, yy] ≤ threshold
                dv[xx, yy] = 0
            end
        end
    end
    return dv
end

ddt_vt_smooth = make_zero(vt, egr , ddt_vt_smooth, grid_vt, grid_egr)
ddt_egr_smooth = make_zero(vt, egr, ddt_egr_smooth, grid_vt, grid_egr)

# Plotting vector plot
grid_2d_vt, grid_2d_egr = meshgrid(grid_vt, grid_egr)

using GR

Plots.plot()
Plots.quiver(collect(grid_2d_vt), collect(grid_2d_egr), quiver = (ddt_vt_smooth, ddt_egr_smooth)) # )color="k", )# linewidth=5.0 .* sqrt(ddt_vt_smooth.^2 .+ ddt_egr_smooth.^2) ./ maximum(sqrt(ddt_vt_smooth.^2 .+ ddt_egr_smooth.^2)))
xlabel("H Fl")
ylabel("BC")
png("624/vector_plot1.png")

# Plotting streamlines
using Makie
using CairoMakie
function f_stream(x, y)
    index_x = argmin(abs.(x .- grid_vt))
    index_y = argmin(abs.(y .- grid_egr))
    return Point2f(ddt_vt_smooth[index_x[1], index_y[1]], ddt_egr_smooth[index_x[1], index_y[1]])
end

# color function is the log of a norm

# function f_color(p::Point)
#     x, y = p[1], p[2]
#     index_x = argmin(abs.(x .- grid_vt))
#     index_y = argmin(abs.(y .- grid_egr))
#     return exp(norm([ddt_vt_smooth[index_x[1], index_y[1]], ddt_egr_smooth[index_x[1], index_y[1]]]))
# end

using StatsBase


fig = Makie.Figure()
ax = Makie.Axis(fig[1, 1], xlabel = "heat flux [K m s^-1]", ylabel = "baroclinicity [day^-1]", title = "phase speed plot")

# add a gray elliptic oval with radii of the standard deviations
function ellipse(x, y, σ_x, σ_y)
    t = range(0, 2π, length=100)
    return x .+ σ_x .* cos.(t), y .+ σ_y .* sin.(t)
end
# plot the ellipse
x, y = ellipse(0.0, 0.0, σ_x * maximum(vt), σ_y * maximum(egr))
Makie.lines!(ax, x, y, color = :gray, linewidth = 2.0)

# add 2d histogram of the data with transparency of 0.5
hist = fit(Histogram, (vt, egr), (grid_vt, grid_egr))
Makie.surface!(ax, hist, alpha = 0.5, colormap = Reverse(:grays))



# add streamlines
# Makie.streamplot!(ax, f_stream, grid_2d_vt, grid_2d_egr) #, color = f_color)
# streamplot with a transparent background
streamplot!(ax, f_stream, grid_2d_vt, grid_2d_egr, linewidth = 2.0)

Makie.save("624/streamlines1.png", fig, px_per_unit = 2)



# fig, ax, pl = Makie.scatter(vt, egr, color = :black)
# Makie.save("624/scatter.png", fig, px_per_unit = 2)





# using PyPlot
# using PyCall
# @pyimport matplotlib.patches as patch

# cfig = figure()
# ax = cfig[:add_subplot](1,1,1)
# ax[:set_aspect]("equal")
# c = patch.Circle([0.5,0.5],0.4,fc="blue",ec="red",linewidth=.5,zorder=0)
# ax[:add_artist](c)
# cfig[:savefig]("circle.png")