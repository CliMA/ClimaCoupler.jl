#=
# Flux snapshot plots

Maps of one coupler flux on every grid it passes through, read from the files
`SimOutput.write_flux_snapshot` writes. Each grid is drawn natively — exchange
polygons as polygons, SE nodes as points, FV cells as quads — so a single bad node
or polygon stays visible instead of being smoothed away by remapping.
=#

import JLD2
import Makie
import Printf
import ClimaCoupler: Plotting

const FLUX_SNAPSHOT_UNITS = Dict(
    "F_sh" => "W m⁻²",
    "F_lh" => "W m⁻²",
    "F_moisture" => "kg m⁻² s⁻¹",
    "F_τu" => "N m⁻²",
    "F_τv" => "N m⁻²",
    "SW_d" => "W m⁻²",
    "LW_d" => "W m⁻²",
    "P_liq" => "kg m⁻² s⁻¹",
    "P_snow" => "kg m⁻² s⁻¹",
)

# Nodes with weighted coverage at or below this got no flux from the scatter
# (the `cov_cutoff` of `scatter_poly_fluxes_to_boundary!`).
const FLUX_SNAPSHOT_COV_CUTOFF = 1e-3

const _Face = Makie.GeometryBasics.GLTriangleFace

"""
    Plotting.plot_flux_snapshot(snapshot_file; fluxes, regions, outdir)

One figure per flux and region from a snapshot written by
`SimOutput.write_flux_snapshot`, with geometry read from
`exchange_grid_geometry.jld2` in the same directory. `regions` holds `nothing`
(global) and/or `(lat, lon, radius_km)` tuples. Returns the PNG paths.

- **Row 1:** the ocean side, grid by grid: exchange polygons, SE nodes, ocean cells.
- **Row 2:** what the atmosphere receives, plus panels reserved for sea ice and land.
- **Row 3:** consistency. How far each node lies outside its polygons' range (L2
  overshoot), area integrals per grid, and a panel reserved for the ocean–ice
  interface.

NaNs are drawn magenta.
"""
function Plotting.plot_flux_snapshot(
    snapshot_file::AbstractString;
    fluxes = ("F_sh", "F_lh", "F_moisture", "F_τu", "F_τv"),
    regions = (nothing,),
    outdir = joinpath(dirname(snapshot_file), "plots"),
)
    geometry_file = joinpath(dirname(snapshot_file), "exchange_grid_geometry.jld2")
    geo = JLD2.jldopen(_read_geometry, geometry_file, "r")
    mkpath(outdir)
    stem = splitext(basename(snapshot_file))[1]
    paths = String[]
    JLD2.jldopen(snapshot_file, "r") do snap
        for flux in fluxes, region in regions
            fig = _flux_snapshot_figure(geo, snap, string(flux), region)
            tag =
                isnothing(region) ? "global" :
                Printf.@sprintf("%.1fN_%.1fE", region[1], region[2])
            path = joinpath(outdir, "$(stem)_$(flux)_$(tag).png")
            Makie.save(path, fig)
            push!(paths, path)
        end
    end
    return paths
end

_read_geometry(file) = (;
    vert_ptr = file["poly/vert_ptr"],
    vert_lon = file["poly/vert_lon"],
    vert_lat = file["poly/vert_lat"],
    poly_area = file["poly/area"],
    node_lon = file["nodes/lon"],
    node_lat = file["nodes/lat"],
    node_Jw = file["nodes/Jw"],
    node_poly_ptr = file["nodes/poly_ptr"],
    node_poly = file["nodes/poly"],
    corner_lon = file["cells/corner_lon"],
    corner_lat = file["cells/corner_lat"],
    cell_lon = file["cells/lon"],
    cell_lat = file["cells/lat"],
    wet_area = file["cells/wet_area"],
)

_read(file, path) = haskey(file, path) ? file[path] : nothing

# The globe, or a box around (lat, lon).
struct _SnapshotView
    is_global::Bool
    lon0::Float64
    limits::NTuple{4, Float64} # lon_min, lon_max, lat_min, lat_max
    markersize::Float64
end

_snapshot_view(::Nothing) = _SnapshotView(true, 0.0, (-180.0, 180.0, -90.0, 90.0), 4.0)
function _snapshot_view(region)
    lat, lon, radius_km = region
    dlat = radius_km / 111.2
    dlon = dlat / max(cosd(lat), 0.1)
    return _SnapshotView(false, lon, (lon - dlon, lon + dlon, lat - dlat, lat + dlat), 7.0)
end

# `lon` shifted by a multiple of 360° to lie within 180° of the view center.
_wrap_lon(view, lon) = view.lon0 + mod(lon - view.lon0 + 180, 360) - 180

function _in_view(view, lon, lat)
    view.is_global && return true
    x = _wrap_lon(view, lon)
    return view.limits[1] <= x <= view.limits[2] && view.limits[3] <= lat <= view.limits[4]
end

# A ring's longitudes made continuous (no ±360° jumps between neighbors) and
# shifted near the view center, plus the longitude the ring returns to when it
# closes. That differs from the first vertex by ±360° when the ring winds around
# a pole.
function _ring_lon(view, lon)
    ring = collect(Float64, lon)
    for v in 2:length(ring)
        ring[v] += 360 * round((ring[v - 1] - ring[v]) / 360)
    end
    ring .+= _wrap_lon(view, ring[1]) - ring[1]
    step = ring[1] - ring[end]
    closing_lon = ring[end] + step - 360 * round(step / 360)
    return ring, closing_lon
end

function _polygon_rings(geo, view)
    ids = Int[]
    rings = Tuple{Vector{Float64}, Vector{Float64}}[]
    for k in 1:(length(geo.vert_ptr) - 1)
        r = geo.vert_ptr[k]:(geo.vert_ptr[k + 1] - 1)
        any(v -> _in_view(view, geo.vert_lon[v], geo.vert_lat[v]), r) || continue
        push!(ids, k)
        push!(rings, (geo.vert_lon[r], geo.vert_lat[r]))
    end
    return ids, rings
end

function _cell_rings(geo, view)
    ids = Int[]
    rings = Tuple{Vector{Float64}, Vector{Float64}}[]
    for c in eachindex(geo.wet_area)
        geo.wet_area[c] > 0 || continue
        _in_view(view, geo.cell_lon[c], geo.cell_lat[c]) || continue
        push!(ids, c)
        push!(rings, (geo.corner_lon[:, c], geo.corner_lat[:, c]))
    end
    return ids, rings
end

# Copies of a shape spanning longitudes `lo`..`hi` needed to cover a global map.
_lon_offsets(view, lo, hi) =
    view.is_global ? filter(o -> lo + o < 180 && hi + o > -180, (-360.0, 0.0, 360.0)) :
    (0.0,)

# Fan-triangulate one convex polygon (vertices as `(lon, lat)`) into `mesh`.
function _push_triangles!(mesh, value, vertices)
    base = length(mesh.points)
    for (x, y) in vertices
        push!(mesh.points, Makie.Point2f(x, y))
        push!(mesh.colors, value)
    end
    for v in 2:(length(vertices) - 1)
        push!(mesh.faces, _Face(base + 1, base + v, base + v + 1))
    end
    return nothing
end

# Flat-colored polygons as one triangle mesh. Convex rings are fan-triangulated.
# A ring around a pole becomes a strip from each edge down to the pole, which is
# what a polar cap looks like in longitude–latitude. Shapes that cross the map's
# edge are drawn again on the other side. Non-finite values go into a second,
# magenta mesh so NaNs stand out.
#
# `inactive` marks polygons whose flux is not applied (zero weight); they are drawn
# light gray rather than colored.
function _draw_polygons!(
    ax,
    view,
    rings,
    values,
    colorrange,
    colormap;
    inactive = falses(length(values)),
)
    good = (; points = Makie.Point2f[], faces = _Face[], colors = Float32[])
    bad = (; points = Makie.Point2f[], faces = _Face[], colors = Float32[])
    off = (; points = Makie.Point2f[], faces = _Face[], colors = Float32[])
    for (((lon, lat), value), skip) in zip(zip(rings, values), inactive)
        mesh = skip ? off : isfinite(value) ? good : bad
        x, closing_lon = _ring_lon(view, lon)
        if abs(closing_lon - x[1]) > 180
            pole = sum(lat) > 0 ? 90.0 : -90.0
            xs, ys = [x; closing_lon], [lat; lat[1]]
            lo, hi = extrema((xs[1], xs[end]))
            for o in _lon_offsets(view, lo, hi), i in 1:(length(xs) - 1)
                quad = [
                    (xs[i] + o, ys[i]),
                    (xs[i + 1] + o, ys[i + 1]),
                    (xs[i + 1] + o, pole),
                    (xs[i] + o, pole),
                ]
                _push_triangles!(mesh, value, quad)
            end
        else
            for o in _lon_offsets(view, extrema(x)...)
                _push_triangles!(mesh, value, [(xv + o, yv) for (xv, yv) in zip(x, lat)])
            end
        end
    end
    isempty(good.faces) || Makie.mesh!(
        ax,
        good.points,
        good.faces;
        color = good.colors,
        colorrange,
        colormap,
        shading = Makie.NoShading,
    )
    isempty(bad.faces) ||
        Makie.mesh!(ax, bad.points, bad.faces; color = :magenta, shading = Makie.NoShading)
    isempty(off.faces) ||
        Makie.mesh!(ax, off.points, off.faces; color = :gray85, shading = Makie.NoShading)
    return nothing
end

function _draw_outlines!(ax, view, rings)
    points = Makie.Point2f[]
    for (lon, lat) in rings
        x, closing_lon = _ring_lon(view, lon)
        append!(points, Makie.Point2f.(x, lat))
        push!(points, Makie.Point2f(closing_lon, lat[1]), Makie.Point2f(NaN, NaN))
    end
    Makie.lines!(ax, points; color = (:black, 0.35), linewidth = 0.5)
    return nothing
end

function _draw_nodes!(ax, view, geo, values, keep, colorrange, colormap; kwargs...)
    ids = [
        n for n in eachindex(values) if
        keep[n] && _in_view(view, geo.node_lon[n], geo.node_lat[n])
    ]
    good = filter(n -> isfinite(values[n]), ids)
    bad = filter(n -> !isfinite(values[n]), ids)
    xs(ns) = [_wrap_lon(view, geo.node_lon[n]) for n in ns]
    isempty(good) || Makie.scatter!(
        ax,
        xs(good),
        geo.node_lat[good];
        color = values[good],
        colorrange,
        colormap,
        markersize = view.markersize,
        strokewidth = 0,
        kwargs...,
    )
    isempty(bad) || Makie.scatter!(
        ax,
        xs(bad),
        geo.node_lat[bad];
        color = :magenta,
        markersize = 2 * view.markersize,
    )
    return nothing
end

# 1–99% limits of the finite values, symmetric if they change sign.
function _robust_limits(values)
    finite = sort!(filter(isfinite, collect(Float64, values)))
    isempty(finite) && return (-1.0, 1.0)
    pick(q) = finite[clamp(round(Int, q * length(finite)), 1, length(finite))]
    lo, hi = pick(0.01), pick(0.99)
    lo < 0 < hi && return (-max(-lo, hi), max(-lo, hi))
    lo == hi && return (lo - 1, hi + 1)
    return (lo, hi)
end

_snapshot_colormap(limits) = limits[1] < 0 < limits[2] ? :balance : :viridis

# How far each ocean node's value lies outside the range of the polygons that
# feed it (NaN where the node got no flux). Positive values are L2 overshoot.
function _node_excess(geo, node_values, cov, poly_values, weight)
    excess = fill(NaN, length(node_values))
    for n in eachindex(node_values)
        cov[n] > FLUX_SNAPSHOT_COV_CUTOFF || continue
        lo, hi = Inf, -Inf
        for p in geo.node_poly_ptr[n]:(geo.node_poly_ptr[n + 1] - 1)
            k = geo.node_poly[p]
            weight[k] > 0 || continue
            lo = min(lo, poly_values[k])
            hi = max(hi, poly_values[k])
        end
        lo <= hi && (excess[n] = max(node_values[n] - hi, lo - node_values[n], 0.0))
    end
    return excess
end

# Open-water-weighted area integral of the flux on each grid that has it.
function _integrals(geo, poly_values, weight, node_values, cov, cell_values)
    labels, values = String[], Float64[]
    if !isnothing(poly_values)
        push!(labels, "polygons")
        push!(values, sum(geo.poly_area .* weight .* poly_values))
    end
    if !isnothing(node_values)
        push!(labels, "SE nodes")
        covered = (
            geo.node_Jw[n] * cov[n] * node_values[n] for
            n in eachindex(node_values) if cov[n] > FLUX_SNAPSHOT_COV_CUTOFF
        )
        push!(values, sum(covered; init = 0.0))
    end
    if !isnothing(cell_values)
        push!(labels, "ocean cells")
        wet = (
            geo.wet_area[c] * cell_values[c] for
            c in eachindex(cell_values) if geo.wet_area[c] > 0
        )
        push!(values, sum(wet; init = 0.0))
    end
    return labels, values
end

function _placeholder!(position, title, message)
    ax = Makie.Axis(position; title)
    Makie.hidedecorations!(ax)
    Makie.text!(
        ax,
        0.5,
        0.5;
        text = message,
        align = (:center, :center),
        space = :relative,
        color = :gray50,
        fontsize = 16,
    )
    return ax
end

function _map_axis(position, view, title)
    ax = Makie.Axis(position; title, xlabel = "longitude", ylabel = "latitude")
    Makie.limits!(ax, view.limits...)
    return ax
end

function _flux_snapshot_figure(geo, snap, flux, region)
    view = _snapshot_view(region)
    units = get(FLUX_SNAPSHOT_UNITS, flux, "")
    poly_values = _read(snap, "ocean/poly/$flux")
    weight = _read(snap, "ocean/poly/weight")
    node_values = _read(snap, "ocean/nodes/$flux")
    cov = _read(snap, "ocean/nodes/cov")
    cell_values = _read(snap, "ocean/cells/$flux")
    atmos_values = _read(snap, "atmos/$flux")

    active_poly = isnothing(poly_values) ? nothing : poly_values[weight .> 0]
    colorrange = _robust_limits(something(active_poly, atmos_values, [0.0]))
    colormap = _snapshot_colormap(colorrange)
    # The atmosphere sees every surface, so it gets its own scale.
    atmos_range = isnothing(atmos_values) ? colorrange : _robust_limits(atmos_values)

    fig = Makie.Figure(size = (1800, 1350))
    place =
        isnothing(region) ? "global" :
        Printf.@sprintf("%.1f°N %.1f°E, radius %.0f km", region...)
    Makie.Label(
        fig[0, 1:4],
        "$flux [$units] at $(snap["date"]), $place";
        fontsize = 22,
        font = :bold,
    )

    # Row 1: the ocean side of the exchange, grid by grid.
    if isnothing(poly_values)
        _placeholder!(fig[1, 1], "exchange grid (ocean)", "not recorded")
    else
        title = "exchange grid: ocean polygons (gray: fully ice-covered)"
        ax = _map_axis(fig[1, 1], view, title)
        ids, rings = _polygon_rings(geo, view)
        inactive = weight[ids] .<= 0
        _draw_polygons!(ax, view, rings, poly_values[ids], colorrange, colormap; inactive)
        view.is_global || _draw_outlines!(ax, view, rings)
    end
    if isnothing(node_values)
        _placeholder!(fig[1, 2], "SE nodes (ocean)", "not recorded")
    else
        ax = _map_axis(fig[1, 2], view, "SE nodes: ocean only")
        covered = cov .> FLUX_SNAPSHOT_COV_CUTOFF
        _draw_nodes!(ax, view, geo, node_values, covered, colorrange, colormap)
    end
    if isnothing(cell_values)
        _placeholder!(fig[1, 3], "ocean cells", "not recorded")
    else
        ax = _map_axis(fig[1, 3], view, "ocean cells (open-water weighted)")
        ids, rings = _cell_rings(geo, view)
        _draw_polygons!(ax, view, rings, cell_values[ids], colorrange, colormap)
        view.is_global || _draw_outlines!(ax, view, rings)
    end
    Makie.Colorbar(fig[1, 4]; limits = colorrange, colormap, label = "$units (ocean)")

    # Row 2: what the atmosphere receives, and the surfaces not captured yet.
    if isnothing(atmos_values)
        _placeholder!(fig[2, 1], "SE nodes (atmosphere)", "not recorded")
    else
        ax = _map_axis(
            fig[2, 1],
            view,
            "SE nodes: all surfaces, as the atmosphere receives it",
        )
        keep = trues(length(atmos_values))
        atmos_colormap = _snapshot_colormap(atmos_range)
        _draw_nodes!(ax, view, geo, atmos_values, keep, atmos_range, atmos_colormap)
        Makie.Colorbar(
            fig[2, 4];
            limits = atmos_range,
            colormap = atmos_colormap,
            label = "$units (atmosphere)",
        )
    end
    _placeholder!(fig[2, 2], "sea ice", "not captured yet")
    _placeholder!(fig[2, 3], "land", "not captured yet")

    # Row 3: consistency between grids.
    if isnothing(poly_values) || isnothing(node_values)
        _placeholder!(fig[3, 1], "SE node outside its polygons' range", "not recorded")
    else
        excess = _node_excess(geo, node_values, cov, poly_values, weight)
        # Roundoff is not overshoot: ignore excess below a millionth of the flux's
        # magnitude. Nodes within range are drawn gray.
        scale = maximum(abs, filter(isfinite, poly_values); init = 0.0)
        floor = max(1e-6 * scale, floatmin(Float64))
        outside = [isfinite(e) && e > floor for e in excess]
        inside = [isfinite(e) && e <= floor for e in excess]
        emax = any(outside) ? maximum(excess[outside]) : 0.0
        title = Printf.@sprintf(
            "SE nodes outside their polygons' range: %d (max %.3g %s)",
            count(outside),
            emax,
            units,
        )
        ax = _map_axis(fig[3, 1], view, title)
        gray = [:gray80, :gray80]
        _draw_nodes!(ax, view, geo, zeros(length(excess)), inside, (0.0, 1.0), gray)
        excess_range = (max(floor, 1e-4 * emax), max(emax, 2 * floor))
        _draw_nodes!(
            ax,
            view,
            geo,
            excess,
            outside,
            excess_range,
            :inferno;
            colorscale = log10,
        )
        Makie.Colorbar(
            fig[3, 4];
            limits = excess_range,
            colormap = :inferno,
            scale = log10,
            label = "outside range ($units)",
        )
    end
    labels, values = _integrals(geo, poly_values, weight, node_values, cov, cell_values)
    if isempty(values)
        _placeholder!(fig[3, 2], "area integrals", "not recorded")
    else
        spread = (maximum(values) - minimum(values)) / max(abs(values[1]), eps())
        ax = Makie.Axis(
            fig[3, 2];
            title = Printf.@sprintf("area integrals (relative spread %.2g)", spread),
            xticks = (1:length(labels), labels),
            ylabel = "$units m²",
        )
        bar_labels = [Printf.@sprintf("%.4g", v) for v in values]
        Makie.barplot!(ax, 1:length(values), values; bar_labels)
    end
    _placeholder!(fig[3, 3], "ocean–ice interface", "not captured yet")
    return fig
end
