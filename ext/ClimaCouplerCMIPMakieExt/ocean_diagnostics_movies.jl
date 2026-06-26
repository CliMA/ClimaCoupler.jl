import Oceananigans.OutputReaders: FieldTimeSeries

const DEFAULT_OCEAN_MOVIE_PROJECTION = "+proj=eqearth"
const OCEAN_LAND_NAN_COLOR = :black

"""
    _ocean_surface_plot_kwargs(colormap, vmin, vmax)

Keyword arguments shared by ocean diagnostic heatmaps and globe surfaces.
"""
function _ocean_surface_plot_kwargs(colormap, vmin, vmax)
    return (; colormap, colorrange = (vmin, vmax), nan_color = OCEAN_LAND_NAN_COLOR)
end

"""
    _underlying_ocean_grid(grid)

Return the underlying Oceananigans grid, stripping an immersed-boundary wrapper if present.
"""
function _underlying_ocean_grid(grid)
    return grid isa OC.ImmersedBoundaries.ImmersedBoundaryGrid ? grid.underlying_grid : grid
end

"""
    _is_latlon_ocean_grid(grid)

Return `true` when the ocean grid is a regular `LatitudeLongitudeGrid`.
"""
function _is_latlon_ocean_grid(grid)
    return _underlying_ocean_grid(grid) isa OC.LatitudeLongitudeGrid
end

const DEFAULT_OCEAN_SURFACE_MOVIE_FIELDS = [
    ("sea_surface_temperature", "Sea surface temperature"),
    ("sea_surface_salinity", "Sea surface salinity"),
    ("sea_surface_height", "Sea surface height"),
    ("sea_surface_zonal_velocity", "Sea surface zonal velocity"),
    ("sea_surface_meridional_velocity", "Sea surface meridional velocity"),
    ("zonal_wind_stress", "Zonal wind stress"),
    ("meridional_wind_stress", "Meridional wind stress"),
    ("surface_heat_flux", "Surface heat flux"),
    ("surface_salinity_flux", "Surface salinity flux"),
]

# The 3-D fields writer stores T, S, u, v at all depths; surface slices duplicate
# `DEFAULT_OCEAN_SURFACE_MOVIE_FIELDS`. Pass `field_fields` explicitly for subsurface
# or non-overlapping variables (e.g. `vertical_velocity`).
const DEFAULT_OCEAN_FIELDS_MOVIE_FIELDS = Tuple{String, String}[]

# Surface field overlaid translucently on bathymetry in static diagnostic plots.
const DEFAULT_BATHYMETRY_OVERLAY_FIELD =
    ("sea_surface_temperature", "Sea surface temperature")

const DEFAULT_BATHYMETRY_OVERLAY_ALPHA = 0.55

const DEFAULT_SEAICE_SURFACE_MOVIE_FIELDS = [
    ("ice_concentration", "Ice concentration"),
    ("ice_thickness", "Ice thickness"),
    ("ice_zonal_velocity", "Ice zonal velocity"),
    ("ice_meridional_velocity", "Ice meridional velocity"),
    ("ice_top_temperature", "Ice top temperature"),
]

const DEFAULT_SEAICE_OVERLAY_BASE_FIELD = ("ice_thickness", "Ice thickness")
const DEFAULT_SEAICE_OVERLAY_FIELD = ("ice_top_temperature", "Ice top temperature")
const DEFAULT_SEAICE_MEAN_SNAPSHOT_FIELD = DEFAULT_SEAICE_OVERLAY_BASE_FIELD

# Colormap and color-range style for each diagnostic variable.
# `:sequential` uses the data extrema; `:divergent` uses a symmetric range about zero.
const OCEAN_MOVIE_PLOT_STYLES = Dict(
    "sea_surface_temperature" => (:thermal, :sequential),
    "temperature" => (:thermal, :sequential),
    "sea_surface_salinity" => (:haline, :sequential),
    "salinity" => (:haline, :sequential),
    "sea_surface_height" => (:balance, :divergent),
    "sea_surface_zonal_velocity" => (:balance, :divergent),
    "sea_surface_meridional_velocity" => (:balance, :divergent),
    "zonal_velocity" => (:balance, :divergent),
    "meridional_velocity" => (:balance, :divergent),
    "zonal_wind_stress" => (:curl, :divergent),
    "meridional_wind_stress" => (:curl, :divergent),
    "surface_heat_flux" => (:RdBu, :divergent),
    "surface_salinity_flux" => (:BrBG, :divergent),
    "bathymetry" => (:deep, :sequential),
    "ice_concentration" => (:ice, :sequential),
    "ice_thickness" => (:dense, :sequential),
    "ice_zonal_velocity" => (:balance, :divergent),
    "ice_meridional_velocity" => (:balance, :divergent),
    "ice_top_temperature" => (:thermal, :sequential),
)

"""
    _ocean_movie_plot_style(variable_name)

Return `(colormap, range_style)` for an ocean diagnostic variable.
Unknown variables default to a sequential viridis scale.
"""
function _ocean_movie_plot_style(variable_name)
    return get(OCEAN_MOVIE_PLOT_STYLES, variable_name, (:viridis, :sequential))
end

"""
    _ocean_part_jld2_paths(collection_path)

Return sorted paths to split JLD2 part files for a diagnostic collection, if any exist.
"""
function _ocean_part_jld2_paths(collection_path)
    dir = dirname(collection_path)
    base = basename(collection_path)
    part_number(f) = parse(Int, match(r"_part(\d+)\.jld2$", f).captures[1])
    part_files = filter(readdir(dir)) do f
        startswith(f, base * "_part") && endswith(f, ".jld2")
    end
    return joinpath.(dir, sort(part_files, by = part_number))
end

"""
    _ocean_diagnostic_field_time_series(collection_path, variable_name)

Load a diagnostic `FieldTimeSeries` for visualization only.

Boundary conditions are omitted (`boundary_conditions = nothing`) because serialized
metadata for some velocity fields (notably meridional velocity on tripolar grids)
can fail JLD2 deserialization and break halo filling when indexing snapshots.
Interior values are sufficient for plotting.

When split part files exist (``<collection>_part1.jld2``, etc.), all parts are merged.
Oceananigans only auto-discovers parts if the base ``<collection>.jld2`` is missing;
writers often leave a stale base file alongside parts, which would otherwise truncate
movies at the first split boundary (typically 10 days for a 15-day split interval).
"""
function _ocean_diagnostic_field_time_series(collection_path, variable_name)
    part_paths = _ocean_part_jld2_paths(collection_path)
    if !isempty(part_paths)
        return JLD2.jldopen(first(part_paths)) do file
            FieldTimeSeries(
                file,
                variable_name;
                boundary_conditions = nothing,
                Nparts = length(part_paths),
                part_paths = part_paths,
                path = first(part_paths),
            )
        end
    end
    return FieldTimeSeries(collection_path, variable_name; boundary_conditions = nothing)
end

"""
    _ocean_jld2_collection_path(output_path, collection, filename_prefix)

Return the JLD2 collection path (without extension) for an Oceananigans diagnostic
writer, or `nothing` if no matching files exist in `output_path`.
"""
function _ocean_jld2_collection_path(output_path, collection, filename_prefix)
    collection_path = joinpath(output_path, filename_prefix * "_" * collection)
    prefix = basename(collection_path)
    has_output =
        any(startswith(f, prefix) && endswith(f, ".jld2") for f in readdir(output_path))
    return has_output ? collection_path : nothing
end

"""
    _ocean_grid_coordinates(grid)

Return `(longitude, latitude)` matrices in degrees on the ocean model grid.
Longitudes are normalized to ``[-180, 180]``.
"""
function _ocean_grid_coordinates(grid)
    longitude, latitude = _ocean_grid_longitude_latitude(grid, OC.Center(), OC.Center())
    longitude = copy(longitude)
    longitude[longitude .> 180] .-= 360
    return longitude, latitude
end

"""
    _ocean_grid_longitude_latitude(grid, ℓx, ℓy)

Return raw `(longitude, latitude)` node matrices in degrees at field location
`(ℓx, ℓy)`, without normalizing longitude.
"""
function _ocean_grid_longitude_latitude(grid, ℓx, ℓy)
    underlying_grid = _underlying_ocean_grid(grid)
    λ = vec(Array(OC.λnodes(underlying_grid, ℓx(), ℓy(), OC.Center())))
    φ = vec(Array(OC.φnodes(underlying_grid, ℓx(), ℓy(), OC.Center())))
    longitude = repeat(λ, 1, length(φ))
    latitude = repeat(φ', length(λ), 1)
    return longitude, latitude
end

"""
    _periodic_longitude_globe(grid)

Return `true` when the ocean grid is periodic in longitude.
"""
function _periodic_longitude_globe(grid)
    return OC.topology(_underlying_ocean_grid(grid), 1) === OC.Periodic
end

"""
    _unwrap_globe_longitudes(longitude)

Unwrap a 2D longitude field so values vary smoothly along the first dimension.

When adjacent columns jump by more than 180°, add or subtract 360° for plotting
so `surface!` on a `GlobeAxis` does not tear at the antimeridian.
"""
function _unwrap_globe_longitudes(longitude)
    lons = copy(longitude)
    Nx, Ny = size(lons)
    for j in 1:Ny
        for i in 2:Nx
            Δλ = lons[i, j] - lons[i - 1, j]
            if Δλ > 180
                lons[i, j] -= 360
            elseif Δλ < -180
                lons[i, j] += 360
            end
        end
    end
    return lons
end

"""
    _prepare_globe_surface_coordinates(longitude, latitude, grid)

Return `(longitude, latitude)` prepared for globe `surface!`.
"""
function _prepare_globe_surface_coordinates(longitude, latitude, grid)
    if _periodic_longitude_globe(grid)
        longitude = _unwrap_globe_longitudes(longitude)
    end
    return longitude, latitude
end

"""
    _ocean_globe_coordinates(field; land_mask)

Return `(longitude, latitude, data)` in degrees for interactive globe plotting.
Longitudes are unwrapped along the periodic dimension when needed.
"""
function _ocean_globe_coordinates(field; land_mask)
    data = _field_to_heatmap_array(field; land_mask)
    ℓx, ℓy, _ = OC.location(field)
    longitude, latitude = _ocean_grid_longitude_latitude(field.grid, ℓx(), ℓy())
    longitude, latitude = _prepare_globe_surface_coordinates(longitude, latitude, field.grid)
    return longitude, latitude, data
end

"""
    _plot_ocean_globe!(ax, grid, longitude, latitude, data; colormap, vmin, vmax)

Plot ocean diagnostic data on a GeoMakie `GlobeAxis` via `surface!` with native
`(λ, φ)` nodes. Longitudes should already be unwrapped for periodic grids.
"""
function _plot_ocean_globe!(ax, grid, longitude, latitude, data; colormap, vmin, vmax)
    plot_kwargs = _ocean_surface_plot_kwargs(colormap, vmin, vmax)
    z = zeros(eltype(longitude), size(longitude))
    return surface!(
        ax,
        longitude,
        latitude,
        z;
        color = data,
        shading = NoShading,
        plot_kwargs...,
    )
end

"""
    _ocean_land_mask(grid)

Return a `BitMatrix` where `true` marks dry land cells on the ocean grid.
Uses the immersed-boundary mask at the top vertical level.
"""
function _ocean_land_mask(grid)
    Nx, Ny, Nz = size(grid)
    land = falses(Nx, Ny)
    for i in 1:Nx, j in 1:Ny
        land[i, j] = OC.ImmersedBoundaries.immersed_cell(i, j, Nz, grid)
    end
    return land
end

"""
    _field_to_heatmap_array(field; land_mask = nothing)

Extract a 2D array from an Oceananigans `Field` for heatmap plotting.
Surface diagnostics use the sole vertical level; 3D fields use the top level.
Land cells are set to `NaN` when `land_mask` is provided.
"""
function _field_to_heatmap_array(field; land_mask = nothing)
    _, _, nz = size(field)
    z = nz == 1 ? 1 : nz
    data = Array(OC.interior(field, :, :, z))
    if !isnothing(land_mask)
        data = copy(data)
        data[land_mask] .= NaN
    end
    return data
end

"""
    _format_simulation_time(t)

Format simulation time (in seconds) for plot titles.
"""
function _format_simulation_time(t)
    days = t / 86400
    return Printf.@sprintf("%.1f days", days)
end

"""
    _global_colorrange(fts; range_style = :sequential)

Return `(vmin, vmax)` over all time snapshots in `fts`.
For `:divergent` fields, return a symmetric range about zero.
"""
function _global_colorrange(fts; range_style = :sequential, land_mask = nothing)
    vmin = Inf
    vmax = -Inf
    for i in eachindex(fts.times)
        data = _field_to_heatmap_array(fts[i]; land_mask)
        for val in data
            isnan(val) && continue
            vmin = min(vmin, val)
            vmax = max(vmax, val)
        end
    end
    if range_style === :divergent
        lim = max(abs(vmin), abs(vmax))
        return (-lim, lim)
    else
        return (vmin, vmax)
    end
end

"""
    _array_colorrange(data; range_style = :sequential)

Return `(vmin, vmax)` for a single 2D array, ignoring `NaN` values.
"""
function _array_colorrange(data; range_style = :sequential)
    vmin = Inf
    vmax = -Inf
    for val in data
        isnan(val) && continue
        vmin = min(vmin, val)
        vmax = max(vmax, val)
    end
    if range_style === :divergent
        lim = max(abs(vmin), abs(vmax))
        return (-lim, lim)
    else
        return (vmin, vmax)
    end
end

"""
    _ocean_bathymetry_depth_array(grid)

Return a 2D array of ocean depth in meters (positive downward) on the model grid.
Dry land cells are `NaN`. Requires an `ImmersedBoundaryGrid` with a fitted bottom.
"""
function _ocean_bathymetry_depth_array(grid)
    grid isa OC.ImmersedBoundaryGrid ||
        error("Bathymetry plotting requires an `ImmersedBoundaryGrid`")
    land_mask = _ocean_land_mask(grid)
    bottom_height = grid.immersed_boundary.bottom_height
    return -_field_to_heatmap_array(bottom_height; land_mask)
end

"""
    _ocean_surface_time_mean_array(fts; land_mask)

Return the time mean of a surface diagnostic, ignoring `NaN` values cell-wise.
"""
function _ocean_surface_time_mean_array(fts; land_mask)
    sum_data = zeros(Float64, size(land_mask))
    counts = zeros(Int, size(land_mask))
    for i in eachindex(fts.times)
        data = _field_to_heatmap_array(fts[i]; land_mask)
        for idx in eachindex(data)
            val = data[idx]
            isnan(val) && continue
            sum_data[idx] += val
            counts[idx] += 1
        end
    end
    mean_data = fill(NaN, size(land_mask))
    for idx in eachindex(mean_data)
        counts[idx] > 0 && (mean_data[idx] = sum_data[idx] / counts[idx])
    end
    return mean_data
end

"""
    _plot_ocean_surface_color_layer!(ax, grid, data; colormap, vmin, vmax, alpha = 1.0)

Plot a flat 2D ocean surface field as colors on `ax`. Land is black when `alpha == 1`
and transparent when `alpha < 1` so an underlying bathymetry layer shows through.
"""
function _plot_ocean_surface_color_layer!(
    ax,
    grid,
    data;
    colormap,
    vmin,
    vmax,
    alpha = 1.0,
    land_nan_color = OCEAN_LAND_NAN_COLOR,
)
    plot_kwargs = (;
        colormap,
        colorrange = (vmin, vmax),
        nan_color = land_nan_color,
        alpha,
    )
    if _is_latlon_ocean_grid(grid)
        longitude, latitude = _ocean_grid_coordinates(grid)
        z = zeros(eltype(longitude), size(longitude))
        return surface!(
            ax,
            longitude,
            latitude,
            z;
            color = data,
            shading = NoShading,
            plot_kwargs...,
        )
    else
        return heatmap!(ax, data; plot_kwargs...)
    end
end

"""
    _plot_ocean_surface_snapshot!(
        fig,
        grid,
        data;
        title,
        colormap,
        vmin,
        vmax,
        dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
        colorbar_label,
    )

Plot a single 2D ocean surface field on `fig` and return the plot object.
"""
function _plot_ocean_surface_snapshot!(
    fig,
    grid,
    data;
    title,
    colormap,
    vmin,
    vmax,
    dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
    colorbar_label,
)
    plot_kwargs = _ocean_surface_plot_kwargs(colormap, vmin, vmax)
    if _is_latlon_ocean_grid(grid)
        longitude, latitude = _ocean_grid_coordinates(grid)
        ax = GeoMakie.GeoAxis(
            fig[1, 1];
            dest,
            title,
            xlabel = "Longitude (°)",
            ylabel = "Latitude (°)",
        )
        hm = CairoMakie.surface!(ax, longitude, latitude, data; plot_kwargs...)
    else
        ax = CairoMakie.Axis(
            fig[1, 1];
            title,
            xlabel = "i",
            ylabel = "j",
        )
        hm = CairoMakie.heatmap!(ax, data; plot_kwargs...)
    end
    CairoMakie.Colorbar(fig[1, 2], hm; label = colorbar_label)
    return hm
end

"""
    make_ocean_bathymetry_plot(
        collection_path::AbstractString,
        output_file::AbstractString;
        title = "Ocean depth",
        colormap = nothing,
        dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
        reference_variable = "sea_surface_temperature",
    )

Create a static map of ocean bathymetry (depth in meters) from the grid serialized
in Oceananigans JLD2 surface diagnostics. Bathymetry is time-invariant, so a single
snapshot is written instead of a movie. Land is masked in black.
"""
function make_ocean_bathymetry_plot(
    collection_path::AbstractString,
    output_file::AbstractString;
    title = "Ocean depth",
    colormap = nothing,
    dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
    reference_variable = "sea_surface_temperature",
)
    default_colormap, default_range_style = _ocean_movie_plot_style("bathymetry")
    colormap = isnothing(colormap) ? default_colormap : colormap

    fts = _ocean_diagnostic_field_time_series(collection_path, reference_variable)
    grid = fts.grid
    if !(grid isa OC.ImmersedBoundaryGrid)
        @warn "Skipping bathymetry plot: grid is not an ImmersedBoundaryGrid"
        return nothing
    end

    depth = _ocean_bathymetry_depth_array(grid)
    vmin, vmax = _array_colorrange(depth; range_style = default_range_style)
    if !isfinite(vmin) || !isfinite(vmax)
        @warn "Skipping bathymetry plot: no wet ocean cells"
        return nothing
    end
    vmin == vmax && (vmax = vmin + 1)

    fig = CairoMakie.Figure(; size = (900, 500))
    _plot_ocean_surface_snapshot!(
        fig,
        grid,
        depth;
        title,
        colormap,
        vmin,
        vmax,
        dest,
        colorbar_label = "depth (m)",
    )
    CairoMakie.save(output_file, fig)
    return output_file
end

"""
    _ocean_bathymetry_overlay_fields(
        collection_path::AbstractString,
        variable_name::AbstractString;
        variable_title = variable_name,
        title = nothing,
        time_index = nothing,
        depth_colormap = nothing,
        overlay_colormap = nothing,
    )

Load bathymetry and a surface diagnostic for overlay plotting. Returns a `NamedTuple`
or `nothing` when the data are unavailable.
"""
function _ocean_bathymetry_overlay_fields(
    collection_path::AbstractString,
    variable_name::AbstractString;
    variable_title = variable_name,
    title = nothing,
    time_index = nothing,
    depth_colormap = nothing,
    overlay_colormap = nothing,
)
    depth_colormap_default, depth_range_style = _ocean_movie_plot_style("bathymetry")
    overlay_colormap_default, overlay_range_style =
        _ocean_movie_plot_style(variable_name)
    depth_colormap =
        isnothing(depth_colormap) ? depth_colormap_default : depth_colormap
    overlay_colormap =
        isnothing(overlay_colormap) ? overlay_colormap_default : overlay_colormap

    fts = _ocean_diagnostic_field_time_series(collection_path, variable_name)
    isempty(fts.times) && begin
        @warn "Skipping bathymetry overlay for `$variable_name`: no time snapshots"
        return nothing
    end

    grid = fts.grid
    if !(grid isa OC.ImmersedBoundaryGrid)
        @warn "Skipping bathymetry overlay: grid is not an ImmersedBoundaryGrid"
        return nothing
    end

    land_mask = _ocean_land_mask(grid)
    depth = _ocean_bathymetry_depth_array(grid)
    depth_vmin, depth_vmax =
        _array_colorrange(depth; range_style = depth_range_style)
    if !isfinite(depth_vmin) || !isfinite(depth_vmax)
        @warn "Skipping bathymetry overlay: no wet ocean cells"
        return nothing
    end
    depth_vmin == depth_vmax && (depth_vmax = depth_vmin + 1)

    overlay_data = if isnothing(time_index)
        _ocean_surface_time_mean_array(fts; land_mask)
    else
        idx = clamp(time_index, 1, length(fts.times))
        _field_to_heatmap_array(fts[idx]; land_mask)
    end

    overlay_vmin, overlay_vmax =
        _global_colorrange(fts; range_style = overlay_range_style, land_mask)
    if !isfinite(overlay_vmin) || !isfinite(overlay_vmax)
        @warn "Skipping bathymetry overlay for `$variable_name`: no finite overlay values"
        return nothing
    end
    overlay_vmin == overlay_vmax && (overlay_vmax = overlay_vmin + 1)

    plot_title = if isnothing(title)
        if isnothing(time_index)
            "Ocean depth with mean $(variable_title)"
        else
            idx = clamp(time_index, 1, length(fts.times))
            "Ocean depth with $(variable_title), t = $(_format_simulation_time(fts.times[idx]))"
        end
    else
        title
    end

    return (;
        grid,
        depth,
        overlay_data,
        depth_vmin,
        depth_vmax,
        overlay_vmin,
        overlay_vmax,
        depth_colormap,
        overlay_colormap,
        variable_title,
        plot_title,
    )
end

"""
    make_ocean_bathymetry_overlay_plot(
        collection_path::AbstractString,
        variable_name::AbstractString,
        output_file::AbstractString;
        variable_title = variable_name,
        title = nothing,
        time_index = nothing,
        overlay_alpha = DEFAULT_BATHYMETRY_OVERLAY_ALPHA,
        depth_colormap = nothing,
        overlay_colormap = nothing,
        dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
    )

Create a static map with ocean bathymetry underneath and a translucent surface
diagnostic on top (default: time-mean sea surface temperature). Use `time_index` to
overlay a single snapshot instead of the time mean. Land is drawn in black on the
bathymetry layer; wet cells with missing overlay values remain unobscured.
"""
function make_ocean_bathymetry_overlay_plot(
    collection_path::AbstractString,
    variable_name::AbstractString,
    output_file::AbstractString;
    variable_title = variable_name,
    title = nothing,
    time_index = nothing,
    overlay_alpha = DEFAULT_BATHYMETRY_OVERLAY_ALPHA,
    depth_colormap = nothing,
    overlay_colormap = nothing,
    dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
)
    fields = _ocean_bathymetry_overlay_fields(
        collection_path,
        variable_name;
        variable_title,
        title,
        time_index,
        depth_colormap,
        overlay_colormap,
    )
    isnothing(fields) && return nothing

    fig = CairoMakie.Figure(; size = (1100, 500))
    if _is_latlon_ocean_grid(fields.grid)
        ax = GeoMakie.GeoAxis(
            fig[1, 1];
            dest,
            title = fields.plot_title,
            xlabel = "Longitude (°)",
            ylabel = "Latitude (°)",
        )
    else
        ax = CairoMakie.Axis(
            fig[1, 1];
            title = fields.plot_title,
            xlabel = "i",
            ylabel = "j",
        )
    end

    depth_plot = _plot_ocean_surface_color_layer!(
        ax,
        fields.grid,
        fields.depth;
        colormap = fields.depth_colormap,
        vmin = fields.depth_vmin,
        vmax = fields.depth_vmax,
        alpha = 1.0,
    )
    overlay_plot = _plot_ocean_surface_color_layer!(
        ax,
        fields.grid,
        fields.overlay_data;
        colormap = fields.overlay_colormap,
        vmin = fields.overlay_vmin,
        vmax = fields.overlay_vmax,
        alpha = overlay_alpha,
        land_nan_color = :transparent,
    )
    CairoMakie.Colorbar(fig[1, 2], depth_plot; label = "depth (m)")
    CairoMakie.Colorbar(fig[1, 3], overlay_plot; label = fields.variable_title)

    CairoMakie.save(output_file, fig)
    return output_file
end

"""
    view_ocean_bathymetry_overlay(
        collection_path::AbstractString,
        variable_name::AbstractString;
        variable_title = variable_name,
        title = nothing,
        time_index = nothing,
        overlay_alpha = DEFAULT_BATHYMETRY_OVERLAY_ALPHA,
        depth_colormap = nothing,
        overlay_colormap = nothing,
        dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
        wait = true,
    )

Open an interactive map with ocean bathymetry underneath and a translucent surface
diagnostic on top. Use the opacity slider to fade the overlay in and out; land is
drawn in black on the bathymetry layer.

Requires `GLMakie`. By default the call blocks until the window is closed
(`wait = false` returns immediately, e.g. from the REPL).

Returns `(fig, ax, depth_plot, overlay_plot)`.
"""
function Plotting.view_ocean_bathymetry_overlay(
    collection_path::AbstractString,
    variable_name::AbstractString;
    variable_title = variable_name,
    title = nothing,
    time_index = nothing,
    overlay_alpha = DEFAULT_BATHYMETRY_OVERLAY_ALPHA,
    depth_colormap = nothing,
    overlay_colormap = nothing,
    dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
    wait = true,
)
    GLM = _glmakie_module()
    GLM.activate!()

    fields = _ocean_bathymetry_overlay_fields(
        collection_path,
        variable_name;
        variable_title,
        title,
        time_index,
        depth_colormap,
        overlay_colormap,
    )
    isnothing(fields) && return nothing

    fig = GLM.Figure(; size = (1150, 580))
    if _is_latlon_ocean_grid(fields.grid)
        ax = GeoMakie.GeoAxis(
            fig[1, 1];
            dest,
            title = fields.plot_title,
            xlabel = "Longitude (°)",
            ylabel = "Latitude (°)",
        )
    else
        ax = GLM.Axis(
            fig[1, 1];
            title = fields.plot_title,
            xlabel = "i",
            ylabel = "j",
        )
    end

    depth_plot = _plot_ocean_surface_color_layer!(
        ax,
        fields.grid,
        fields.depth;
        colormap = fields.depth_colormap,
        vmin = fields.depth_vmin,
        vmax = fields.depth_vmax,
        alpha = 1.0,
    )
    overlay_alpha_obs = Observable(overlay_alpha)
    overlay_plot = _plot_ocean_surface_color_layer!(
        ax,
        fields.grid,
        fields.overlay_data;
        colormap = fields.overlay_colormap,
        vmin = fields.overlay_vmin,
        vmax = fields.overlay_vmax,
        alpha = overlay_alpha_obs,
        land_nan_color = :transparent,
    )
    GLM.Colorbar(fig[1, 2], depth_plot; label = "depth (m)")
    GLM.Colorbar(fig[1, 3], overlay_plot; label = fields.variable_title)

    alpha_slider = GLM.Slider(
        fig[2, 1],
        range = 0:0.01:1,
        startvalue = overlay_alpha,
        tellheight = false,
    )
    GLM.Label(
        fig[2, 1],
        lift(alpha_slider.value) do alpha
            "overlay opacity: $(Printf.@sprintf("%.2f", alpha))"
        end,
        tellwidth = false,
    )
    on(alpha_slider.value) do alpha
        overlay_alpha_obs[] = alpha
    end

    screen = display(fig)
    _wait_for_glmakie_screen(
        GLM,
        screen;
        wait,
        message = "Interactive bathymetry overlay open — close the window to exit (Ctrl+C also works)",
    )
    return fig, ax, depth_plot, overlay_plot
end

"""
    make_ocean_field_movie(
        collection_path::AbstractString,
        variable_name::AbstractString,
        output_file::AbstractString;
        title = variable_name,
        framerate = 4,
        format = "mp4",
    )

Create a movie for a single Oceananigans diagnostic variable stored in JLD2 output.
`collection_path` is the path prefix (without `.jld2`) passed to `FieldTimeSeries`,
for example `joinpath(output_dir, "ocean_surface")`.

Colormap and color-range style default to values in `OCEAN_MOVIE_PLOT_STYLES` for the
given `variable_name`, but can be overridden with `colormap` and `range_style`.

Tripolar and other curvilinear grids are plotted in native model-index space
(`heatmap!` with ``i``/``j`` axes). Regular latitude-longitude grids use
`GeoMakie.GeoAxis` with `surface!`. For an interactive 3D globe view of the same data, use `Plotting.view_ocean_field_globe` or
`experiments/CMIP/view_ocean_globe.jl`. All ocean diagnostic plots mask land
in black and include a colorbar.
"""
function make_ocean_field_movie(
    collection_path::AbstractString,
    variable_name::AbstractString,
    output_file::AbstractString;
    title = variable_name,
    colormap = nothing,
    range_style = nothing,
    dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
    framerate = 4,
    format = "mp4",
)
    default_colormap, default_range_style = _ocean_movie_plot_style(variable_name)
    colormap = isnothing(colormap) ? default_colormap : colormap
    range_style = isnothing(range_style) ? default_range_style : range_style

    fts = _ocean_diagnostic_field_time_series(collection_path, variable_name)
    length(fts.times) < 2 && begin
        @warn "Skipping movie for `$variable_name`: fewer than 2 time snapshots"
        return nothing
    end

    land_mask = _ocean_land_mask(fts.grid)
    vmin, vmax = _global_colorrange(fts; range_style, land_mask)
    if vmin == vmax
        @warn "Skipping movie for `$variable_name`: field is spatially constant"
        return nothing
    end

    data = Observable(_field_to_heatmap_array(fts[1]; land_mask))
    fig = CairoMakie.Figure(; size = (900, 500))
    plot_title = "$title, t = $(_format_simulation_time(fts.times[1]))"

    if _is_latlon_ocean_grid(fts.grid)
        longitude, latitude = _ocean_grid_coordinates(fts.grid)
        ax = GeoMakie.GeoAxis(
            fig[1, 1];
            dest,
            title = plot_title,
            xlabel = "Longitude (°)",
            ylabel = "Latitude (°)",
        )
        hm = CairoMakie.surface!(
            ax,
            longitude,
            latitude,
            data;
            _ocean_surface_plot_kwargs(colormap, vmin, vmax)...,
        )
    else
        ax = CairoMakie.Axis(
            fig[1, 1];
            title = plot_title,
            xlabel = "i",
            ylabel = "j",
        )
        hm = CairoMakie.heatmap!(
            ax,
            data;
            _ocean_surface_plot_kwargs(colormap, vmin, vmax)...,
        )
    end

    CairoMakie.Colorbar(fig[1, 2], hm; label = title)

    time_indices = eachindex(fts.times)
    Makie.record(
        fig,
        output_file,
        time_indices;
        framerate,
        format,
    ) do i
        data[] = _field_to_heatmap_array(fts[i]; land_mask)
        ax.title[] = "$title, t = $(_format_simulation_time(fts.times[i]))"
    end

    return output_file
end

"""
    Plotting.debug_plot!(ax, fig, field::OC.Field, i, j)

Plot a heatmap of the provided Oceananigans field with a colorbar.
Land cells on immersed-boundary grids are masked in black.
"""
function Plotting.debug_plot!(ax, fig, field::OC.Field, i, j)
    land_mask = _ocean_land_mask(field.grid)
    data = _field_to_heatmap_array(field; land_mask)
    hm = CairoMakie.heatmap!(ax, data; nan_color = OCEAN_LAND_NAN_COLOR)
    Makie.Colorbar(fig[i, j * 2], hm)
    return nothing
end

function Plotting.debug_plot!(ax, fig, field::OC.AbstractOperations.AbstractOperation, i, j)
    evaluated_field = OC.Field(field)
    OC.compute!(evaluated_field)
    Plotting.debug_plot!(ax, fig, evaluated_field, i, j)
    return nothing
end

"""
    _glmakie_module()

Return the loaded `GLMakie` module, or throw if it is not available.
"""
function _glmakie_module()
    if isdefined(Main, :GLMakie)
        return getfield(Main, :GLMakie)
    end
    for (pkgid, mod) in Base.loaded_modules
        if pkgid.name == :GLMakie
            return mod
        end
    end
    error("Run `using GLMakie` before calling interactive ocean diagnostic viewers.")
end

"""
    _wait_for_glmakie_screen(GLM, screen; wait = true, message = ...)

Block until the GLMakie window is closed. Required when running from a script;
otherwise the Julia process exits immediately and the window disappears.
"""
function _wait_for_glmakie_screen(
    GLM,
    screen;
    wait = true,
    message = "Interactive viewer open — close the window to exit (Ctrl+C also works)",
)
    wait || return screen
    @info message
    while GLM.isopen(screen)
        sleep(0.05)
    end
    return screen
end

"""
    view_ocean_field_globe(
        collection_path::AbstractString,
        variable_name::AbstractString;
        title = variable_name,
        time_index = 1,
        colormap = nothing,
        range_style = nothing,
    )

Open an interactive 3D globe of an Oceananigans diagnostic time series.

Requires `GLMakie`. Data are plotted on a GeoMakie `GlobeAxis` with `surface!` at
native `(λ, φ)` nodes; longitudes are unwrapped across the periodic boundary before
plotting. Immersed-boundary land is drawn in black, and a colorbar is shown. Use
the time slider to scrub snapshots; rotate and zoom with the mouse. By default the
call blocks until the window is closed (set `wait = false` to return immediately,
e.g. from the REPL).

Returns `(fig, ax, plt)`.
"""
function Plotting.view_ocean_field_globe(
    collection_path::AbstractString,
    variable_name::AbstractString;
    title = variable_name,
    time_index = 1,
    colormap = nothing,
    range_style = nothing,
    wait = true,
)
    GLM = _glmakie_module()
    GLM.activate!()

    default_colormap, default_range_style = _ocean_movie_plot_style(variable_name)
    colormap = isnothing(colormap) ? default_colormap : colormap
    range_style = isnothing(range_style) ? default_range_style : range_style

    fts = _ocean_diagnostic_field_time_series(collection_path, variable_name)
    Nt = length(fts.times)
    time_index = clamp(time_index, 1, Nt)

    land_mask = _ocean_land_mask(fts.grid)
    vmin, vmax = _global_colorrange(fts; range_style, land_mask)

    longitude, latitude, data0 = _ocean_globe_coordinates(fts[time_index]; land_mask)

    time_idx = Observable(time_index)
    data_obs = Observable(data0)
    plot_title = lift(time_idx) do idx
        "$(title), t = $(_format_simulation_time(fts.times[idx]))"
    end

    fig = GLM.Figure(; size = (950, 780))
    ax = GlobeAxis(
        fig[1, 1];
        show_axis = false,
        backgroundvisible = false,
        title = plot_title,
    )
    plt = _plot_ocean_globe!(
        ax,
        fts.grid,
        longitude,
        latitude,
        data_obs;
        colormap,
        vmin,
        vmax,
    )
    GLM.Colorbar(fig[1, 2], plt; label = title)

    slider = GLM.Slider(
        fig[2, 1],
        range = 1:Nt,
        startvalue = time_index,
        tellheight = false,
    )
    GLM.Label(
        fig[2, 1],
        lift(time_idx) do idx
            "snapshot $(idx) / $(Nt)"
        end,
        tellwidth = false,
    )
    on(slider.value) do val
        idx = Int(round(val))
        time_idx[] = idx
        _, _, data_obs[] = _ocean_globe_coordinates(fts[idx]; land_mask)
    end

    screen = display(fig)
    _wait_for_glmakie_screen(
        GLM,
        screen;
        wait,
        message = "Interactive globe open — close the window to exit (Ctrl+C also works)",
    )
    return fig, ax, plt
end

"""
    make_ocean_diagnostics_movies(
        output_path::AbstractString,
        plot_path::AbstractString;
        output_prefix = "ocean_",
        filename_prefix = "ocean",
        surface_fields = DEFAULT_OCEAN_SURFACE_MOVIE_FIELDS,
        field_fields = DEFAULT_OCEAN_FIELDS_MOVIE_FIELDS,
        framerate = 4,
        format = "mp4",
    )

Create movies for Oceananigans ocean diagnostics written by `add_ocean_diagnostics!`.

Reads JLD2 collections in `output_path`:
- `<filename_prefix>_surface[...].jld2` for 2D surface diagnostics (default movies)
- `<filename_prefix>_fields[...].jld2` only when `field_fields` is non-empty

Movies are saved to `plot_path` as `<output_prefix><variable_name>.<format>`.
Static bathymetry maps are also saved when surface diagnostics are available:
- `<output_prefix>bathymetry.png`
- `<output_prefix>bathymetry_<variable_name>.png` for a translucent overlay of the
  time-mean surface field on bathymetry (default: sea surface temperature). For an
  interactive overlay with an opacity slider, use `Plotting.view_ocean_bathymetry_overlay`
  or `experiments/CMIP/view_ocean_bathymetry_overlay.jl`.
Returns a vector of paths to the created movie files and static plots.
"""
function Plotting.make_ocean_diagnostics_movies(
    output_path::AbstractString,
    plot_path::AbstractString;
    output_prefix = "ocean_",
    filename_prefix = "ocean",
    surface_fields = DEFAULT_OCEAN_SURFACE_MOVIE_FIELDS,
    field_fields = DEFAULT_OCEAN_FIELDS_MOVIE_FIELDS,
    dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
    framerate = 4,
    format = "mp4",
)
    mkpath(plot_path)
    output_files = String[]

    surface_path = _ocean_jld2_collection_path(output_path, "surface", filename_prefix)
    if !isnothing(surface_path)
        bathymetry_file = joinpath(plot_path, output_prefix * "bathymetry.png")
        try
            result = make_ocean_bathymetry_plot(surface_path, bathymetry_file; dest)
            !isnothing(result) && push!(output_files, result)
        catch err
            @warn "Failed to create bathymetry plot" exception = err
        end

        overlay_variable, overlay_title = DEFAULT_BATHYMETRY_OVERLAY_FIELD
        overlay_file = joinpath(
            plot_path,
            output_prefix * "bathymetry_" * overlay_variable * ".png",
        )
        try
            result = make_ocean_bathymetry_overlay_plot(
                surface_path,
                overlay_variable,
                overlay_file;
                variable_title = overlay_title,
                dest,
            )
            !isnothing(result) && push!(output_files, result)
        catch err
            @warn "Failed to create bathymetry overlay plot for `$overlay_variable`" exception =
                err
        end

        for (variable_name, title) in surface_fields
            output_file =
                joinpath(plot_path, output_prefix * variable_name * "." * format)
            try
                result = make_ocean_field_movie(
                    surface_path,
                    variable_name,
                    output_file;
                    title,
                    dest,
                    framerate,
                    format,
                )
                !isnothing(result) && push!(output_files, result)
            catch err
                @warn "Failed to create movie for surface field `$variable_name`" exception =
                    err
            end
        end
    end

    fields_path = _ocean_jld2_collection_path(output_path, "fields", filename_prefix)
    if !isnothing(fields_path)
        for (variable_name, title) in field_fields
            output_file =
                joinpath(plot_path, output_prefix * variable_name * "." * format)
            try
                result = make_ocean_field_movie(
                    fields_path,
                    variable_name,
                    output_file;
                    title,
                    dest,
                    framerate,
                    format,
                )
                !isnothing(result) && push!(output_files, result)
            catch err
                @warn "Failed to create movie for 3D field `$variable_name`" exception = err
            end
        end
    end

    if isempty(output_files)
        @info "No ocean diagnostic movies or bathymetry plot created in `$output_path`"
    else
        @info "Created $(length(output_files)) ocean diagnostic output file(s) in `$plot_path`"
    end

    return output_files
end

"""
    _field_time_aggregate(fts, land_mask, time_index = nothing)

Return a 2D field array from `fts`, using the time mean when `time_index` is `nothing`.
"""
function _field_time_aggregate(fts, land_mask, time_index = nothing)
    if isnothing(time_index)
        return _ocean_surface_time_mean_array(fts; land_mask)
    end
    idx = clamp(time_index, 1, length(fts.times))
    return _field_to_heatmap_array(fts[idx]; land_mask)
end

"""
    _surface_field_overlay_fields(
        collection_path::AbstractString,
        base_variable::AbstractString,
        overlay_variable::AbstractString;
        base_title = base_variable,
        overlay_title = overlay_variable,
        title = nothing,
        base_time_index = nothing,
        overlay_time_index = nothing,
        base_colormap = nothing,
        overlay_colormap = nothing,
    )

Load two surface diagnostics from the same JLD2 collection for overlay plotting.
"""
function _surface_field_overlay_fields(
    collection_path::AbstractString,
    base_variable::AbstractString,
    overlay_variable::AbstractString;
    base_title = base_variable,
    overlay_title = overlay_variable,
    title = nothing,
    base_time_index = nothing,
    overlay_time_index = nothing,
    base_colormap = nothing,
    overlay_colormap = nothing,
)
    base_fts = _ocean_diagnostic_field_time_series(collection_path, base_variable)
    overlay_fts = _ocean_diagnostic_field_time_series(collection_path, overlay_variable)

    isempty(base_fts.times) && begin
        @warn "Skipping overlay: no snapshots for base field `$base_variable`"
        return nothing
    end
    isempty(overlay_fts.times) && begin
        @warn "Skipping overlay: no snapshots for overlay field `$overlay_variable`"
        return nothing
    end

    grid = base_fts.grid
    land_mask = _ocean_land_mask(grid)

    base_colormap_default, base_range_style = _ocean_movie_plot_style(base_variable)
    overlay_colormap_default, overlay_range_style =
        _ocean_movie_plot_style(overlay_variable)
    base_colormap = isnothing(base_colormap) ? base_colormap_default : base_colormap
    overlay_colormap =
        isnothing(overlay_colormap) ? overlay_colormap_default : overlay_colormap

    base_data = _field_time_aggregate(base_fts, land_mask, base_time_index)
    overlay_data = _field_time_aggregate(overlay_fts, land_mask, overlay_time_index)

    base_vmin, base_vmax =
        _array_colorrange(base_data; range_style = base_range_style)
    if !isfinite(base_vmin) || !isfinite(base_vmax)
        @warn "Skipping overlay: no finite values for base field `$base_variable`"
        return nothing
    end
    base_vmin == base_vmax && (base_vmax = base_vmin + 1)

    overlay_vmin, overlay_vmax = _global_colorrange(
        overlay_fts;
        range_style = overlay_range_style,
        land_mask,
    )
    if !isfinite(overlay_vmin) || !isfinite(overlay_vmax)
        @warn "Skipping overlay: no finite values for overlay field `$overlay_variable`"
        return nothing
    end
    overlay_vmin == overlay_vmax && (overlay_vmax = overlay_vmin + 1)

    plot_title = if isnothing(title)
        if isnothing(base_time_index) && isnothing(overlay_time_index)
            "$(base_title) with mean $(overlay_title)"
        elseif isnothing(base_time_index) && !isnothing(overlay_time_index)
            "Mean $(base_title) with $(overlay_title), t = $(_format_simulation_time(overlay_fts.times[clamp(overlay_time_index, 1, length(overlay_fts.times))]))"
        elseif !isnothing(base_time_index) && isnothing(overlay_time_index)
            "$(base_title) with mean $(overlay_title), t = $(_format_simulation_time(base_fts.times[clamp(base_time_index, 1, length(base_fts.times))]))"
        else
            base_idx = clamp(base_time_index, 1, length(base_fts.times))
            overlay_idx = clamp(overlay_time_index, 1, length(overlay_fts.times))
            "$(base_title) and $(overlay_title), t = $(_format_simulation_time(base_fts.times[base_idx]))"
        end
    else
        title
    end

    return (;
        grid,
        base_data,
        overlay_data,
        base_vmin,
        base_vmax,
        overlay_vmin,
        overlay_vmax,
        base_colormap,
        overlay_colormap,
        base_title,
        overlay_title,
        plot_title,
    )
end

"""
    _surface_overlay_axis(fig, grid, plot_title; dest)

Create a map axis for a two-layer surface overlay figure.
"""
function _surface_overlay_axis(fig, grid, plot_title; dest = DEFAULT_OCEAN_MOVIE_PROJECTION)
    if _is_latlon_ocean_grid(grid)
        return GeoMakie.GeoAxis(
            fig[1, 1];
            dest,
            title = plot_title,
            xlabel = "Longitude (°)",
            ylabel = "Latitude (°)",
        )
    else
        return Axis(
            fig[1, 1];
            title = plot_title,
            xlabel = "i",
            ylabel = "j",
        )
    end
end

"""
    _plot_surface_overlay_layers!(
        ax,
        fields;
        overlay_alpha,
        overlay_alpha_obs = nothing,
    )

Plot base and overlay layers on `ax`. Returns `(base_plot, overlay_plot)`.
"""
function _plot_surface_overlay_layers!(
    ax,
    fields;
    overlay_alpha,
    overlay_alpha_obs = nothing,
)
    base_plot = _plot_ocean_surface_color_layer!(
        ax,
        fields.grid,
        fields.base_data;
        colormap = fields.base_colormap,
        vmin = fields.base_vmin,
        vmax = fields.base_vmax,
        alpha = 1.0,
    )
    overlay_plot = _plot_ocean_surface_color_layer!(
        ax,
        fields.grid,
        fields.overlay_data;
        colormap = fields.overlay_colormap,
        vmin = fields.overlay_vmin,
        vmax = fields.overlay_vmax,
        alpha = something(overlay_alpha_obs, overlay_alpha),
        land_nan_color = :transparent,
    )
    return base_plot, overlay_plot
end

"""
    make_surface_field_overlay_plot(
        collection_path::AbstractString,
        base_variable::AbstractString,
        overlay_variable::AbstractString,
        output_file::AbstractString;
        base_title = base_variable,
        overlay_title = overlay_variable,
        kwargs...
    )

Create a static map with a base surface diagnostic and translucent overlay field.
"""
function make_surface_field_overlay_plot(
    collection_path::AbstractString,
    base_variable::AbstractString,
    overlay_variable::AbstractString,
    output_file::AbstractString;
    base_title = base_variable,
    overlay_title = overlay_variable,
    overlay_alpha = DEFAULT_BATHYMETRY_OVERLAY_ALPHA,
    dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
    kwargs...,
)
    fields = _surface_field_overlay_fields(
        collection_path,
        base_variable,
        overlay_variable;
        base_title,
        overlay_title,
        kwargs...,
    )
    isnothing(fields) && return nothing

    fig = CairoMakie.Figure(; size = (1100, 500))
    ax = _surface_overlay_axis(fig, fields.grid, fields.plot_title; dest)
    base_plot, overlay_plot =
        _plot_surface_overlay_layers!(ax, fields; overlay_alpha = overlay_alpha)
    CairoMakie.Colorbar(fig[1, 2], base_plot; label = fields.base_title)
    CairoMakie.Colorbar(fig[1, 3], overlay_plot; label = fields.overlay_title)
    CairoMakie.save(output_file, fig)
    return output_file
end

"""
    make_mean_field_snapshot_plot(
        collection_path::AbstractString,
        variable_name::AbstractString,
        output_file::AbstractString;
        variable_title = variable_name,
        title = nothing,
        colormap = nothing,
        dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
    )

Create a static map of the time-mean of a surface diagnostic.
"""
function make_mean_field_snapshot_plot(
    collection_path::AbstractString,
    variable_name::AbstractString,
    output_file::AbstractString;
    variable_title = variable_name,
    title = nothing,
    colormap = nothing,
    dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
)
    default_colormap, range_style = _ocean_movie_plot_style(variable_name)
    colormap = isnothing(colormap) ? default_colormap : colormap

    fts = _ocean_diagnostic_field_time_series(collection_path, variable_name)
    isempty(fts.times) && begin
        @warn "Skipping mean snapshot for `$variable_name`: no time snapshots"
        return nothing
    end

    land_mask = _ocean_land_mask(fts.grid)
    data = _ocean_surface_time_mean_array(fts; land_mask)
    vmin, vmax = _global_colorrange(fts; range_style, land_mask)
    if !isfinite(vmin) || !isfinite(vmax)
        @warn "Skipping mean snapshot for `$variable_name`: no finite values"
        return nothing
    end
    vmin == vmax && (vmax = vmin + 1)

    plot_title = isnothing(title) ? "Mean $(variable_title)" : title
    fig = CairoMakie.Figure(; size = (900, 500))
    _plot_ocean_surface_snapshot!(
        fig,
        fts.grid,
        data;
        title = plot_title,
        colormap,
        vmin,
        vmax,
        dest,
        colorbar_label = variable_title,
    )
    CairoMakie.save(output_file, fig)
    return output_file
end

"""
    view_surface_field_overlay(
        collection_path::AbstractString,
        base_variable::AbstractString,
        overlay_variable::AbstractString;
        base_title = base_variable,
        overlay_title = overlay_variable,
        overlay_alpha = DEFAULT_BATHYMETRY_OVERLAY_ALPHA,
        dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
        wait = true,
        kwargs...
    )

Open an interactive two-layer surface diagnostic map with an opacity slider.
"""
function Plotting.view_surface_field_overlay(
    collection_path::AbstractString,
    base_variable::AbstractString,
    overlay_variable::AbstractString;
    base_title = base_variable,
    overlay_title = overlay_variable,
    overlay_alpha = DEFAULT_BATHYMETRY_OVERLAY_ALPHA,
    dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
    wait = true,
    kwargs...,
)
    GLM = _glmakie_module()
    GLM.activate!()

    fields = _surface_field_overlay_fields(
        collection_path,
        base_variable,
        overlay_variable;
        base_title,
        overlay_title,
        kwargs...,
    )
    isnothing(fields) && return nothing

    fig = GLM.Figure(; size = (1150, 580))
    ax = _surface_overlay_axis(fig, fields.grid, fields.plot_title; dest)
    overlay_alpha_obs = Observable(overlay_alpha)
    base_plot, overlay_plot = _plot_surface_overlay_layers!(
        ax,
        fields;
        overlay_alpha,
        overlay_alpha_obs,
    )
    GLM.Colorbar(fig[1, 2], base_plot; label = fields.base_title)
    GLM.Colorbar(fig[1, 3], overlay_plot; label = fields.overlay_title)

    alpha_slider = GLM.Slider(
        fig[2, 1],
        range = 0:0.01:1,
        startvalue = overlay_alpha,
        tellheight = false,
    )
    GLM.Label(
        fig[2, 1],
        lift(alpha_slider.value) do alpha
            "overlay opacity: $(Printf.@sprintf("%.2f", alpha))"
        end,
        tellwidth = false,
    )
    on(alpha_slider.value) do alpha
        overlay_alpha_obs[] = alpha
    end

    screen = display(fig)
    _wait_for_glmakie_screen(
        GLM,
        screen;
        wait,
        message = "Interactive overlay open — close the window to exit (Ctrl+C also works)",
    )
    return fig, ax, base_plot, overlay_plot
end

"""
    view_seaice_field_overlay(collection_path; kwargs...)

Interactive sea-ice overlay viewer with default base `ice_thickness` and overlay
`ice_top_temperature`.
"""
function Plotting.view_seaice_field_overlay(
    collection_path::AbstractString;
    base_variable = first(DEFAULT_SEAICE_OVERLAY_BASE_FIELD),
    overlay_variable = first(DEFAULT_SEAICE_OVERLAY_FIELD),
    base_title = last(DEFAULT_SEAICE_OVERLAY_BASE_FIELD),
    overlay_title = last(DEFAULT_SEAICE_OVERLAY_FIELD),
    kwargs...,
)
    return Plotting.view_surface_field_overlay(
        collection_path,
        base_variable,
        overlay_variable;
        base_title,
        overlay_title,
        kwargs...,
    )
end

"""
    make_seaice_diagnostics_movies(
        output_path::AbstractString,
        plot_path::AbstractString;
        output_prefix = "seaice_",
        filename_prefix = "seaice",
        surface_fields = DEFAULT_SEAICE_SURFACE_MOVIE_FIELDS,
        framerate = 4,
        format = "mp4",
    )

Create movies and static maps for ClimaSeaIce diagnostics written by `add_seaice_diagnostics!`.

Reads `seaice_surface*.jld2` in `output_path` and writes:
- `<output_prefix><variable_name>.<format>` movies for each surface field
- `<output_prefix>mean_<base_variable>.png` time-mean snapshot of the base overlay field
- `<output_prefix><base_variable>_<overlay_variable>.png` translucent overlay plot
"""
function Plotting.make_seaice_diagnostics_movies(
    output_path::AbstractString,
    plot_path::AbstractString;
    output_prefix = "seaice_",
    filename_prefix = "seaice",
    surface_fields = DEFAULT_SEAICE_SURFACE_MOVIE_FIELDS,
    dest = DEFAULT_OCEAN_MOVIE_PROJECTION,
    framerate = 4,
    format = "mp4",
)
    mkpath(plot_path)
    output_files = String[]

    surface_path = _ocean_jld2_collection_path(output_path, "surface", filename_prefix)
    if isnothing(surface_path)
        @info "No sea-ice diagnostic movies created in `$output_path`"
        return output_files
    end

    mean_variable, mean_title = DEFAULT_SEAICE_MEAN_SNAPSHOT_FIELD
    mean_file = joinpath(plot_path, output_prefix * "mean_" * mean_variable * ".png")
    try
        result = make_mean_field_snapshot_plot(
            surface_path,
            mean_variable,
            mean_file;
            variable_title = mean_title,
            dest,
        )
        !isnothing(result) && push!(output_files, result)
    catch err
        @warn "Failed to create mean sea-ice snapshot for `$mean_variable`" exception = err
    end

    base_variable, base_title = DEFAULT_SEAICE_OVERLAY_BASE_FIELD
    overlay_variable, overlay_title = DEFAULT_SEAICE_OVERLAY_FIELD
    overlay_file = joinpath(
        plot_path,
        output_prefix * base_variable * "_" * overlay_variable * ".png",
    )
    try
        result = make_surface_field_overlay_plot(
            surface_path,
            base_variable,
            overlay_variable,
            overlay_file;
            base_title,
            overlay_title,
            dest,
        )
        !isnothing(result) && push!(output_files, result)
    catch err
        @warn "Failed to create sea-ice overlay plot for `$overlay_variable` on `$base_variable`" exception =
            err
    end

    for (variable_name, title) in surface_fields
        output_file = joinpath(plot_path, output_prefix * variable_name * "." * format)
        try
            result = make_ocean_field_movie(
                surface_path,
                variable_name,
                output_file;
                title,
                dest,
                framerate,
                format,
            )
            !isnothing(result) && push!(output_files, result)
        catch err
            @warn "Failed to create movie for sea-ice field `$variable_name`" exception = err
        end
    end

    if isempty(output_files)
        @info "No sea-ice diagnostic movies created in `$output_path`"
    else
        @info "Created $(length(output_files)) sea-ice diagnostic output file(s) in `$plot_path`"
    end

    return output_files
end
