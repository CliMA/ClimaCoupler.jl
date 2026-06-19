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

const DEFAULT_OCEAN_FIELDS_MOVIE_FIELDS = [
    ("temperature", "Temperature (surface)"),
    ("salinity", "Salinity (surface)"),
    ("zonal_velocity", "Zonal velocity (surface)"),
    ("meridional_velocity", "Meridional velocity (surface)"),
]

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
    underlying_grid = _underlying_ocean_grid(grid)
    longitude = Array(OC.λnodes(underlying_grid, OC.Center(), OC.Center()))
    latitude = Array(OC.φnodes(underlying_grid, OC.Center(), OC.Center()))
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
    longitude = Array(OC.λnodes(underlying_grid, ℓx, ℓy))
    latitude = Array(OC.φnodes(underlying_grid, ℓx, ℓy))
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

    fts = FieldTimeSeries(collection_path, variable_name)
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
    error("Run `using GLMakie` before calling `view_ocean_field_globe`.")
end

"""
    _wait_for_glmakie_screen(GLM, screen; wait = true)

Block until the GLMakie window is closed. Required when running from a script;
otherwise the Julia process exits immediately and the window disappears.
"""
function _wait_for_glmakie_screen(GLM, screen; wait = true)
    wait || return screen
    @info "Interactive globe open — close the window to exit (Ctrl+C also works)"
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

    fts = FieldTimeSeries(collection_path, variable_name)
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
    _wait_for_glmakie_screen(GLM, screen; wait)
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
- `<filename_prefix>_surface[...].jld2` for 2D surface diagnostics
- `<filename_prefix>_fields[...].jld2` for 3D fields (surface level plotted)

Movies are saved to `plot_path` as `<output_prefix><variable_name>.<format>`.
Returns a vector of paths to the created movie files.
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
        @info "No ocean diagnostic movies created in `$output_path`"
    else
        @info "Created $(length(output_files)) ocean diagnostic movie(s) in `$plot_path`"
    end

    return output_files
end
