import Dates

"""
    _open_series(path, name)

Open the Oceananigans `FieldTimeSeries` named `name` from the JLD2 output at
`path` (on disk, not loaded into memory), returning `nothing` if it cannot be
opened (e.g. the file or variable does not exist).
"""
function _open_series(path, name)
    try
        return OC.FieldTimeSeries(path, name; backend = OC.OnDisk())
    catch
        return nothing
    end
end

"""
    _geo_field_plot!(fig, field, short_name, title; p_loc = (1, 1))

Plot a single 2-D Oceananigans surface `field` on a global map at grid position
`p_loc`, using the shared plot styling: the Robinson projection
([`Plotting.PROJECTION`](@ref)), coastlines, and the per-variable colormap
([`Plotting.colormap_for`](@ref)). The field's own longitude/latitude nodes are
used as the surface coordinates, so this works on curvilinear (e.g. tripolar)
grids.
"""
function _geo_field_plot!(fig, field, short_name, title; p_loc = (1, 1))
    # `nodes` resolves the longitude/latitude coordinates at the field's own
    # location (which varies: tracers are cell-centered, velocity components are
    # face-centered), and returns 2-D arrays for curvilinear (e.g. tripolar) grids.
    λ, φ, _ = OC.nodes(field)
    data = Array(OC.interior(field, :, :, 1))

    finite = filter(isfinite, data)
    isempty(finite) && return nothing

    # Broadcast the coordinates to 2-D arrays matching `data`. This makes the plot
    # work uniformly for both regular lat/lon grids (1-D nodes) and curvilinear
    # grids (2-D nodes), and avoids a Makie `surface!` indexing bug that is
    # triggered on projected axes when the coordinate vectors differ in length.
    nx, ny = size(data)
    λ2d = ndims(λ) == 2 ? Array(λ) : [Array(λ)[i] for i in 1:nx, _ in 1:ny]
    φ2d = ndims(φ) == 2 ? Array(φ) : [Array(φ)[j] for _ in 1:nx, j in 1:ny]

    colormap = Plotting.colormap_for(short_name, data)
    plot_kwargs = Dict{Symbol, Any}(:colormap => colormap, :shading => Makie.NoShading)
    if Plotting.is_divergent(short_name, data)
        crange = Plotting.symmetric_colorrange(data)
        isnothing(crange) || (plot_kwargs[:colorrange] = crange)
    end

    row, col = p_loc
    ax = GeoMakie.GeoAxis(
        fig[row, 2col - 1];
        dest = Plotting.PROJECTION,
        title,
        xgridvisible = false,
        ygridvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        xticksvisible = false,
        yticksvisible = false,
    )
    hm = Makie.surface!(ax, λ2d, φ2d, data; plot_kwargs...)
    # Fill the continents in gray (these are ocean/sea-ice fields, so land carries
    # no data), then draw the coastlines on top.
    Makie.poly!(ax, GeoMakie.land(); color = :gray70)
    Makie.lines!(ax, GeoMakie.coastlines(); color = :black)
    # Narrow, label-free colorbar (the variable name and units are in the title).
    Makie.Colorbar(fig[row, 2col], hm; label = "", width = Plotting.COLORBAR_WIDTH)
    return nothing
end

"""
    _plot_field_series(specs, plot_path, title_prefix; plot_diagnostics = :all, times, dates)

Plot the Oceananigans surface diagnostics described by `specs` on global maps,
one figure per time step, to `plot_path`.

`specs` is a vector of `(short_name, title, getter)` tuples, where `getter(n)`
returns the field to plot for time index `n` (this indirection lets us build
derived fields such as surface speed from two velocity components). `times` are
the series times and `dates` their calendar dates (used to name the files).

`plot_diagnostics` selects which time steps are plotted:
- `:all` (default): one figure per time step, `<title_prefix>_summary_<date>.png`.
- `:last`: only the final time step, `<title_prefix>_summary.png`.
"""
function _plot_field_series(
    specs,
    plot_path,
    title_prefix;
    display_name = title_prefix,
    plot_diagnostics = :all,
    times = nothing,
    dates = nothing,
)
    isempty(specs) && return nothing

    (isnothing(times) || isempty(times)) && return nothing
    n_times = length(times)

    indices =
        plot_diagnostics == :last ? (n_times:n_times) :
        plot_diagnostics == :all ? (1:n_times) :
        error("`plot_diagnostics` must be `:all` or `:last`, got `$(plot_diagnostics)`")

    n = length(specs)
    ncols = ceil(Int, sqrt(n))
    nrows = ceil(Int, n / ncols)

    for i in indices
        date = isnothing(dates) ? nothing : dates[i]
        fig = CairoMakie.Figure(size = (700 * ncols, 450 * nrows))
        for (k, (short_name, title, getter)) in enumerate(specs)
            field = getter(i)
            isnothing(field) && continue
            row = div(k - 1, ncols) + 1
            col = mod(k - 1, ncols) + 1
            _geo_field_plot!(fig, field, short_name, title; p_loc = (row, col))
        end
        date_str = isnothing(date) ? "" : " at $(_date_str(date; sep = "-"))"
        label = "$(display_name) diagnostics$(date_str)"
        CairoMakie.Label(fig[0, :], label, fontsize = 20, font = :bold)

        suffix =
            plot_diagnostics == :last ? "" :
            "_" * (isnothing(date) ? string(i) : _date_str(date))
        output_file = joinpath(plot_path, "$(title_prefix)_summary$(suffix).png")
        CairoMakie.save(output_file, fig)
    end
    return nothing
end

"""
    _date_str(date; sep = "_")

Format a `FieldTimeSeries` time (a `Dates.DateTime`, or a number of seconds) as a
`yyyy<sep>mm<sep>dd` date, falling back to the integer seconds. The default
underscore separator is filename-friendly; pass `sep = "-"` for titles.
"""
_date_str(date::Dates.DateTime; sep = "_") = Dates.format(date, "yyyy$(sep)mm$(sep)dd")
_date_str(date; sep = "_") = string(round(Int, date)) * "s"

"""
    _series_dates(fts)

Return the calendar dates of the `FieldTimeSeries` `fts` if its times are
`DateTime`s, otherwise return its raw `times` (seconds).
"""
_series_dates(fts) = fts.times

"""
    Plotting.make_ocean_diagnostics_plots(output_path, plot_path; output_prefix = "", plot_diagnostics = :all)

Plot the Oceananigans ocean surface diagnostics found under `output_path` (the
`ocean_surface` JLD2 series) to `plot_path`, one global-map figure per time step.

The fields plotted, if available, are surface current speed (from the zonal and
meridional surface velocities), sea surface temperature, sea surface salinity,
and sea surface height.

`plot_diagnostics` (`:all` by default) selects whether every time step or only
the last is plotted (see [`Plotting.make_diagnostics_plots`](@ref)).
"""
function Plotting.make_ocean_diagnostics_plots(
    output_path::AbstractString,
    plot_path::AbstractString;
    output_prefix = "",
    plot_diagnostics = :all,
)
    path = joinpath(output_path, "ocean_surface")
    u = _open_series(path, "sea_surface_zonal_velocity")
    v = _open_series(path, "sea_surface_meridional_velocity")
    sst = _open_series(path, "sea_surface_temperature")
    sss = _open_series(path, "sea_surface_salinity")
    ssh = _open_series(path, "sea_surface_height")

    # Use any available series to define the time axis
    ref = something(sst, sss, ssh, u, v, Some(nothing))
    isnothing(ref) && return nothing

    specs = Tuple{String, String, Function}[]
    if !isnothing(u) && !isnothing(v)
        push!(
            specs,
            ("surface_speed", "Sea surface speed (m s⁻¹)", n -> _speed(u[n], v[n])),
        )
    end
    isnothing(sst) ||
        push!(specs, ("sea_surface_temperature", "SST (°C)", n -> sst[n]))
    isnothing(sss) || push!(specs, ("sea_surface_salinity", "SSS (PSU)", n -> sss[n]))
    isnothing(ssh) || push!(specs, ("sea_surface_height", "SSH (m)", n -> ssh[n]))

    _plot_field_series(
        specs,
        plot_path,
        output_prefix * "ocean";
        display_name = "Ocean",
        plot_diagnostics,
        times = ref.times,
        dates = _series_dates(ref),
    )
    return nothing
end

"""
    Plotting.make_seaice_diagnostics_plots(output_path, plot_path; output_prefix = "", plot_diagnostics = :all)

Plot the Oceananigans sea ice surface diagnostics found under `output_path` (the
`seaice_surface` JLD2 series) to `plot_path`, one global-map figure per time step.

The fields plotted, if available, are ice drift speed (from the zonal and
meridional ice velocities), ice concentration, ice thickness, and ice top-surface
temperature.

`plot_diagnostics` (`:all` by default) selects whether every time step or only
the last is plotted.
"""
function Plotting.make_seaice_diagnostics_plots(
    output_path::AbstractString,
    plot_path::AbstractString;
    output_prefix = "",
    plot_diagnostics = :all,
)
    path = joinpath(output_path, "seaice_surface")
    ui = _open_series(path, "ice_zonal_velocity")
    vi = _open_series(path, "ice_meridional_velocity")
    conc = _open_series(path, "ice_concentration")
    thickness = _open_series(path, "ice_thickness")
    temp = _open_series(path, "ice_top_temperature")

    ref = something(conc, thickness, temp, ui, vi, Some(nothing))
    isnothing(ref) && return nothing

    specs = Tuple{String, String, Function}[]
    if !isnothing(ui) && !isnothing(vi)
        push!(
            specs,
            ("surface_speed", "Sea-ice speed (m s⁻¹)", n -> _speed(ui[n], vi[n])),
        )
    end
    isnothing(conc) ||
        push!(specs, ("ice_concentration", "Ice concentration", n -> conc[n]))
    isnothing(thickness) ||
        push!(specs, ("ice_thickness", "Ice thickness (m)", n -> thickness[n]))
    isnothing(temp) ||
        push!(specs, ("surface_temperature", "Ice top temperature (°C)", n -> temp[n]))

    _plot_field_series(
        specs,
        plot_path,
        output_prefix * "seaice";
        display_name = "Sea ice",
        plot_diagnostics,
        times = ref.times,
        dates = _series_dates(ref),
    )
    return nothing
end

"""
    _speed(u, v)

Return the magnitude field `sqrt(u^2 + v^2)` of two Oceananigans velocity
component fields, computed into a concrete `Field`.
"""
function _speed(u, v)
    speed = OC.Field(sqrt(u^2 + v^2))
    OC.compute!(speed)
    return speed
end
