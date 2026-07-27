import CairoMakie
import ClimaAnalysis as CAN
import ClimaComms
import ClimaCore as CC
import Dates

# `OrderedDict` (from OrderedCollections) is used to build `OutputVar` dimensions
# with a guaranteed order. We reach it through ClimaAnalysis rather than adding a
# direct dependency on OrderedCollections.
const OrderedDict = CAN.Var.OrderedDict

import ClimaCoupler: Interfacer, Plotting

"""
    _field_to_outputvar(field, short_name; hresolution = 180)

Interpolate a horizontal `ClimaCore` `field` (defined on a global spectral
element space) onto a regular lat/lon grid and wrap the result in a
`ClimaAnalysis.OutputVar` so that it can be plotted with the same styled global
plotting machinery used by the diagnostics postprocessing.

Returns `nothing` if the field is not defined on a global (sphere) horizontal
space, or on any non-root process (interpolation only returns data on root).

`hresolution` sets the number of points along each of the latitude and
longitude dimensions.
"""
function _field_to_outputvar(field, short_name; hresolution = 180)
    cpu_field = CC.to_cpu(field)
    space = axes(cpu_field)

    # Only global (sphere) spectral element spaces can be plotted on a globe.
    space isa CC.Spaces.SpectralElementSpace2D || return nothing
    topology = CC.Spaces.topology(space)
    CC.Meshes.domain(topology.mesh) isa CC.Domains.SphereDomain || return nothing

    # `default_target_hcoords_as_vectors` returns (latitudes, longitudes) and
    # `interpolate` returns an array indexed as [latitude, longitude].
    lats, lons =
        CC.Remapping.default_target_hcoords_as_vectors(space; hresolution)
    remapper = CC.Remapping.Remapper(space; target_hcoords = CC.Remapping.default_target_hcoords(space; hresolution))
    data = CC.Remapping.interpolate(remapper, cpu_field)
    # `interpolate` returns `nothing` off the root process
    isnothing(data) && return nothing

    # `ClimaAnalysis`' global plotting expects the data indexed as
    # [longitude, latitude] (it calls `surface!(ax, lon, lat, data)`), so order
    # the dimensions as (longitude, latitude) and transpose the [lat, lon] data
    # from `interpolate` accordingly. Getting this wrong transposes every plot.
    dims = OrderedDict("longitude" => Array(lons), "latitude" => Array(lats))
    dim_attribs = OrderedDict(
        "longitude" => Dict("units" => "degrees_east"),
        "latitude" => Dict("units" => "degrees_north"),
    )
    attribs = Dict{String, Any}(
        "short_name" => string(short_name),
        "long_name" => _prettify_name(short_name),
        "units" => "",
    )
    return CAN.OutputVar(attribs, dims, dim_attribs, permutedims(Array(data)))
end

# Human-readable titles for short names that would otherwise prettify to a
# cryptic single letter (the signed surface velocity components).
const _PRETTY_NAME_OVERRIDES =
    Dict("u" => "Zonal velocity", "v" => "Meridional velocity", "w" => "Vertical velocity")

"""
    _prettify_name(short_name)

Turn a variable's short name (e.g. `:surface_temperature`) into a human-readable
title (e.g. `"Surface temperature"`) for use as a plot title. A few short names
(the velocity components `u`/`v`/`w`) are mapped to spelled-out titles via
`_PRETTY_NAME_OVERRIDES`.
"""
function _prettify_name(short_name)
    key = lowercase(string(short_name))
    haskey(_PRETTY_NAME_OVERRIDES, key) && return _PRETTY_NAME_OVERRIDES[key]
    words = split(replace(string(short_name), "_" => " "))
    isempty(words) && return string(short_name)
    return uppercasefirst(join(words, " "))
end

"""
    Plotting.snapshot_plot(place, var::CAN.OutputVar; p_loc = (1, 1))

Plot a single 2-D lat/lon `var` on a global map at grid position `p_loc`, using
the shared plot styling: a per-variable colormap ([`Plotting.colormap_for`](@ref)),
the Robinson projection ([`Plotting.PROJECTION`](@ref)), and coastlines. This is
the common rendering path shared by the snapshot callback and the diagnostics
postprocessing, so both produce consistent-looking plots.
"""
function Plotting.snapshot_plot(place, var::CAN.OutputVar; p_loc = (1, 1), mask_land = false)
    var = _ensure_square_latlon(var)
    title = _var_title(var)
    more_kwargs = Plotting.geo_plot_kwargs(CAN.short_name(var), var.data; title)
    # `geo_plot_kwargs` lives in the (Makie-free) base package, so the colorbar
    # height (a `Makie.Relative`) has to be filled in here.
    more_kwargs[:cb][:height] = CairoMakie.Relative(Plotting.COLORBAR_HEIGHT_FRACTION)
    # For ocean/sea-ice fields, blank out and gray-fill the continents (they carry
    # no data there). `landmask()` both masks the land data and draws the polygons,
    # whose color we set via the `:mask` kwargs.
    mask = mask_land ? CAN.Visualize.landmask() : nothing
    mask_land && (more_kwargs[:mask] = Dict(:color => :gray70))
    CAN.Visualize.heatmap2D_on_globe!(place, var; p_loc, mask, more_kwargs)
    return nothing
end

"""
    _var_title(var)

Build a panel title for `var` from its long name (falling back to its short
name) and its units, e.g. `"Sea surface temperature [K]"`. The units are omitted
if unknown/empty.
"""
function _var_title(var)
    name = get(var.attributes, "long_name", "")
    isempty(name) && (name = string(CAN.short_name(var)))
    units = CAN.units(var)
    return isempty(units) ? name : "$name [$units]"
end

"""
    _ensure_square_latlon(var)

Return `var` resampled onto a square (equal number of latitude and longitude
points) grid if it is a 2-D lat/lon variable whose latitude and longitude
dimensions have different lengths, otherwise return `var` unchanged.

This works around an upstream limitation: the projected global heatmap
(`Makie.surface!` on a GeoMakie `GeoAxis`) mis-indexes and errors when the two
horizontal dimensions have different lengths. Snapshot variables are already
square (they are remapped onto an `n × n` grid), but diagnostics read from
NetCDF are typically non-square (e.g. 90 × 180), so we square them here. The
target resolution is the larger of the two dimensions, so no detail is lost.
"""
function _ensure_square_latlon(var)
    (CAN.has_longitude(var) && CAN.has_latitude(var)) || return var
    lats = CAN.latitudes(var)
    lons = CAN.longitudes(var)
    length(lats) == length(lons) && return var

    n = max(length(lats), length(lons))
    new_lats = collect(range(first(lats), last(lats), length = n))
    new_lons = collect(range(first(lons), last(lons), length = n))
    return CAN.resampled_as(var; latitude = new_lats, longitude = new_lons)
end

"""
    Plotting.plot_snapshots(cs::Interfacer.CoupledSimulation, date)

Create instantaneous ("snapshot") plots of each component model's key surface
fields, plus the coupler exchange fields, and save them (one figure each) to
`artifacts/<component>/snapshot_<date>.png` and
`artifacts/coupler/snapshot_<date>.png`.

Snapshots are cheap, read directly from the live state via `Interfacer.get_field`
(so they are produced even for runs that later crash), and are meant to be
generated periodically during a run by a callback (see the `plot_interval`
configuration flag).

The fields plotted for each component are given by
[`Plotting.snapshot_plot_fields`](@ref), which can be specialized per component.
All plotting is performed on the root process only.
"""
function Plotting.plot_snapshots(cs::Interfacer.CoupledSimulation, date)
    ClimaComms.iamroot(ClimaComms.context(cs)) || return nothing
    boundary_space = Interfacer.boundary_space(cs)
    for (sim_name, sim) in pairs(cs.model_sims)
        Plotting.plot_snapshots(
            sim,
            sim_name,
            boundary_space,
            cs.dir_paths.artifacts_dir,
            date,
        )
    end
    # Also plot the coupler exchange fields (fluxes, area fractions, etc.)
    Plotting.plot_snapshots(cs.fields, cs.dir_paths.artifacts_dir, date)
    return nothing
end

"""
    Plotting.plot_snapshots(sim, sim_name, boundary_space, artifacts_dir, date)

Create a single figure of snapshot plots for one component `sim` (identified by
its coupler key `sim_name`, e.g. `:atmos_sim`) and save it to
`artifacts_dir/<sim_name>/snapshot_<date>.png`. One panel is drawn per field
returned by [`Plotting.snapshot_plot_fields`](@ref).

Each field is remapped onto the coupler `boundary_space` before plotting, so
this works uniformly for all component models regardless of their native grid
(e.g. the atmosphere's spectral element grid or the ocean's finite volume grid).
"""
function Plotting.plot_snapshots(
    sim::Interfacer.AbstractComponentSimulation,
    sim_name,
    boundary_space,
    artifacts_dir,
    date,
)
    field_names = Plotting.snapshot_plot_fields(sim)
    isempty(field_names) && return nothing

    vars = [
        (field_name, _snapshot_var(sim, boundary_space, field_name)) for
        field_name in field_names
    ]
    title = _snapshot_figure_title(sim_name, date)
    return _save_snapshot_figure(vars, title, artifacts_dir, sim_name, date)
end

"""
    _snapshot_file_date(date)

Format `date` as a filename-friendly ISO-8601 timestamp, e.g.
`"2010-01-01T12:23:01"`, for use in `snapshot_<date>.png`. Unlike the date-only
figure title, this keeps the time of day so that sub-daily snapshots do not
collide.
"""
_snapshot_file_date(date) = Dates.format(date, "yyyy-mm-ddTHH:MM:SS")

"""
    _snapshot_figure_title(sim_name, date)

Build the overall figure title for a snapshot, e.g. `"Ocean snapshots at
2010-01-01"`. `sim_name` is the coupler key (`:ocean_sim`, `:atmos_sim`, ...) or
`"coupler"`.
"""
function _snapshot_figure_title(sim_name, date)
    labels = Dict(
        "atmos_sim" => "Atmosphere",
        "ocean_sim" => "Ocean",
        "land_sim" => "Land",
        "ice_sim" => "Sea ice",
        "coupler" => "Coupler",
    )
    label = get(labels, string(sim_name), string(sim_name))
    return "$(label) snapshots at $(Dates.format(date, "yyyy-mm-dd"))"
end

"""
    Plotting.plot_snapshots(cs_fields::CC.Fields.Field, artifacts_dir, date)

Create a single figure of snapshot plots of the coupler exchange fields
(`cs.fields`: fluxes, area fractions, boundary-layer quantities, etc.) and save
it to `artifacts_dir/coupler/snapshot_<date>.png`. The coupler fields are already
defined on the boundary space, so no remapping is needed.
"""
function Plotting.plot_snapshots(cs_fields::CC.Fields.Field, artifacts_dir, date)
    field_names = propertynames(cs_fields)
    vars = [
        (field_name, _field_to_outputvar(getproperty(cs_fields, field_name), field_name))
        for field_name in field_names
    ]
    title = _snapshot_figure_title("coupler", date)
    return _save_snapshot_figure(vars, title, artifacts_dir, "coupler", date)
end

"""
    _save_snapshot_figure(vars, title, artifacts_dir, sim_name, date)

Given `vars`, a vector of `(field_name, var_or_nothing)` pairs, lay out one map
panel per plottable variable in a roughly square grid and save the figure to
`artifacts_dir/<sim_name>/snapshot_<date>.png`.

Fields that cannot be plotted on a globe (e.g. in single-column mode, where
`var` is `nothing`) are skipped. Fields that are entirely NaN/Inf are skipped
with a warning (a NaN can indicate a misbehaving run). Returns the output path,
or `nothing` if there was nothing to plot.
"""
function _save_snapshot_figure(vars, title, artifacts_dir, sim_name, date)
    plottable = CAN.OutputVar[]
    for (field_name, var) in vars
        isnothing(var) && continue
        if all(x -> isnan(x) || isinf(x), var.data)
            @warn "Snapshot field $(field_name) of $(sim_name) is entirely NaN/Inf - skipping"
            continue
        end
        any(x -> isnan(x) || isinf(x), var.data) &&
            @warn "Snapshot field $(field_name) of $(sim_name) contains NaN/Inf values"
        push!(plottable, var)
    end
    isempty(plottable) && return nothing

    # Lay out panels in a roughly square grid. Each panel uses two columns so
    # that `heatmap2D_on_globe!` can place its colorbar in the adjacent column.
    # Ocean and sea-ice fields carry no data over land, so gray out the continents.
    mask_land = string(sim_name) in ("ocean_sim", "ice_sim")

    elapsed = @elapsed begin
        n = length(plottable)
        ncols = ceil(Int, sqrt(n))
        nrows = ceil(Int, n / ncols)
        fig = CairoMakie.Figure(size = (600 * ncols, 400 * nrows))
        for (i, var) in enumerate(plottable)
            row = div(i - 1, ncols) + 1
            col = mod(i - 1, ncols) + 1
            Plotting.snapshot_plot(fig, var; p_loc = (row, 2 * col - 1), mask_land)
        end
        CairoMakie.Label(fig[0, :], title, fontsize = 20, font = :bold)

        dir = Plotting.component_artifacts_dir(artifacts_dir, sim_name)
        output_file = joinpath(dir, "snapshot_$(_snapshot_file_date(date)).png")
        CairoMakie.save(output_file, fig)
    end
    @info "Saved $(output_file) in $(round(elapsed, digits = 2)) s"
    return output_file
end

"""
    _snapshot_var(sim, boundary_space, field_name)

Return a lat/lon `CAN.OutputVar` for the snapshot field `field_name` of `sim`,
remapped onto `boundary_space`, or `nothing` if it cannot be plotted.

Special-cases the pseudo-field `:surface_speed`, which is computed as
`hypot(u, v)` from the component's near-surface horizontal velocity. The
velocity components are read from `:u_int`/`:v_int` (atmosphere) or `:u`/`:v`
(ocean), whichever the component provides.
"""
function _snapshot_var(sim, boundary_space, field_name)
    if field_name === :surface_speed
        speed = _surface_speed(sim, boundary_space)
        isnothing(speed) && return nothing
        return _field_to_outputvar(speed, :surface_speed)
    end
    field = _try_get_field(boundary_space, sim, field_name)
    isnothing(field) && return nothing
    return _field_to_outputvar(field, field_name)
end

"""
    _surface_speed(sim, boundary_space)

Return the near-surface horizontal wind/current speed of `sim` on
`boundary_space`, or `nothing` if the component does not expose horizontal
velocity components. Tries `:u_int`/`:v_int` (used by the atmosphere) first,
then `:u`/`:v` (used by the ocean).
"""
function _surface_speed(sim, boundary_space)
    for (u_name, v_name) in ((:u_int, :v_int), (:u, :v))
        u = _try_get_field(boundary_space, sim, u_name)
        v = _try_get_field(boundary_space, sim, v_name)
        (isnothing(u) || isnothing(v)) && continue
        return @. hypot(u, v)
    end
    return nothing
end

"""
    _try_get_field(boundary_space, sim, field_name)

Like `Interfacer.get_field(boundary_space, sim, Val(field_name))`, but return
`nothing` instead of erroring if the field is not defined for `sim`.
"""
function _try_get_field(boundary_space, sim, field_name)
    try
        return Interfacer.get_field(boundary_space, sim, Val(field_name))
    catch e
        @warn "Snapshot: could not get field $(field_name) for $(nameof(typeof(sim)))" exception =
            (e, catch_backtrace())
        return nothing
    end
end

"""
    Plotting.snapshot_plot_fields(sim::Interfacer.AbstractSurfaceSimulation)

Return the default set of fields to include in instantaneous snapshot plots for
a surface model. This is intended to be a small, cheap, curated set of the most
useful fields, and can be specialized per component model.
"""
Plotting.snapshot_plot_fields(sim::Interfacer.AbstractSurfaceSimulation) =
    (:area_fraction, :surface_temperature)
