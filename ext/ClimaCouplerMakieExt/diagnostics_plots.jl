import CairoMakie
import CairoMakie.Makie
import ClimaAnalysis as CAN
import Dates
using Poppler_jll: pdfunite

const LARGE_NUM = typemax(Int)
const LAST_SNAP = LARGE_NUM
const FIRST_SNAP = -LARGE_NUM
const BOTTOM_LVL = -LARGE_NUM
const TOP_LVL = LARGE_NUM

function Makie.get_tickvalues(yticks::Int, ymin, ymax)
    return range(ymin, ymax, yticks)
end

YLINEARSCALE =
    Dict(:axis => CAN.Utils.kwargs(dim_on_y = true, yticks = 10, ytickformat = "{:.3e}"))

long_name(var) = var.attributes["long_name"]
short_name(var) = var.attributes["short_name"]

"""
    _styled_plot!(grid_loc, var; kwargs...)

Plot a single diagnostic `var` at `grid_loc` using the shared plot styling.

Any keyword arguments (e.g. `time = ...`) are treated as slices and applied to
`var` first (by nearest value), so that after slicing the variable is reduced to
the dimensions being plotted. This mirrors how `ClimaAnalysis`'s `plot!` consumes
its slicing keywords, but lets us inspect the *sliced* data to choose styling.

For 2-D lat/lon variables this applies the per-variable colormap, Robinson
projection, and coastlines (see [`Plotting.geo_plot_kwargs`](@ref)), so that the
diagnostics plots match the instantaneous snapshot plots. Other variables (e.g.
vertical profiles or zonally-averaged fields) fall back to the default
`ClimaAnalysis` plotting.
"""
function _styled_plot!(grid_loc, var; kwargs...)
    # Apply the slicing keywords (e.g. `time`) up front so we style and plot the
    # already-sliced variable. Only slice by dimensions the variable actually has
    # (`slice` slices by nearest value by default, so `time = LAST_SNAP` selects
    # the last time step).
    for (dim_name, val) in kwargs
        haskey(var.dims, string(dim_name)) || continue
        var = CAN.slice(var; NamedTuple{(dim_name,)}((val,))...)
    end

    is_latlon_2d =
        length(var.dims) == 2 &&
        (CAN.has_longitude(var) && CAN.has_latitude(var))
    if is_latlon_2d
        # Reuse the shared global-map renderer (Robinson projection + coastlines +
        # per-variable colormap) so diagnostics match the snapshot plots.
        Plotting.snapshot_plot(grid_loc, var)
    else
        CAN.Visualize.plot!(grid_loc, var)
    end
    return nothing
end

"""
    make_diagnostics_plots(
        output_path::AbstractString,
        plot_path::AbstractString;
        output_prefix = "",
        plot_diagnostics = :all,
    )
Create plots for diagnostics. The plots are saved to `plot_path` (typically a
per-component subdirectory of `artifacts`, e.g. `artifacts/atmos_sim`).
This function will plot all variables that have been saved in `output_path`.

`plot_diagnostics` controls which time steps are plotted:
- `:all` (default): plot every saved time step, producing one summary file per
  time step (named with the date, e.g. `summary_2D_2010_01_01.pdf`).
- `:last`: plot only the last saved time step (i.e. the final averaging window),
  producing a single summary file (e.g. `summary_2D.pdf`).

When plotting diagnostics, diagnostics that are daily averages in z coordinates
(if available) are prioritized first.

Two-dimensional lat/lon variables are plotted on a global map using the shared
plot styling (per-variable colormap, Robinson projection, and coastlines; see
[`Plotting.geo_plot_kwargs`](@ref)), so that these diagnostics match the
instantaneous snapshot plots.

For column / box-grid output (no lat/lon dimensions), degenerate horizontal
dimensions are sliced out and 3-D variables are plotted as vertical profiles
(value vs z) instead of being zonally averaged.
"""
function Plotting.make_diagnostics_plots(
    output_path::AbstractString,
    plot_path::AbstractString;
    output_prefix = "",
    plot_diagnostics = :all,
)
    simdir = CAN.SimDir(output_path)
    short_names = CAN.available_vars(simdir)

    # Return if there are no variables to plot
    isempty(short_names) && return

    # Create a CAN.OutputVar for each input field
    vars = Array{CAN.OutputVar}(undef, length(short_names))
    for (i, short_name) in enumerate(short_names)
        # Use "average" if available, otherwise use the first reduction
        reductions = CAN.available_reductions(simdir; short_name)
        reduction = "average" in reductions ? "average" : first(reductions)
        periods = CAN.available_periods(simdir; short_name, reduction)
        period = "1d" in periods ? "1d" : first(periods)
        coord_types = CAN.available_coord_types(simdir; short_name, reduction, period)
        coord_type = nothing in coord_types ? nothing : first(coord_types)
        vars[i] = get(simdir; short_name, reduction, period, coord_type)
    end

    is_3d = var -> CAN.has_altitude(var) || CAN.has_pressure(var)
    has_latlon = any(v -> CAN.has_longitude(v) || CAN.has_latitude(v), vars)

    if has_latlon
        # Global mode: zonally-averaged 3-D fields on lat-z, 2-D fields on lon-lat
        vars_3D = map(var_3D -> CAN.average_lon(var_3D), filter(is_3d, vars))
        vars_2D = filter(var -> !is_3d(var), vars)

        _plot_var_group(
            output_path,
            plot_path,
            vars_3D,
            output_prefix * "summary_3D",
            plot_diagnostics;
            more_kwargs = YLINEARSCALE,
        )
        _plot_var_group(
            output_path,
            plot_path,
            vars_2D,
            output_prefix * "summary_2D",
            plot_diagnostics;
            plot_fn = _styled_plot!,
        )
    else
        # Column / box-grid mode: slice out degenerate horizontal dims
        vars_profile = _slice_to_column.(filter(is_3d, vars))
        vars_surface = _slice_to_column.(filter(var -> !is_3d(var), vars))

        _plot_var_group(
            output_path,
            plot_path,
            vars_profile,
            output_prefix * "summary_profiles",
            plot_diagnostics;
            more_kwargs = YLINEARSCALE,
        )
        _plot_var_group(
            output_path,
            plot_path,
            vars_surface,
            output_prefix * "summary_surface",
            plot_diagnostics,
        )
    end
end

"""
    _plot_var_group(output_path, plot_path, vars, base_name, plot_diagnostics; kwargs...)

Plot a group of diagnostic `vars` (e.g. all 2-D lat/lon fields) to `plot_path`,
dispatching on `plot_diagnostics`:
- `:last`: one summary file `<base_name>.pdf` at the last time step.
- `:all`: one summary file per saved time step, `<base_name>_<date>.pdf`.

For `:all`, the set of time steps is taken from the union of the `vars`' time
dimensions; variables without a time dimension are plotted once (at their single
value) in every time step's file. Extra `kwargs` (e.g. `plot_fn`, `more_kwargs`)
are forwarded to [`make_plots_generic`](@ref).
"""
function _plot_var_group(
    output_path,
    plot_path,
    vars,
    base_name,
    plot_diagnostics;
    kwargs...,
)
    isempty(vars) && return nothing

    if plot_diagnostics == :last
        make_plots_generic(
            output_path,
            plot_path,
            vars,
            time = LAST_SNAP,
            output_name = base_name;
            kwargs...,
        )
        return nothing
    end

    plot_diagnostics == :all || error(
        "`plot_diagnostics` must be `:all` or `:last`, got `$(plot_diagnostics)`",
    )

    # Collect the union of all time steps across the variables that have a time
    # dimension, then produce one summary file per time step.
    times = Float64[]
    for var in vars
        CAN.has_time(var) && union!(times, CAN.times(var))
    end
    sort!(times)

    # If no variable has a time dimension, there is a single (timeless) plot.
    if isempty(times)
        make_plots_generic(
            output_path,
            plot_path,
            vars,
            output_name = base_name;
            kwargs...,
        )
        return nothing
    end

    for t in times
        make_plots_generic(
            output_path,
            plot_path,
            vars,
            time = t,
            output_name = "$(base_name)_$(_time_suffix(vars, t))";
            kwargs...,
        )
    end
    return nothing
end

"""
    _time_suffix(vars, t)

Return a filename-friendly suffix for the time step `t` (seconds): the date
`yyyy_mm_dd` if any of `vars` carries a `start_date` (so a calendar date can be
computed), otherwise the integer number of seconds.
"""
function _time_suffix(vars, t)
    for var in vars
        if CAN.has_time(var) && haskey(var.attributes, "start_date")
            date = Dates.DateTime(var.attributes["start_date"]) + Dates.Second(round(Int, t))
            return Dates.format(date, "yyyy_mm_dd")
        end
    end
    return string(round(Int, t)) * "s"
end



"""
    _slice_to_column(var)

Remove horizontal dimensions from a `CAN.OutputVar` by slicing at the first
value. Used for column / box-grid diagnostics saved by the atmosphere in
column mode to get 1D vertical profiles.
"""
function _slice_to_column(var)
    CAN.has_longitude(var) && (var = CAN.slice(var, by = CAN.Index(), longitude = 1))
    CAN.has_latitude(var) && (var = CAN.slice(var, by = CAN.Index(), latitude = 1))

    dim_names = collect(keys(var.dims))
    "x" in dim_names && (var = CAN.slice(var, by = CAN.Index(), x = 1))
    "y" in dim_names && (var = CAN.slice(var, by = CAN.Index(), y = 1))
    return var
end

"""
    make_plots_generic(
        file_path::Union{<:AbstractString, Vector{<:AbstractString}},
        plot_path,
        vars,
        args...;
        plot_fn = nothing,
        output_name = "summary",
        summary_files = String[],
        MAX_NUM_COLS = 1,
        MAX_NUM_ROWS = min(4, length(vars)),
        kwargs...,
    )
Create plots for each variable in `vars` and save them to `plot_path`. The number of plots per
page is determined by `MAX_NUM_COLS` and `MAX_NUM_ROWS`. The `plot_fn` function is used to create the
plots. If `plot_fn` is not provided, a default plotting function is used. The default plotting function
is determined by the keyword arguments `kwargs`.
"""
function make_plots_generic(
    file_path::Union{<:AbstractString, Vector{<:AbstractString}},
    plot_path,
    vars,
    args...;
    plot_fn = nothing,
    output_name = "summary",
    summary_files = String[],
    MAX_NUM_COLS = 1,
    MAX_NUM_ROWS = min(4, length(vars)),
    kwargs...,
)
    # When file_path is a Vector with multiple elements, this means that this function is
    # being used to produce a comparison plot. In that case, we modify the output name, and
    # the number of columns (to match how many simulations we are comparing).
    is_comparison = file_path isa Vector
    #
    # However, we don't want to do this when the vector only contains one element.
    if is_comparison && length(file_path) == 1
        # Fallback to the "file_path isa String" case
        file_path = file_path[1]
        is_comparison = false
    end

    if is_comparison
        MAX_NUM_COLS = length(file_path)
        plot_path = file_path[1]
        output_name *= "_comparison"
    end

    # Both the default and any provided `plot_fn` need access to the slicing
    # `args`/`kwargs` (e.g. `time = ...`), so close over them here. A provided
    # `plot_fn` (e.g. `_styled_plot!`) receives them as keyword arguments.
    user_plot_fn = plot_fn
    if isnothing(user_plot_fn)
        plot_fn = (grid_loc, var) -> CAN.Visualize.plot!(grid_loc, var, args...; kwargs...)
    else
        plot_fn = (grid_loc, var) -> user_plot_fn(grid_loc, var; kwargs...)
    end

    MAX_PLOTS_PER_PAGE = MAX_NUM_ROWS * MAX_NUM_COLS
    vars_left_to_plot = length(vars)

    # Define fig, grid, and grid_pos, used below. (Needed for scope)
    function makefig()
        fig = CairoMakie.Figure(; size = (900, 300 * MAX_NUM_ROWS))
        if is_comparison
            for (col, path) in (file_path)
                # CairoMakie seems to use this Label to determine the width of the figure.
                # Here we normalize the length so that all the columns have the same width.
                LABEL_LENGTH = 40
                path = convert(Vector{Float64}, path)
                normalized_path =
                    lpad(path, LABEL_LENGTH + 1, " ")[(end - LABEL_LENGTH):end]

                CairoMakie.Label(fig[0, col], path)
            end
        end
        return fig
    end

    # Standardizes grid layout
    gridlayout() = map(1:MAX_PLOTS_PER_PAGE) do i
        row = mod(div(i - 1, MAX_NUM_COLS), MAX_NUM_ROWS) + 1
        col = mod(i - 1, MAX_NUM_COLS) + 1
        return fig[row, col] = CairoMakie.GridLayout()
    end

    fig = makefig()
    grid = gridlayout()
    page = 1
    grid_pos = 1

    for var in vars
        if all(isnan, var.data)
            @warn "$(short_name(var)) diagnostic is entirely NaN - skipping plot"
            vars_left_to_plot -= 1
            continue
        end
        if minimum(var.data) == maximum(var.data)
            @warn "$(short_name(var)) diagnostic is spatially constant - skipping plot"
            vars_left_to_plot -= 1
            continue
        end

        if grid_pos > MAX_PLOTS_PER_PAGE
            fig = makefig()
            grid = gridlayout()
            grid_pos = 1
        end

        plot_fn(grid[grid_pos], var)
        grid_pos += 1

        # Flush current page
        if grid_pos > min(MAX_PLOTS_PER_PAGE, vars_left_to_plot)
            # Save current page as a separate PDF file
            file_path = joinpath(plot_path, "$(output_name)_$(page).pdf")
            CairoMakie.resize_to_layout!(fig)
            CairoMakie.save(file_path, fig)
            push!(summary_files, file_path)
            vars_left_to_plot -= MAX_PLOTS_PER_PAGE
            page += 1
        end
    end

    # Return early if there are no plots to save
    isempty(summary_files) && return nothing

    # Save plots
    output_file = joinpath(plot_path, "$(output_name).pdf")

    pdfunite() do unite
        run(Cmd([unite, summary_files..., output_file]))
    end

    # Cleanup
    rm.(summary_files, force = true)
    return output_file
end

function map_comparison(func, simdirs, args)
    return vcat([[func(simdir, arg) for simdir in simdirs] for arg in args]...)
end
