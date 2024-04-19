# this file follows the ClimaAtmos.Diagnostics and ci_plots interfaces
import CairoMakie
import CairoMakie.Makie
import ClimaAnalysis as CAN
using Poppler_jll: pdfunite
import Base.Filesystem

const LARGE_NUM = typemax(Int)
const LAST_SNAP = LARGE_NUM
const FIRST_SNAP = -LARGE_NUM
const BOTTOM_LVL = -LARGE_NUM
const TOP_LVL = LARGE_NUM

function Makie.get_tickvalues(yticks::Int, ymin, ymax)
    return range(max(ymin, 0), ymax, yticks)
end

YLINEARSCALE = Dict(:axis => CAN.Utils.kwargs(dim_on_y = true, yticks = 10, ytickformat = "{:.3e}"))

long_name(var) = var.attributes["long_name"]
short_name(var) = var.attributes["short_name"]

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

    # Default plotting function needs access to kwargs
    if isnothing(plot_fn)
        plot_fn = (grid_loc, var) -> CAN.Visualize.plot!(grid_loc, var, args...; kwargs...)
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
                normalized_path = lpad(path, LABEL_LENGTH + 1, " ")[(end - LABEL_LENGTH):end]

                CairoMakie.Label(fig[0, col], path)
            end
        end
        return fig
    end

    # Standardizes grid layout
    gridlayout() =
        map(1:MAX_PLOTS_PER_PAGE) do i
            row = mod(div(i - 1, MAX_NUM_COLS), MAX_NUM_ROWS) + 1
            col = mod(i - 1, MAX_NUM_COLS) + 1
            return fig[row, col] = CairoMakie.GridLayout()
        end

    fig = makefig()
    grid = gridlayout()
    page = 1
    grid_pos = 1

    for var in vars
        if grid_pos > MAX_PLOTS_PER_PAGE
            fig = makefig()
            grid = gridlayout()
            grid_pos = 1
        end

        plot_fn(grid[grid_pos], var)
        grid_pos += 1

        # Flush current page
        if grid_pos > min(MAX_PLOTS_PER_PAGE, vars_left_to_plot)
            file_path = joinpath(plot_path, "$(output_name)_$page.pdf")
            CairoMakie.resize_to_layout!(fig)
            CairoMakie.save(file_path, fig)
            push!(summary_files, file_path)
            vars_left_to_plot -= MAX_PLOTS_PER_PAGE
            page += 1
        end
    end

    # Save plots
    output_file = joinpath(plot_path, "$(output_name).pdf")

    pdfunite() do unite
        run(Cmd([unite, summary_files..., output_file]))
    end

    # Cleanup
    Filesystem.rm.(summary_files, force = true)
    return output_file
end

function map_comparison(func, simdirs, args)
    return vcat([[func(simdir, arg) for simdir in simdirs] for arg in args]...)
end

"""
    make_plots(
        ::Union{Val{:general_ci_plots}},
        output_paths::Vector{<:AbstractString},
        plot_path::AbstractString;
        reduction::String = "average",
    )
Create plots for the general CI diagnostics. The plots are saved to `plot_path`.
This is the default plotting function for the CI diagnostics and it can be extended
to include additional diagnostics.
The `reduction` keyword argument should be consistent with the reduction used to save the diagnostics.
"""
function make_plots(
    ::Union{Val{:general_ci_plots}},
    output_paths::Vector{<:AbstractString},
    plot_path::AbstractString;
    reduction::String = "average",
)
    simdirs = CAN.SimDir.(output_paths)

    # Default output diagnostics
    short_names_3D = ["mse", "lr", "ediff"]
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
