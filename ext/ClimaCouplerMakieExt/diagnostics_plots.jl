import CairoMakie
import CairoMakie.Makie
import ClimaAnalysis as CAN
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
    plot_spectrum_with_line!(grid_loc, spectrum; exponent = -3.0)

Plots the given spectrum alongside a line that identifies a power law.
Assumes 1D spectrum (e.g. spectrum sliced at one level).
"""
function plot_spectrum_with_line!(grid_loc, spectrum; exponent = -3.0)
    ndims(spectrum.data) == 1 || error("Can only work with 1D spectrum")
    CAN.Visualize.plot!(grid_loc, spectrum)

    dim_name = spectrum.index2dim[begin]
    ax = CairoMakie.current_axis()

    # Ignore below wavenumber of 10
    spectrum_10 = CAN.window(spectrum, dim_name; left = log10(10))

    # Add reference line
    wavenumbers = spectrum_10.dims[dim_name]
    max_spectrum_10 = maximum(spectrum_10.data)
    wavenumber_at_max = wavenumbers[argmax(spectrum_10.data)]
    intercept = 1.2 * (max_spectrum_10 - exponent * wavenumber_at_max)
    reference_line(k) = exponent * k + intercept

    color = :orange
    CairoMakie.lines!(ax, wavenumbers, reference_line.(wavenumbers); color)
    CairoMakie.text!(
        ax,
        wavenumber_at_max,
        reference_line(wavenumber_at_max),
        text = "k^$exponent";
        color,
    )
    return nothing
end

"""
    make_atmos_spectra_plots(
        output_path::AbstractString,
        plot_path::AbstractString;
        output_prefix = "",
    )

Create spectra plots for atmosphere diagnostics (currently `ua`, `ta`, `hus`).

Reads variables from `output_path` via ClimaAnalysis, computes spectra using
`ClimaCouplerClimaCoreExt.compute_spectrum`, slices at `z = 1500`, and saves the result
to ``output_prefix * "spectra.pdf"`` under `plot_path`. Returns `nothing` (no plots) when
ClimaCoreSpectra is not loaded; load it (e.g. `using ClimaCoreSpectra`) to enable spectra plots.
"""
function Plotting.make_atmos_spectra_plots(
    output_path::AbstractString,
    plot_path::AbstractString;
    output_prefix = "",
)
    simdir = CAN.SimDir(output_path)
    short_names = CAN.available_vars(simdir)
    isempty(short_names) && return nothing

    short_names_spectra = ["ua", "ta", "hus"]
    vars_spectra = CAN.OutputVar[]
    # ClimaCouplerClimaCoreExt loads only when the user does `using ClimaCoreSpectra` (it is
    # not triggered when ClimaCoreSpectra is loaded transitively as a dep of ClimaCouplerMakieExt).
    ext = Base.get_extension(ClimaCoupler, :ClimaCouplerClimaCoreExt)
    isnothing(ext) && return nothing
    compute_spectrum_fn = ext.compute_spectrum

    for sn in short_names_spectra
        sn in short_names || continue
        reductions = CAN.available_reductions(simdir; short_name = sn)
        reduction = "average" in reductions ? "average" : first(reductions)
        periods = CAN.available_periods(simdir; short_name = sn, reduction)
        period = "1d" in periods ? "1d" : first(periods)

        var = get(simdir; short_name = sn, reduction, period)
        var = CAN.slice(var; time = LAST_SNAP)
        var = compute_spectrum_fn(var)
        var = CAN.slice(var; z = 1500)
        push!(vars_spectra, var)
    end

    isempty(vars_spectra) && return nothing

    return make_plots_generic(
        output_path,
        plot_path,
        vars_spectra;
        output_name = output_prefix * "spectra",
        plot_fn = plot_spectrum_with_line!,
    )
end

"""
    make_diagnostics_plots(
        output_path::AbstractString,
        plot_path::AbstractString;
        output_prefix = "",
    )
Create plots for diagnostics. The plots are saved to `plot_path`.
This function will plot all variables that have been saved in `output_path`.
The `reduction` keyword argument should be consistent with the reduction used to save the diagnostics.
"""
function Plotting.make_diagnostics_plots(
    output_path::AbstractString,
    plot_path::AbstractString;
    output_prefix = "",
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
        "average" in reductions ? (reduction = "average") : (reduction = first(reductions))
        periods = CAN.available_periods(simdir; short_name, reduction)
        "1d" in periods ? (period = "1d") : (period = first(periods))
        vars[i] = get(simdir; short_name, reduction, period)
    end

    # Filter vars into 2D and 3D variable diagnostics vectors
    # 3D fields are zonally averaged platted on the lat-z plane
    # 2D fields are plotted on the lon-lat plane
    vars_3D =
        map(var_3D -> CAN.average_lon(var_3D), filter(var -> CAN.has_altitude(var), vars))
    vars_2D = filter(var -> !CAN.has_altitude(var), vars)

    # Generate plots and save in `plot_path`
    !isempty(vars_3D) && make_plots_generic(
        output_path,
        plot_path,
        vars_3D,
        time = LAST_SNAP,
        output_name = output_prefix * "summary_3D",
        more_kwargs = YLINEARSCALE,
    )
    !isempty(vars_2D) && make_plots_generic(
        output_path,
        plot_path,
        vars_2D,
        time = LAST_SNAP,
        output_name = output_prefix * "summary_2D",
    )
end

"""
    make_ocean_diagnostics_plots(output_path::AbstractString, plot_path::AbstractString; output_prefix = "")

Create plots for diagnostics. The plots are saved to `ocean_summary_2D.pdf` in `plot_path`.
This function will plot the following variables, if they have been saved in `output_path`:
    - Temperature (`T`)
    - Salinity (`S`)
    - Zonal velocity (`u`)
    - Meridional velocity (`v`)

For each variable, take the surface level (top level) of the variable
and create a 2D plot. The plots will be saved in a single PDF file.
"""
function Plotting.make_ocean_diagnostics_plots(
    output_path::AbstractString,
    plot_path::AbstractString;
    output_prefix = "",
)
    expected_output_path = joinpath(output_path, "ocean_diagnostics.nc")
    isfile(expected_output_path) || return nothing

    # Create an OutputVar for each diagnostic, so we can use ClimaAnalysis to plot
    var_names = ["T", "S", "u", "v"]
    vars = Array{Union{CAN.OutputVar, Nothing}}(undef, length(var_names))
    for (i, var_name) in enumerate(var_names)
        # Create an OutputVar if the variable is available in the output file
        output_var = CAN.OutputVar(expected_output_path, var_name)
        output_var.attributes["short_name"] = var_name

        # Take the top level (surface) of the variable
        output_var = CAN.slice(output_var, z_aac = output_var.dims["z_aac"][1])

        vars[i] = output_var
    end

    # Filter out any variables that are not available
    vars = filter(!isnothing, vars)

    # Make plots for each variable, saved in one PDF file
    !isempty(vars) && make_plots_generic(
        expected_output_path, # file_path
        plot_path,
        vars,
        time = LAST_SNAP,
        output_name = output_prefix * "summary_2D",
    )
    return nothing
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
                normalized_path =
                    lpad(path, LABEL_LENGTH + 1, " ")[(end - LABEL_LENGTH):end]

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
