import Printf
import ClimaCore as CC
import Makie
import ClimaCoreMakie
import CairoMakie
import ClimaCoupler: Interfacer, ConservationChecker, Plotting
import StaticArrays

"""
    plot_conservation(
        cc::ConservationChecker.ConservationCheck,
        coupler_sim::Interfacer.CoupledSimulation;
        figname = "conservation.png"
    )

Plot globally integrated conserved quantity with time.

The first panel is the evolution of the total conserved quantity and its
components relative to their initial values, while the second panel plots
the log of the relative error of the total conserved quantity.
"""
function Plotting.plot_conservation(
    cc::ConservationChecker.ConservationCheck,
    coupler_sim::Interfacer.CoupledSimulation;
    figname = "conservation.png",
)
    components = cc.components
    name = ConservationChecker.name(cc.conserved_quantity)
    units = ConservationChecker.units(cc.conserved_quantity)

    days = collect(1:length(components[1])) * float(coupler_sim.Δt_cpl) / 86400
    total = zeros(length(components[1]))
    abs_sum = 0.0 # this will be ∑ᵢ|xᵢ| at the final time step, used for relative error
    f = Makie.Figure(size = (800, 500))
    ax1 = Makie.Axis(f[1, 1], xlabel = "time (days)", ylabel = "Δ$name ($units)")
    for (component, timeseries) in pairs(components)
        Makie.lines!(ax1, days, timeseries .- timeseries[1], label = string(component))
        abs_sum += abs(timeseries[end])
        total .+= timeseries
    end
    Makie.lines!(
        ax1,
        days,
        total .- total[1],
        label = "total",
        color = :black,
        linestyle = :dash,
    )
    Makie.Legend(f[1, 2], ax1, framevisible = false)

    # use ∑ᵢ|xᵢ| at the final time step as a reference for the error
    err = abs.((total .- total[1]) ./ abs_sum)
    l_err = log10.(err)
    ax2 = Makie.Axis(
        f[2, 1],
        xlabel = "time (days)",
        ylabel = "log( |Δtotal| / ∑|components| )",
    )
    Makie.lines!(ax2, days, l_err, color = :black)
    l_err_valid = filter(x -> !isinf(x) && !isnan(x), l_err)
    if !isempty(l_err_valid)
        y_min = minimum(l_err_valid)
        y_max = maximum(l_err_valid)
        if y_min != y_max
            Makie.ylims!(ax2, y_min, y_max)
        else
            # If all values are the same, add a small padding to avoid Makie error
            padding = max(abs(y_min) * 0.01, 0.1)
            Makie.ylims!(ax2, y_min - padding, y_max + padding)
        end
    end
    Makie.resize_to_layout!(f)
    Makie.save(figname, f)

    # Summarize the budget at the end of the simulation with a table
    header = ("Component", "Value ($units)", "Change ($units)", "Change (%)")
    rows = [header]
    for (component, timeseries) in pairs((; components..., total))
        value = timeseries[end]
        change = timeseries[end] - timeseries[1]
        push!(
            rows,
            (
                string(component),
                Printf.@sprintf("%.3e", value),
                Printf.@sprintf("%.3e", change),
                iszero(value) ? "n/a" : Printf.@sprintf("%.2f", change / value * 100),
            ),
        )
    end
    widths = [maximum(length(row[i]) for row in rows) for i in eachindex(header)]
    format_row(row) =
        "\n  " *
        rpad(row[1], widths[1]) *
        join("  " * lpad(row[i], widths[i]) for i in 2:lastindex(widths))
    info_msg = "Total $name at end of simulation:"
    info_msg *= format_row(header)
    info_msg *= "\n  " * repeat("-", sum(widths) + 2 * (length(widths) - 1))
    for row in rows[2:end]
        info_msg *= format_row(row)
    end
    @info info_msg
end

"""
    debug(cs::Interfacer.CoupledSimulation, dir = "debug", cs_fields_ref = nothing)

Plot the fields of a coupled simulation and save plots to a directory.
"""
function Plotting.debug(
    cs::Interfacer.CoupledSimulation,
    dir = "debug",
    cs_fields_ref = nothing,
)
    isdir(dir) || mkpath(dir)
    @info "plotting debug in $dir"
    for sim in cs.model_sims
        Plotting.debug(sim, dir)
    end
    Plotting.debug(cs.fields, dir, cs_fields_ref)
end

"""
    debug(cs_fields::CC.Fields.Field, dir, cs_fields_ref = nothing)

Plot useful coupler fields (in `field_names`) and save plots to a directory.

If `cs_fields_ref` is provided (e.g., using a copy of cs.fields from the initialization),
plot the anomalies of the fields with respect to `cs_fields_ref`.

For vector fields which are not defined on a Cartesian basis, rotate them to the Cartesian
basis before plotting so they can be interpreted physically.

For single-column (PointSpace) fields, each subplot shows only the title with extrema;
no heatmap or table is drawn.
"""
function Plotting.debug(cs_fields::CC.Fields.Field, dir, cs_fields_ref = nothing)
    field_names = propertynames(cs_fields)

    fig = Makie.Figure(size = (1500, 800))
    min_square_len = ceil(Int, sqrt(length(field_names)))

    # Set up to rotate vector fields to the Cartesian basis
    local_geometry = CC.Fields.local_geometry_field(getproperty(cs_fields, field_names[1]))

    has_nan = false
    for i in 1:min_square_len, j in 1:min_square_len
        field_index = (i - 1) * min_square_len + j

        # Rotate vector fields to the Cartesian basis so they can be interpreted physically
        rotated_vectors = Dict{Symbol, CC.Fields.Field}()
        if :u_int in field_names && :v_int in field_names
            u_int_cartesian, v_int_cartesian = _get_cartesian_vector_components(
                local_geometry,
                getproperty(cs_fields, :u_int),
                getproperty(cs_fields, :v_int),
            )

            rotated_vectors[:u_int] = u_int_cartesian
            rotated_vectors[:v_int] = v_int_cartesian
        end
        if :F_turb_ρτxz in field_names && :F_turb_ρτyz in field_names
            F_turb_ρτxz_cartesian, F_turb_ρτyz_cartesian = _get_cartesian_vector_components(
                local_geometry,
                getproperty(cs_fields, :F_turb_ρτxz),
                getproperty(cs_fields, :F_turb_ρτyz),
            )

            rotated_vectors[:F_turb_ρτxz] = F_turb_ρτxz_cartesian
            rotated_vectors[:F_turb_ρτyz] = F_turb_ρτyz_cartesian
        end

        if field_index <= length(field_names)
            field_name = field_names[field_index]

            if field_name in keys(rotated_vectors)
                field = rotated_vectors[field_name]
            else
                field = getproperty(cs_fields, field_name)
            end

            extrema_str, field_has_nan = Plotting.print_extrema(field)
            has_nan = has_nan || field_has_nan
            title = string(field_name) * extrema_str
            ax = Makie.Axis(fig[i, j * 2 - 1]; title)
            Plotting.debug_plot!(ax, fig, field, i, j)
        end
    end
    mkpath(dir)
    Makie.save(joinpath(dir, "debug_coupler.png"), fig)

    # plot anomalies if a reference cs.fields, `cs_fields_ref`, are provided
    if !isnothing(cs_fields_ref)
        for i in 1:min_square_len, j in 1:min_square_len
            field_index = (i - 1) * min_square_len + j
            if field_index <= length(field_names)
                field_name = field_names[field_index]
                field = getproperty(cs_fields, field_name)

                extrema_str, field_has_nan = Plotting.print_extrema(field)
                has_nan = has_nan || field_has_nan
                title = string(field_name) * extrema_str
                ax = Makie.Axis(fig[i, j * 2 - 1]; title)
                Plotting.debug_plot!(ax, fig, field, i, j)
            end
        end
        Makie.save(joinpath(dir, "debug_coupler_anomalies.png"), fig)
    end

    # Check for NaN errors after plots are saved
    if has_nan
        @warn "NaN values found in coupler fields extrema"
    end
end

function _get_cartesian_vector_components(
    local_geometry,
    u_component::CC.Fields.Field,
    v_component::CC.Fields.Field,
)
    # Get the vector components in the CT1 and CT2 directions
    xz = @. CT12(CT1(_unit_basis_vector_data(CT1, local_geometry)), local_geometry)
    yz = @. CT12(CT2(_unit_basis_vector_data(CT2, local_geometry)), local_geometry)

    # Convert the vector components to a UVVector on the Cartesian basis
    uv_cartesian =
        @. CC.Geometry.UVVector(u_component * xz + v_component * yz, local_geometry)
    return uv_cartesian.components.data.:1, uv_cartesian.components.data.:2
end

"""
    debug(sim::Interfacer.AbstractComponentSimulation, dir)

Plot the fields of a component model simulation and save plots to a directory.
"""
function Plotting.debug(sim::Interfacer.AbstractComponentSimulation, dir)
    field_names = Plotting.debug_plot_fields(sim)
    fig = Makie.Figure(size = (1500, 800))
    min_square_len = ceil(Int, sqrt(length(field_names)))
    has_nan = false
    for i in 1:min_square_len, j in 1:min_square_len
        field_index = (i - 1) * min_square_len + j
        if field_index <= length(field_names)
            field_name = field_names[field_index]
            field = Interfacer.get_field(sim, Val(field_name))
            extrema_str, field_has_nan = Plotting.print_extrema(field)
            has_nan = has_nan || field_has_nan
            title = string(field_name) * extrema_str
            ax = Makie.Axis(fig[i, j * 2 - 1]; title)
            Plotting.debug_plot!(ax, fig, field, i, j)
        end
    end
    Makie.save(joinpath(dir, "debug_$(nameof(sim)).png"), fig)

    # Check for NaN errors after plots are saved
    if has_nan
        error("NaN values found in field extrema of $(nameof(sim))")
    end
end

"""
    Plotting.debug_plot!(ax, fig, field::CC.Fields.Field, i, j)

Helper function to plot a ClimaCore field in the given figure at position (i, j).
Dispatches on the field's space type via `_debug_plot_field!`.
"""
function Plotting.debug_plot!(ax, fig, field::CC.Fields.Field, i, j)
    cpu_field = CC.to_cpu(field)
    _debug_plot_field!(ax, fig, cpu_field, axes(cpu_field), i, j)
end

"""
    _debug_plot_field!(ax, fig, cpu_field, space::CC.Spaces.FiniteDifferenceSpace, i, j)

Plot a 1D column field (center or face FD space) as a line with z on the y-axis.
"""
function _debug_plot_field!(
    ax,
    fig,
    cpu_field,
    space::CC.Spaces.FiniteDifferenceSpace,
    i,
    j,
)
    _debug_plot_column!(ax, cpu_field)
    return nothing
end

"""
    _debug_plot_field!(ax, fig, cpu_field, space::CC.Spaces.ExtrudedFiniteDifferenceSpace, i, j)

Plot a 3D extruded field: if horizontal space is cubed-sphere, take level 1 and heatmap;
otherwise extract a single column and plot as a line (for BoxGrid column mode, used by atmos).
"""
function _debug_plot_field!(
    ax,
    fig,
    cpu_field,
    space::CC.Spaces.ExtrudedFiniteDifferenceSpace,
    i,
    j,
)
    hspace = CC.Spaces.horizontal_space(space)
    if hspace isa CC.Spaces.CubedSphereSpectralElementSpace2D
        # Cubed sphere case: make heatmap of level 1
        _debug_plot_heatmap!(ax, fig, CC.Fields.level(cpu_field, 1), i, j)
    else
        # BoxGrid column mode: extract a single column and plot as a line
        col_field = CC.Fields.column(cpu_field, 1, 1, 1)
        _debug_plot_column!(ax, col_field)
    end
    return nothing
end

"""
    _debug_plot_field!(ax, fig, cpu_field, space, i, j)

Fallback: plot a 2D spatial field as a heatmap (e.g. `SpectralElementSpace2D`).
"""
function _debug_plot_field!(ax, fig, cpu_field, space, i, j)
    _debug_plot_heatmap!(ax, fig, cpu_field, i, j)
    return nothing
end

"""
    _debug_plot_column!(ax, field)

Plot a 1D column field as a line with z on the y-axis.
"""
function _debug_plot_column!(ax, field)
    z = vec(Array(parent(CC.Fields.coordinate_field(field).z)))
    vals = vec(Array(parent(field)))

    valid = @. !(isnan(vals) || isinf(vals))
    any(valid) || return nothing

    Makie.lines!(ax, vals[valid], z[valid])
    ax.ylabel = "z"
    return nothing
end

"""
    _debug_plot_heatmap!(ax, fig, cpu_field, i, j)

Plot a 2D spatial field as a heatmap with colorbar.
"""
function _debug_plot_heatmap!(ax, fig, cpu_field, i, j)
    FT = CC.Spaces.undertype(axes(cpu_field))
    isinvalid = (x) -> isnan(x) || isinf(x)
    field_valid_min, field_valid_max =
        extrema(map(x -> isinvalid(x) ? FT(0) : x, parent(cpu_field)))
    map!(x -> isinvalid(x) ? 100field_valid_max : x, parent(cpu_field), parent(cpu_field))

    # If the values are too small, `isapprox` can't be computed accurately because of floating point precision issues.
    is_toosmall = (x) -> log10(abs(x)) < log10(floatmin(Float64)) / 2

    # If the field is constant, skip plotting it to avoid heatmap errors.
    if isapprox(field_valid_min, field_valid_max) ||
       (is_toosmall(field_valid_min) && is_toosmall(field_valid_max))
        return nothing
    else
        colorrange = (field_valid_min, field_valid_max)
        hm = ClimaCoreMakie.fieldheatmap!(ax, cpu_field, colorrange = colorrange)
        Makie.Colorbar(fig[i, j * 2], hm)
    end
    return nothing
end

"""
    Plotting.debug_plot!(ax, fig, field, i, j)

Make a line plot of the provided array.
"""
function Plotting.debug_plot!(ax, fig, field::AbstractArray, i, j)
    return Makie.lines!(ax, Array(field)) # This isn't really a heatmap, but it's okay for debugging
end
Plotting.debug_plot!(ax, fig, field, i, j) = nothing # fallback method

"""
    Plotting.print_extrema(field::Union{CC.Fields.Field, Vector, StaticArrays.SVector, Number})

Return a tuple `(string, has_nan::Bool)` where:
- `string` is the minimum and maximum values of a field formatted as a string
- `has_nan` is true if the extrema contain NaN values
"""
function Plotting.print_extrema(
    field::Union{CC.Fields.Field, Vector, StaticArrays.SVector, Number},
)
    ext_vals = extrema(field)
    min_val = ext_vals[1]
    max_val = ext_vals[2]

    # Check for NaN values
    has_nan = isnan(min_val) || isnan(max_val)

    min = Printf.@sprintf("%.2E", min_val)
    max = Printf.@sprintf("%.2E", max_val)
    return (" [$min, $max]", has_nan)
end

"""
    Plotting.debug_plot_fields(sim::Interfacer.AbstractSurfaceSimulation)

Return the default fields to include in debug plots for a surface model.
This should be extended for any atmosphere model, and any surface model
that has additional fields to plot.
"""
Plotting.debug_plot_fields(sim::Interfacer.AbstractSurfaceSimulation) =
    (:area_fraction, :surface_temperature)

# Define shorthands for ClimaCore types
const CT1 = CC.Geometry.Contravariant1Vector
const CT2 = CC.Geometry.Contravariant2Vector
const CT12 = CC.Geometry.Contravariant12Vector

"""
    _unit_basis_vector_data(type, local_geometry)

The component of the vector of the specified type with length 1 in physical units.
The type should correspond to a vector with only one component, i.e., a basis vector.
"""
function _unit_basis_vector_data(::Type{V}, local_geometry) where {V}
    FT = CC.Geometry.undertype(typeof(local_geometry))
    return FT(1) / CC.Geometry._norm(V(FT(1)), local_geometry)
end
