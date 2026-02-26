import Printf
import ClimaCore as CC
import Makie
import ClimaCoreMakie
import CairoMakie
import ClimaCoupler: Interfacer, ConservationChecker, Plotting
import StaticArrays

"""
    plot_global_conservation(
        cc::AbstractConservationCheck,
        coupler_sim::Interfacer.CoupledSimulation,
        softfail = false;
        figname1 = "total.png",
        figname2 = "total_log.png",
    )

Creates two plots of the globally integrated quantity (energy, ``\\rho e``):
1. global quantity of each model component as a function of time,
relative to the initial value;
2. fractional change in the sum of all components over time on a log scale.

Conservation checks are available for energy and water, and can be enabled by
running a Slabplanet simulation with `energy_check` set to true.

If `softfail` is false, asserts that the relative error in conservation of the
provided quantity is less than a pre-determined threshold. This argument
is controlled by the `conservation_softfail` simulation flag.
"""
function Plotting.plot_global_conservation(
    cc::ConservationChecker.AbstractConservationCheck,
    coupler_sim::Interfacer.CoupledSimulation,
    softfail = false;
    figname1 = "total.png",
    figname2 = "total_log.png",
)

    model_sims = coupler_sim.model_sims
    ccs = cc.sums

    days = collect(1:length(ccs[1])) * float(coupler_sim.Δt_cpl) / 86400

    # evolution of energy of each component relative to initial value
    total = ccs.total  # total

    var_name = nameof(cc)
    cum_total = [0.0]
    f = Makie.Figure()
    ax = Makie.Axis(f[1, 1], xlabel = "time [days]", ylabel = "$var_name: (t) - (t=0)")
    Makie.lines!(ax, days, total .- total[1], label = "total"; linewidth = 3)
    for sim in model_sims
        sim_name = nameof(sim)
        global_field = getproperty(ccs, Symbol(sim_name))
        diff_global_field = (global_field .- global_field[1])
        Makie.lines!(ax, days, diff_global_field[1:length(days)], label = sim_name)
        cum_total .+= abs.(global_field[end])
    end
    if cc isa ConservationChecker.EnergyConservationCheck
        global_field = ccs.toa_net_source
        diff_global_field = (global_field .- global_field[1])
        Makie.lines!(ax, days, diff_global_field[1:length(days)], label = "toa_net")
        cum_total .+= abs.(global_field[end])
    end
    Makie.axislegend(ax, position = :lb)
    Makie.save(figname1, f)

    # use the cumulative global sum at the final time step as a reference for the error calculation
    rse = abs.((total .- total[1]) ./ cum_total)
    l_rse = log.(rse)
    # evolution of log error of total
    lp = Makie.lines(days, l_rse, label = "rs error")
    lp.axis.xlabel = "time [days]"
    lp.axis.ylabel = "log( |x(t) - x(t=0)| / Σx(t=T) )"
    l_rse_valid = filter(x -> !isinf(x) && !isnan(x), l_rse)
    if !isempty(l_rse_valid)
        y_min = minimum(l_rse_valid)
        y_max = maximum(l_rse_valid)
        if y_min != y_max
            Makie.ylims!(y_min, y_max)
        else
            # If all values are the same, add a small padding to avoid Makie error
            padding = max(abs(y_min) * 0.01, 0.1)
            Makie.ylims!(y_min - padding, y_max + padding)
        end
    end
    Makie.axislegend(position = :lt)
    Makie.save(figname2, lp)

    # check that the relative error is small (TODO: reduce this to sqrt(eps(FT)))
    if !softfail
        @info typeof(cc)
        @info rse[end]
        @assert rse[end] < 0.035
    end
end


# plotting functions for the coupled simulation
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
"""
function Plotting.debug(cs_fields::CC.Fields.Field, dir, cs_fields_ref = nothing)
    field_names = propertynames(cs_fields)
    fig = Makie.Figure(size = (1500, 800))
    min_square_len = ceil(Int, sqrt(length(field_names)))

    # Set up to rotate vector fields to the Cartesian basis
    local_geometry = CC.Fields.local_geometry_field(getproperty(cs_fields, field_names[1]))

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

            title = string(field_name) * Plotting.print_extrema(field)
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

                title = string(field_name) * Plotting.print_extrema(field)
                ax = Makie.Axis(fig[i, j * 2 - 1]; title)
                Plotting.debug_plot!(ax, fig, field, i, j)
            end
        end
        Makie.save(joinpath(dir, "debug_coupler_anomalies.png"), fig)
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
    for i in 1:min_square_len, j in 1:min_square_len
        field_index = (i - 1) * min_square_len + j
        if field_index <= length(field_names)
            field_name = field_names[field_index]
            field = Interfacer.get_field(sim, Val(field_name))
            title = string(field_name) * Plotting.print_extrema(field)
            ax = Makie.Axis(fig[i, j * 2 - 1]; title)
            Plotting.debug_plot!(ax, fig, field, i, j)
        end
    end
    Makie.save(joinpath(dir, "debug_$(nameof(sim)).png"), fig)
end

"""
    Plotting.debug_plot!(ax, fig, field::CC.Fields.Field, i, j)

Helper function to plot a heatmap of a ClimaCore field in the given figure at position (i, j).
If the field is constant, skip plotting it to avoid heatmap errors.
"""
function Plotting.debug_plot!(ax, fig, field::CC.Fields.Field, i, j)
    # Copy field onto cpu space if necessary
    cpu_field = CC.to_cpu(field)
    if cpu_field isa CC.Fields.ExtrudedCubedSphereSpectralElementField3D
        cpu_field = CC.Fields.level(cpu_field, 1)
    end

    # ClimaCoreMakie doesn't support NaNs/Infs, so we substitute them with 100max
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

Return the minimum and maximum values of a field as a string.
"""
function Plotting.print_extrema(
    field::Union{CC.Fields.Field, Vector, StaticArrays.SVector, Number},
)
    ext_vals = extrema(field)
    min = Printf.@sprintf("%.2E", ext_vals[1])
    max = Printf.@sprintf("%.2E", ext_vals[2])
    return " [$min, $max]"
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
