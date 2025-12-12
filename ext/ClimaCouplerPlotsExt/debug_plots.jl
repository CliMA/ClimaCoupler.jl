import Printf
import ClimaCore as CC
import Makie
import ClimaCoreMakie
import CairoMakie
import ClimaCoupler: Interfacer, ConservationChecker
import ClimaAtmos as CA
import Oceananigans as OC
import StaticArrays

export debug

"""
    debug(cs::Interfacer.CoupledSimulation, dir = "debug", cs_fields_ref = nothing)

Plot the fields of a coupled simulation and save plots to a directory.
"""
function debug(cs::Interfacer.CoupledSimulation, dir = "debug", cs_fields_ref = nothing)
    isdir(dir) || mkpath(dir)
    @info "plotting debug in $dir"
    for sim in cs.model_sims
        debug(sim, dir)
    end
    debug(cs.fields, dir, cs_fields_ref)
end

"""
    debug(cs_fields::CC.Fields.Field, dir, cs_fields_ref = nothing)

Plot useful coupler fields (in `field_names`) and save plots to a directory.

If `cs_fields_ref` is provided (e.g., using a copy of cs.fields from the initialization),
plot the anomalies of the fields with respect to `cs_fields_ref`.
"""
function debug(cs_fields::CC.Fields.Field, dir, cs_fields_ref = nothing)
    field_names = propertynames(cs_fields)
    fig = Makie.Figure(size = (1500, 800))
    min_square_len = ceil(Int, sqrt(length(field_names)))
    for i in 1:min_square_len, j in 1:min_square_len
        field_index = (i - 1) * min_square_len + j
        if field_index <= length(field_names)
            field_name = field_names[field_index]
            field = getproperty(cs_fields, field_name)
            _heatmap_cc_field!(fig, field, i, j, field_name)
        end
    end
    Makie.save(joinpath(dir, "debug_coupler.png"), fig)

    # plot anomalies if a reference cs.fields, `cs_fields_ref`, are provided
    if !isnothing(cs_fields_ref)
        for i in 1:min_square_len, j in 1:min_square_len
            field_index = (i - 1) * min_square_len + j
            if field_index <= length(field_names)
                field_name = field_names[field_index]
                field = getproperty(cs_fields, field_name)
                _heatmap_cc_field!(fig, field, i, j, field_name)
            end
        end
        Makie.save(joinpath(dir, "debug_coupler_anomalies.png"), fig)
    end
end

"""
    debug(sim::Interfacer.ComponentModelSimulation, dir)

Plot the fields of a component model simulation and save plots to a directory.
"""
function debug(sim::Interfacer.ComponentModelSimulation, dir)
    field_names = plot_field_names(sim)
    fig = Makie.Figure(size = (1500, 800))
    min_square_len = ceil(Int, sqrt(length(field_names)))
    for i in 1:min_square_len, j in 1:min_square_len
        field_index = (i - 1) * min_square_len + j
        if field_index <= length(field_names)
            field_name = field_names[field_index]
            field = Interfacer.get_field(sim, Val(field_name))
            # If field is a ClimaCore field, then _heatmap_cc_field! will add a
            # title to the axis
            title =
                field isa CC.Fields.Field ? "" : string(field_name) * print_extrema(field)
            ax = Makie.Axis(fig[i, j * 2 - 1]; title)
            if field isa OC.Field || field isa OC.AbstractOperations.AbstractOperation
                if field isa OC.AbstractOperations.AbstractOperation
                    field = OC.Field(field)
                    OC.compute!(field)
                end
                grid = field.grid
                hm = CairoMakie.heatmap!(ax, view(field, :, :, grid.Nz))
                Makie.Colorbar(fig[i, j * 2], hm)
            elseif field isa CC.Fields.Field
                _heatmap_cc_field!(fig, field, i, j, field_name)
            elseif field isa AbstractArray
                lin = Makie.lines!(ax, Array(field))
            end
        end
    end
    Makie.save(joinpath(dir, "debug_$(nameof(sim)).png"), fig)
end

"""
    _heatmap_cc_field!(fig, field::CC.Fields.Field, i, j, field_name)

Helper function to plot a heatmap of a ClimaCore field in the given figure at position (i, j).
If the field is constant, skip plotting it to avoid heatmap errors.
"""
function _heatmap_cc_field!(fig, field::CC.Fields.Field, i, j, field_name)
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
        ax = Makie.Axis(
            fig[i, j * 2 - 1],
            title = string(field_name) * print_extrema(cpu_field),
        )
    else
        colorrange = (field_valid_min, field_valid_max)
        ax = Makie.Axis(
            fig[i, j * 2 - 1],
            title = string(field_name) * print_extrema(cpu_field),
        )
        hm = ClimaCoreMakie.fieldheatmap!(ax, cpu_field, colorrange = colorrange)
        Makie.Colorbar(fig[i, j * 2], hm)
    end
    return nothing
end

"""
    print_extrema(field::CC.Fields.Field)

Return the minimum and maximum values of a field as a string.
"""
function print_extrema(field::Union{CC.Fields.Field, Vector, StaticArrays.SVector})
    ext_vals = extrema(field)
    min = Printf.@sprintf("%.2E", ext_vals[1])
    max = Printf.@sprintf("%.2E", ext_vals[2])
    return " [$min, $max]"
end

function print_extrema(num::Number)
    min = Printf.@sprintf("%.2E", num)
    max = Printf.@sprintf("%.2E", num)
    return " [$min, $max]"
end

function print_extrema(operation::OC.AbstractOperations.AbstractOperation)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    return print_extrema(evaluated_field)
end

function print_extrema(field::OC.Field)
    min = Printf.@sprintf("%.2E", minimum(field))
    max = Printf.@sprintf("%.2E", maximum(field))
    return " [$min, $max]"
end

# currently selected plot fields
plot_field_names(sim::Interfacer.SurfaceModelSimulation) =
    (:area_fraction, :surface_temperature)
