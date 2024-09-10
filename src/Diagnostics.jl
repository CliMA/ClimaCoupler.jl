"""
    Diagnostics

This module contains functions for defining, gathering and outputting online model diagnostics from the Coupler.
"""
module Diagnostics
import Dates
import ClimaCore as CC
import ClimaUtilities: CallbackManager
import ..Interfacer

export get_var, init_diagnostics, accumulate_diagnostics!, save_diagnostics, TimeMean

"""
    AbstractOutputGroup

Abstract type for ClimaCoupler's output diagnostics groups. Each diagnostic group should
contain fields that are of the same type and size, so the extended methods for the group's
operation functions work in the same way for all the fields.
"""
abstract type AbstractOutputGroup end

"""
    DiagnosticsGroup{S, NTO <: NamedTuple}

Defines a concrete diagnostics group type with fields `field_vector`, `operations`, `save`,
`output_dir` and `name_tag`.
"""
struct DiagnosticsGroup{S, NTO <: NamedTuple} <: AbstractOutputGroup
    field_vector::CC.Fields.FieldVector
    operations::NTO
    save::S
    output_dir::String
    name_tag::String
end

"""
    AbstractDiagnosticsOperations

Abstract type for operations to be performed on ClimaCoupler's diagnostics.
"""
abstract type AbstractDiagnosticsOperations end

"""
    TimeMean{C}

Defines a concrete operation type for time-averaged diagnostics. The counter `ct`
is used to accumulate the sum of the diagnostics.
"""
struct TimeMean{C} <: AbstractDiagnosticsOperations
    ct::C
end
# TODO: TemporalVariances, SpatialMean, SpatialVariances

"""
    function init_diagnostics(
        names::Tuple,
        space::CC.Spaces.AbstractSpace;
        save = CallbackManager.EveryTimestep(),
        operations = (;),
        output_dir = "",
        name_tag = "",
    )

Initializes diagnostics groups.
"""
function init_diagnostics(
    names::Tuple,
    space::CC.Spaces.AbstractSpace;
    save = CallbackManager.EveryTimestep(),
    operations = (;),
    output_dir = "",
    name_tag = "",
)
    data = NamedTuple{names}(ntuple(i -> CC.Fields.zeros(space), length(names)))
    return DiagnosticsGroup(CC.Fields.FieldVector(; data...), operations, save, output_dir, name_tag)
end

"""
    get_var(cs::Interfacer.CoupledSimulation, x)

Defines variable extraction from the coupler simulation. User specific diagnostics
should extend this function in the experiments folder.

Example:

get_var(cs, ::Val{:T_sfc}) = cs.fields.T_S

"""
get_var(::Interfacer.CoupledSimulation, x) = @warn "Variable $x is not defined."

"""
    accumulate_diagnostics!(cs::Interfacer.CoupledSimulation)

Accumulates user-defined diagnostics listed in the in the `field_vector` of each `dg`.
"""
function accumulate_diagnostics!(cs::Interfacer.CoupledSimulation)
    for dg in cs.diagnostics
        if dg.operations.accumulate !== nothing
            # TODO: avoid collecting at each timestep where not needed
            iterate_operations(cs, dg, collect_diags(cs, dg))
        end
    end
end

"""
    collect_diags(cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup)

Collects diagnostics in diags names.
"""
function collect_diags(cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup)
    FT = Interfacer.float_type(cs)
    diags = (;)

    diag_names = propertynames(dg.field_vector)
    for name in diag_names
        diags = (; diags..., zip((name,), (FT.(get_var(cs, Val(name))),))...)
    end

    return CC.Fields.FieldVector(; diags...)
end

"""
    iterate_operations(cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup, diags::CC.Fields.FieldVector)

Applies iteratively all specified diagnostics operations.
"""
function iterate_operations(cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup, new_diags::CC.Fields.FieldVector)
    for op in dg.operations
        operation(cs, dg, new_diags, op)
    end
end

"""
    operation(cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup, new_diags::CC.Fields.FieldVector, ::TimeMean)

Accumulates in time all entries in `new_diags` and saves the result in `dg.field_vector`, while
increasing the `dg.ct` counter.
"""
function operation(::Interfacer.CoupledSimulation, dg::DiagnosticsGroup, new_diags::CC.Fields.FieldVector, ::TimeMean)
    dg.field_vector .+= new_diags
    dg.operations.accumulate.ct[1] += Int(1)

    return nothing
end

"""
    operation(cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup, new_diags::CC.Fields.FieldVector, ::Nothing)

Accumulates in time all entries in `new_diags` and saves the result in `dg.field_vector`, while
increasing the `dg.ct` counter.
"""
function operation(::Interfacer.CoupledSimulation, dg::DiagnosticsGroup, new_diags::CC.Fields.FieldVector, ::Nothing)
    dg.field_vector .= new_diags
    return nothing
end

"""
    save_diagnostics(cs::Interfacer.CoupledSimulation)

    save_diagnostics(cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup, output_dir::String)

Saves all entries in `dg` in separate HDF5 files per variable in `output_dir`.
"""
function save_diagnostics(cs::Interfacer.CoupledSimulation)
    for dg in cs.diagnostics

        # Check if the date is greater than the next date to save
        if cs.dates.date[1] >= cs.dates.date1[1]
            pre_save(dg.operations.accumulate, cs, dg)
            save_diagnostics(cs, dg)
            post_save(dg.operations.accumulate, cs, dg)
        end
    end
end
function save_diagnostics(cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup)

    date = cs.dates.date[1]
    tag = dg.name_tag
    diag_save = dg.save
    output_dir = dg.output_dir

    !isdir(output_dir) && mkpath(output_dir)

    diag_names = propertynames(dg.field_vector)
    diag_values = map(x -> getproperty(dg.field_vector, x), diag_names)

    @info "Saving coupler diagnostics:"

    for (name, values) in zip(diag_names, diag_values)
        output_file = joinpath(output_dir, "$name.$tag." * string(date) * ".hdf5")
        @info "    $output_file"
        hdfwriter = CC.InputOutput.HDF5Writer(output_file, cs.comms_ctx)
        CC.InputOutput.HDF5.write_attribute(hdfwriter.file, "unix time", save_time_format(date, diag_save))
        CC.InputOutput.write!(hdfwriter, values, string(name))
        Base.close(hdfwriter)
    end
    return nothing

end

"""
    save_time_format(date::Dates.DateTime, ::CallbackManager.Monthly)

Converts the DateTime `date` to the conventional Unix format (seconds elapsed since 00:00:00 UTC on 1 January 1970).
"""
function save_time_format(date::Dates.DateTime, ::CallbackManager.Monthly)
    date_m1 = date - Dates.Day(1) # obtain previous month
    datetime = Dates.DateTime(Dates.yearmonth(date_m1)[1], Dates.yearmonth(date_m1)[2])
    Dates.datetime2unix(datetime)
end

save_time_format(date::Dates.DateTime, ::CallbackManager.EveryTimestep) = Dates.datetime2unix(date)

"""
    pre_save(::TimeMean, cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup)

Divides the accumulated sum by 'ct' to form the mean, before saving the diagnostics.
"""
function pre_save(::TimeMean, cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup)
    dg.field_vector .= dg.field_vector / dg.operations.accumulate.ct[1]
end

"""
    pre_save(::Nothing, cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup

Collects variables and performs all specified operations before saving the snapshot diagnostics.
"""
pre_save(::Nothing, cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup) =
    iterate_operations(cs, dg, collect_diags(cs, dg))

"""
    post_save(::TimeMean, cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup)

Resets accumulating fields and counts after saving the diagnostics.
"""
function post_save(::TimeMean, cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup)
    FT = eltype(dg.field_vector)
    dg.field_vector .= FT(0.0)
    dg.operations.accumulate.ct .= FT(0)
end

post_save(::Nothing, cs::Interfacer.CoupledSimulation, dg::DiagnosticsGroup) = nothing

end # module
