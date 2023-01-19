"""
    Diagnostics

This module contains functions for defining, gathering and outputting online model diagnostics from the Coupler.
"""
module Diagnostics

using ClimaCore: Spaces, Fields, InputOutput
using ClimaCoupler.Utilities: CoupledSimulation
using Dates
using ClimaCoupler.TimeManager: AbstractFrequency, Monthly, EveryTimestep, trigger_callback
using ClimaComms

export get_var, init_diagnostics, accumulate_diagnostics!, save_diagnostics, TimeMean

abstract type AbstractOutputGroup end

struct DiagnosticsGroup{S, NTO <: NamedTuple} <: AbstractOutputGroup
    field_vector::Fields.FieldVector
    operations::NTO
    save::S
    output_dir::String
    name_tag::String
end

abstract type AbstractDiagnosticsOperations end

struct TimeMean{C} <: AbstractDiagnosticsOperations
    ct::C
end
# TODO: TemporalVariances, SpatialMean, SpatialVariances

"""
    init_diagnostics(names::Tuple, space::Spaces.Space; save = EveryTimestep(), operations = (;), output_dir = "", name_tag = "" )

Initializes diagnostics groups.
"""
function init_diagnostics(
    names::Tuple,
    space::Spaces.AbstractSpace;
    save = EveryTimestep(),
    operations = (;),
    output_dir = "",
    name_tag = "",
)
    data = NamedTuple{names}(ntuple(i -> Fields.zeros(space), length(names)))
    return DiagnosticsGroup(Fields.FieldVector(; data...), operations, save, output_dir, name_tag)
end

"""
    get_var(cs::CoupledSimulation, x)

Defines variable extraction from the coupler simulation. User specific diagnostics should extend this function in the experiments folder.

Example:

get_var(cs, ::Val{:T_sfc}) = cs.fields.T_S

"""
get_var(::CoupledSimulation, x) = @warn "Variable $x is not defined."

"""
    accumulate_diagnostics!(cs::CoupledSimulation)

Accumulates user-defined diagnostics listed in the in the `field_vector` of each `dg`.
"""
function accumulate_diagnostics!(cs::CoupledSimulation)
    for dg in cs.diagnostics
        if dg.operations.accumulate !== nothing
            iterate_operations(cs, dg, collect_diags(cs, dg)) # TODO: avoid collecting at each timestep where not needed
        end
    end
end

"""
    collect_diags(cs::CoupledSimulation, dg::DiagnosticsGroup)

Collects diagnostics in diags names.
"""
function collect_diags(cs::CoupledSimulation, dg::DiagnosticsGroup)

    diags = (;)

    diag_names = propertynames(dg.field_vector)
    for name in diag_names
        diags = (; diags..., zip((name,), (get_var(cs, Val(name)),))...)
    end

    return Fields.FieldVector(; diags...)
end

"""
    iterate_operations(cs::CoupledSimulation, dg::DiagnosticsGroup, diags::Fields.FieldVector)

Applies iteratively all specified diagnostics operations.
"""

function iterate_operations(cs::CoupledSimulation, dg::DiagnosticsGroup, new_diags::Fields.FieldVector)
    for op in dg.operations
        operation(cs, dg, new_diags, op)
    end
end

"""
    operation(cs::CoupledSimulation, dg::DiagnosticsGroup, new_diags::Fields.FieldVector, ::TimeMean)

Accumulates in time all entries in `new_diags` and saves the result in `dg.field_vector`, while increasing the `dg.ct` counter.
"""
function operation(::CoupledSimulation, dg::DiagnosticsGroup, new_diags::Fields.FieldVector, ::TimeMean)
    dg.field_vector .+= new_diags
    dg.operations.accumulate.ct[1] += Int(1)

    return nothing
end

"""
    operation(cs::CoupledSimulation, dg::DiagnosticsGroup, new_diags::Fields.FieldVector, ::Nothing)

Accumulates in time all entries in `new_diags` and saves the result in `dg.field_vector`, while increasing the `dg.ct` counter.
"""
function operation(::CoupledSimulation, dg::DiagnosticsGroup, new_diags::Fields.FieldVector, ::Nothing)
    dg.field_vector .= new_diags
    return nothing
end

"""
    save_diagnostics(cs::CoupledSimulation)

    save_diagnostics(cs::CoupledSimulation, dg::DiagnosticsGroup, output_dir::String)

Saves all entries in `dg` in separate HDF5 files per variable in `output_dir`.
"""
function save_diagnostics(cs::CoupledSimulation)
    for dg in cs.diagnostics
        if trigger_callback(cs, dg.save)
            pre_save(dg.operations.accumulate, cs, dg)
            save_diagnostics(cs, dg)
            post_save(dg.operations.accumulate, cs, dg)
        end
    end
end
function save_diagnostics(cs::CoupledSimulation, dg::DiagnosticsGroup)

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
        hdfwriter = InputOutput.HDF5Writer(output_file, cs.comms_ctx)
        InputOutput.HDF5.write_attribute(hdfwriter.file, "unix time", save_time_format(date, diag_save))
        InputOutput.write!(hdfwriter, values, string(name))
        Base.close(hdfwriter)
    end
    return nothing

end

"""
    save_time_format(date::Dates.DateTime, ::Monthly)

Converts the DateTime `date` to the conventional Unix format (seconds elapsed since 00:00:00 UTC on 1 January 1970).
"""
function save_time_format(date::Dates.DateTime, ::Monthly)
    date_m1 = date - Dates.Day(1) # obtain previous month
    datetime = Dates.DateTime(Dates.yearmonth(date_m1)[1], Dates.yearmonth(date_m1)[2])
    Dates.datetime2unix(datetime)
end

save_time_format(date::Dates.DateTime, ::EveryTimestep) = Dates.datetime2unix(date)

"""
    pre_save(::TimeMean, cs::CoupledSimulation, dg::DiagnosticsGroup)

Divides the accumulated sum by 'ct' to form the mean, before saving the diagnostics.
"""
function pre_save(::TimeMean, cs::CoupledSimulation, dg::DiagnosticsGroup)
    dg.field_vector .= dg.field_vector / dg.operations.accumulate.ct[1]
end

"""
    pre_save(::Nothing, cs::CoupledSimulation, dg::DiagnosticsGroup

Collects variables and performs all specified operations before saving the snapshot diagnostics.
"""
pre_save(::Nothing, cs::CoupledSimulation, dg::DiagnosticsGroup) = iterate_operations(cs, dg, collect_diags(cs, dg))

"""
    post_save(::TimeMean, cs::CoupledSimulation, dg::DiagnosticsGroup)

Resets accumulating fields and counts after saving the diagnostics.
"""
function post_save(::TimeMean, cs::CoupledSimulation, dg::DiagnosticsGroup)
    FT = eltype(dg.field_vector)
    dg.field_vector .= FT(0.0)
    dg.operations.accumulate.ct .= FT(0)
end

post_save(::Nothing, cs::CoupledSimulation, dg::DiagnosticsGroup) = nothing

end # module
