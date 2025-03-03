# Define shared methods to allow reading back a saved cache

import ClimaComms
import ClimaCore
import ClimaCore: DataLayouts, Fields, Geometry
import ClimaCore.Fields: Field, FieldVector, field_values
import ClimaCore.DataLayouts: AbstractData
import ClimaCore.Geometry: AxisTensor
import ClimaCore.Spaces: AbstractSpace
import ClimaUtilities.TimeVaryingInputs: AbstractTimeVaryingInput
import StaticArrays
import NCDatasets

"""
    restore!(v1, v2, comms_ctx; ignore)

Recursively traverse `v1` and `v2`, setting each field of `v1` with the
corresponding field in `v2`. In this, ignore all the properties that have name
within the `ignore` iterable.

`ignore` is useful when there are stateful properties, such as live pointers.
"""
function restore!(v1::T1, v2::T2, comms_ctx; name = "", ignore) where {T1, T2}
    # We pick fieldnames(T2) because v2 tend to be simpler (Array as opposed
    # to CuArray)
    fields = filter(x -> !(x in ignore), fieldnames(T2))
    if isempty(fields)
        if !Base.issingletontype(typeof(v1))
            restore!(v1, v2, comms_ctx; name, ignore)
        else
            v1 == v2 || error("$v1 != $v2")
        end
    else
        # Recursive case
        for p in fields
            restore!(getfield(v1, p), getfield(v2, p), comms_ctx; name = "$(name).$(p)", ignore)
        end
    end
    return nothing
end

# Ignoring certain types that don't need to be restored
# UnionAll and DataType are infinitely recursive, so we also ignore those
function restore!(
    v1::Union{AbstractTimeVaryingInput, ClimaComms.AbstractCommsContext, ClimaComms.AbstractDevice, UnionAll, DataType},
    v2::Union{AbstractTimeVaryingInput, ClimaComms.AbstractCommsContext, ClimaComms.AbstractDevice, UnionAll, DataType},
    _comms_ctx;
    name,
    ignore,
)
    return nothing
end

function restore!(
    v1::T1,
    v2::T2,
    comms_ctx;
    name,
    ignore,
) where {T1 <: Union{AbstractData, AbstractArray}, T2 <: Union{AbstractData, AbstractArray}}
    ArrayType = parent(v1) isa Array ? Array : ClimaComms.array_type(ClimaComms.device(comms_ctx))
    moved_to_device = ArrayType(parent(v2))

    parent(v1) .= moved_to_device
    return nothing
end

function restore!(
    v1::T1,
    v2::T2,
    comms_ctx;
    name,
    ignore,
) where {
    T1 <: Union{StaticArrays.StaticArray, Number, UnitRange, Symbol},
    T2 <: Union{StaticArrays.StaticArray, Number, UnitRange, Symbol},
}
    v1 == v2 || error("$name is a immutable but it inconsistent")
    return nothing
end

function restore!(v1::T1, v2::T2, comms_ctx; name, ignore) where {T1 <: Dict, T2 <: Dict}
    # RRTGMP has some internal dictionaries
    v1 == v2 || error("$name is inconsistent")
    return nothing
end
