import ClimaComms
import ClimaAtmos as CA
import ClimaCore
import ClimaCore: DataLayouts, Fields, Geometry
import ClimaCore.Fields: Field, FieldVector, field_values
import ClimaCore.DataLayouts: AbstractData
import ClimaCore.Geometry: AxisTensor
import ClimaCore.Spaces: AbstractSpace
import ClimaUtilities.TimeVaryingInputs: AbstractTimeVaryingInput
import NCDatasets

# Allow cache to be moved on CPU
ClimaCore.Adapt.@adapt_structure CA.AtmosCache

function restore!(
    v1::T1,
    v2::T2;
    name = "",
    ignore = Set([:rc, :graph_context]),
) where {T1 <: Union{NamedTuple, CA.AtmosCache}, T2 <: Union{NamedTuple, CA.AtmosCache}}
    _restore!(v1, v2; name, ignore)
    return nothing
end

function _restore!(v1::T1, v2::T2; name, ignore) where {T1, T2}
    properties = filter(x -> !(x in ignore), propertynames(v1))
    if isempty(properties)
        if !Base.issingletontype(typeof(v1))
            _restore_base!(v1, v2; name, ignore)
        else
            v1 == v2 || error("$v1 != $v2")
        end
    else
        # Recursive case
        for p in properties
            _restore!(getproperty(v1, p), getproperty(v2, p); name = "$(name).$(p)", ignore)
        end
    end
    return nothing
end

# Ignoring certain types that don't need to be restored
function _restore!(v1::AbstractTimeVaryingInput, v2::AbstractTimeVaryingInput; name, ignore)
    return nothing
end

function _restore!(v1::ClimaComms.MPICommsContext, v2::ClimaComms.MPICommsContext; name, ignore)
    return nothing
end

function _restore_base!(v1::T, v2::T; name, ignore) where {T <: Number}
    # To account for NaN
    v1 === v2 || error("$v1 != $v2")
    return nothing
end

function _restore_base!(
    v1::T1,
    v2::T2;
    name,
    ignore,
) where {
    T1 <: Union{Field, FieldVector, AbstractData, AbstractArray},
    T2 <: Union{Field, FieldVector, AbstractData, AbstractArray},
}
    # HACK: Here we are making the assumption that the default context is being
    # used. This might not be true! This is currently needed because we don't
    # have access to CUDA.jl in this file
    ArrayType = parent(v1) isa Array ? Array : ClimaComms.array_type(ClimaComms.device())
    moved_to_device = ArrayType(parent(v2))

    parent(v1) .= moved_to_device
    return nothing
end
