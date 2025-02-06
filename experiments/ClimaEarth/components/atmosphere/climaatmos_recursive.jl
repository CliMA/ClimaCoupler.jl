import ClimaComms
import ClimaAtmos as CA
import ClimaCore
import ClimaCore: DataLayouts, Fields, Geometry
import ClimaCore.Fields: Field, FieldVector, field_values
import ClimaCore.DataLayouts: AbstractData
import ClimaCore.Geometry: AxisTensor
import ClimaCore.Spaces: AbstractSpace
import NCDatasets

function restore!(
    v1::T,
    v2::T;
    name = "",
    ignore = Set([:rc]),
) where {T <: Union{NamedTuple, CA.AtmosCache}}
    _restore!(v1, v2; name, ignore)
    return nothing
end

function _restore!(v1::T, v2::T; name, ignore) where {T}
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

# function _restore!(v1::T, v2::T; name, ignore) where {T}
#     v1 .= v2
#     return nothing
# end

function _restore_base!(v1::T, v2::T; name, ignore) where {T <: Union{AbstractString, Symbol, CA.AtmosModel, Nothing}}
    v1 == v2 || error("$v1 != $v2")
    return nothing
end

function _restore_base!(v1::T, v2::T; name, ignore) where {T <: Number}
    # To account for NaN
    v1 === v2 || error("$v1 != $v2")
    return nothing
end

# # Ignore NCDatasets
# function _restore_base!(v1::T, v2::T; name, ignore) where {T <: NCDatasets.NCDataset}
#     return nothing
# end


function _restore_base!(v1::T, v2::T; name, ignore) where {T <: Union{Field, FieldVector, AbstractData, AbstractArray}}
    parent(v1) .= parent(v2)
    return nothing
end


# # We ignore NCDatasets. They contain a lot of state-ful information
# function _restore!(pass, v1::T, v2::T; name, ignore) where {T <: NCDatasets.NCDataset}
#     return nothing
# end

function _restore!(v1::T1, v2::T2; name, ignore) where {T1, T2}
    error("v1 and v2 have different types")
end
