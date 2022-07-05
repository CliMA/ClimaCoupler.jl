using ClimaCore
using ClimaCore: Operators, Fields
using PrettyTables

export CouplerState
export coupler_push!, coupler_pull!, coupler_put!, coupler_get, coupler_get!
export coupler_add_field!, coupler_add_map!

mutable struct CplFieldInfo{DT, MD}
    # the coupled data
    data::DT
    # the name of the model that has permission to write to data
    write_sim::Symbol
    # catch-all metadata container
    metadata::MD
end

mutable struct CouplerState{CF, RO}
    coupled_fields::CF
    remap_operators::RO
end

_fields(coupler::CouplerState) = getfield(coupler, :coupled_fields)

"""
    CouplerState()

Type for holding coupler "state". This is the namespace through which coupled components
communicate. Its role is to provide a level of indirection so that components remain modular
and so that any data communication, interpolation, reindexing/unit conversions and filtering 
etc... can be embeded in the intermdediate coupling layer.

A field is exported by one component and imported by one or more other components.
"""
function CouplerState()
    return CouplerState(Dict{Symbol, CplFieldInfo}(), Dict{Symbol, Operators.LinearRemap}())
end

"""
    coupler_add_field!(
            coupler::CouplerState,
            fieldname::Symbol,
            fieldvalue,
        )

Add a field to the coupler that is accessible with key `fieldname`. 

# Arguments
- `coupler`: coupler object the field is added to.
- `fieldname`: key to access the field in the coupler.
- `fieldvalue`: data array of field values.
"""
function coupler_add_field!(
    coupler::CouplerState,
    fieldname::Symbol,
    fieldvalue;
    write_sim::AbstractSimulation,
    metadata = nothing,
)
    push!(coupler.coupled_fields, fieldname => CplFieldInfo(fieldvalue, name(write_sim), metadata))
end

"""
    coupler_add_map!(
            coupler::CouplerState,
            map_name::Symbol,
            map::Operators.LinearRemap
        )

Add a map to the coupler that is accessible with key `mapname`. 

# Arguments
- `coupler`: coupler object the field is added to.
- `mapname`: key to access the map in the coupler's map list.
- `map`: a remap operator.
"""
function coupler_add_map!(coupler::CouplerState, map_name::Symbol, map::Operators.LinearRemap)
    push!(coupler.remap_operators, map_name => map)
end

"""
    coupler_put!(coupler::CouplerState, fieldname::Symbol, fieldvalue)

Sets coupler field `fieldname` to `fieldvalue`.
"""
function coupler_put!(coupler::CouplerState, fieldname::Symbol, fieldvalue, source_sim::AbstractSimulation)
    cplfield = coupler.coupled_fields[fieldname]
    @assert cplfield.write_sim == name(source_sim) "$fieldname can only be written to by $(cplfield.write_sim)."

    cplfield.data .= fieldvalue

    return nothing
end

"""
    coupler_push!(coupler::CouplerState, model)

Update coupler with fields retrieved from the coupler.

`coupler_push!` is an adapter function to be implemented for each
model component using the coupler. It should send coupling fields via
`coupler_put!` calls and perform any operations on these fields to prepare
them for the coupler.
"""
function coupler_push!(coupler::CouplerState, model) end

"""
    coupler_get!(target_field::ClimaCore.Fields.Field, coupler::CouplerState, fieldname::Symbol, target_sim::AbstractSimulation)

Retrieve data array corresponding to `fieldname`, remap and store in `target_field`.
"""
function coupler_get!(
    target_field::ClimaCore.Fields.Field,
    coupler::CouplerState,
    fieldname::Symbol,
    target_sim::AbstractSimulation,
)
    cplfield = coupler.coupled_fields[fieldname]
    map = get_remap_operator(coupler, name(target_sim), cplfield.write_sim)
    Operators.remap!(target_field, map, cplfield.data)
end

"""
    coupler_get(coupler::CouplerState, fieldname::Symbol [, target_sim::AbstractSimulation])

Retrieve data array corresponding to `fieldname`.

If a `target_sim` is passed, the field is remapped to that simulation's boundary space.
"""
function coupler_get(coupler::CouplerState, fieldname::Symbol)
    cplfield = coupler.coupled_fields[fieldname]
    return cplfield.data
end

function coupler_get(coupler::CouplerState, fieldname::Symbol, target_sim::AbstractSimulation)
    cplfield = coupler.coupled_fields[fieldname]
    map = get_remap_operator(coupler, name(target_sim), cplfield.write_sim)
    return Operators.remap(map, cplfield.data)
end

"""
    coupler_pull!(model, coupler::CouplerState)

Update model with fields retrieved from the coupler.

`coupler_pull!` is an adapter function to be implemented for each
model component using the coupler. It should get coupling fields via
`coupler_get` calls and perform any operations on these fields to prepare
them for use in the component model.
"""
function coupler_pull!(model, coupler::CouplerState) end

function get_remap_operator(coupler, target_sim_name::Symbol, source_sim_name::Symbol)
    op_name = Symbol(source_sim_name, "_to_", target_sim_name)
    return coupler.remap_operators[op_name]
end

# display table of registered fields & their info
function Base.show(io::IO, coupler::CouplerState)
    #=
    ----------------
    |  Field Name  |
    ----------------
    |  :OceanSST   |
    |  :Field2     |
    ----------------
    =#
    fields = _fields(coupler)
    data = Array{Any}(undef, length(fields), 1)
    for (i, k) in enumerate(keys(fields))
        entry = fields[k]
        data[i, :] = [k]
    end
    header = (["Field Name"], [""])
    pretty_table(data, header = header)
end
