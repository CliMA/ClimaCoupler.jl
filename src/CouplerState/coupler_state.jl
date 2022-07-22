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

mutable struct CouplerState{FT, CF, RO}
    # A dictionary of fields added to the coupler
    coupled_fields::CF
    # A dictionary of remap operators between components
    remap_operators::RO
    # The coupled timestep size
    Δt_coupled::FT
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
function CouplerState(Δt_coupled)
    return CouplerState(Dict{Symbol, CplFieldInfo}(), Dict{Tuple{ClimaCore.Spaces.AbstractSpace, ClimaCore.Spaces.AbstractSpace}, Operators.LinearRemap}(), Δt_coupled)
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
    exchange_space = nothing,
    metadata = nothing,
)
    if exchange_space === nothing
        exchange_space = axes(fieldvalue) # will (eventually) build ClimaCore.Operators.IdentityRemap
    end
    map = ClimaCore.Operators.LinearRemap(exchange_space, axes(fieldvalue))
    coupler_add_map!(coupler, map)
    fieldvalue = Operators.remap(map, fieldvalue)

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
- `map`: a remap operator.
"""
function coupler_add_map!(coupler::CouplerState, map::Operators.LinearRemap)
    push!(coupler.remap_operators, (map.target, map.source) => map)
end

"""
    coupler_put!(coupler::CouplerState, fieldname::Symbol, fieldvalue)

Sets coupler field `fieldname` to `fieldvalue`.
"""
function coupler_put!(coupler::CouplerState, fieldname::Symbol, fieldvalue, source_sim::AbstractSimulation)
    cplfield = coupler.coupled_fields[fieldname]
    @assert cplfield.write_sim == name(source_sim) "$fieldname can only be written to by $(cplfield.write_sim)."

    map = coupler_get_map(coupler, axes(cplfield.data), axes(fieldvalue))
    Operators.remap!(cplfield.data, map, fieldvalue)

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
)
    cplfield = coupler.coupled_fields[fieldname]
    map = coupler_get_map(coupler, axes(target_field), axes(cplfield.data))
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

function coupler_get(coupler::CouplerState, fieldname::Symbol, target_space::ClimaCore.Spaces.AbstractSpace)
    cplfield = coupler.coupled_fields[fieldname]
    map = coupler_get_map(coupler, target_space, axes(cplfield.data))
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

function coupler_get_map(coupler, target_space::ClimaCore.Spaces.AbstractSpace, source_space::ClimaCore.Spaces.AbstractSpace)
    return coupler.remap_operators[target_space, source_space]
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
