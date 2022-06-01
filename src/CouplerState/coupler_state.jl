using PrettyTables

export CouplerState
export coupler_push!, coupler_pull!, coupler_put!, coupler_get
export coupler_add_field!

mutable struct CplFieldInfo{AT, MD}
    data::AT
    metadata::MD
end

mutable struct CouplerState{DT}
    coupled_fields::DT
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
    return CouplerState(Dict{Symbol, CplFieldInfo}())
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
function coupler_add_field!(coupler::CouplerState, fieldname::Symbol, fieldvalue, metadata = nothing)
    push!(coupler.coupled_fields, fieldname => CplFieldInfo(fieldvalue, metadata))
end

"""
    coupler_put!(coupler::CouplerState, fieldname::Symbol, fieldvalue)

Sets coupler field `fieldname` to `fieldvalue`.
"""
function coupler_put!(coupler::CouplerState, fieldname::Symbol, fieldvalue)
    cplfield = coupler.coupled_fields[fieldname]

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
    coupler_get(coupler::CouplerState, fieldname::Symbol)

Retrieve data array corresponding to `fieldname`.

Returns data on the grid specified by `gridinfo` and in the units of `units`. Checks that
the coupler data field is the state at time `datetime`.
"""
function coupler_get(coupler::CouplerState, fieldname::Symbol)
    cplfield = coupler.coupled_fields[fieldname]

    return cplfield.data
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
