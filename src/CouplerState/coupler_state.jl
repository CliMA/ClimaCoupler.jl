using Unitful, Dates
using PrettyTables
using ClimaCore

export CouplerState
export coupler_push!, coupler_pull!, coupler_put!, coupler_get
export coupler_add_field!

mutable struct CplFieldInfo{AT}
    data::AT
end

mutable struct CouplerState{DT}
    CplStateDict::DT
end

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
            grid,
            datetime::DateTime,
            units::Unitful.Units = Unitful.NoUnits, 
        )

Add a field to the coupler that is accessible with key `fieldname`. 

# Arguments
- `coupler`: coupler object the field is added to.
- `fieldname`: key to access the field in the coupler.
- `fieldvalue`: data array of field values.
- `grid`: grid the field is stored on.
- `datetime`: time associated with the field state.
- `units`: units associated with the field values. Dimensionless by default.
"""
function coupler_add_field!(
    coupler::CouplerState,
    fieldname::Symbol,
    fieldvalue,
    # grid = nothing, # by default, the coupler grid is the grid of the registering model. Can otherwise specify the CC.space the field should live on
)
    push!(coupler.CplStateDict, fieldname => CplFieldInfo(fieldvalue))
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

"""
    coupler_get(coupler::CouplerState, fieldname::Symbol, gridinfo, datetime::DateTime, units::Unitful.Units)

Retrieve data array corresponding to `fieldname`.

Returns data on the grid specified by `gridinfo` and in the units of `units`. Checks that
the coupler data field is the state at time `datetime`.
"""
# function coupler_get(coupler::CouplerState, fieldname::Symbol, gridinfo, datetime::DateTime, units::Unitful.Units)
#     cplfield = coupler.CplStateDict[fieldname]

#     # check that retrieving component and coupler are at same time
#     datetime != cplfield.datetime &&
#         throw(ErrorException("Retrieval time ($datetime) != coupler field time ($(cplfield.datetime))"))

#     regriddata = regrid(cplfield.data, gridinfo, cplfield.gridinfo)

#     return uconvert(regriddata, units, cplfield.units)
# end
function coupler_get(coupler::CouplerState, fieldname::Symbol, regrid_space = nothing)
    cplfield = coupler.CplStateDict[fieldname]
    
    # don't regrid by default
    regrid_space = (regrid_space === nothing) ? axes(cplfield.data) : regrid_space

    # call to climacore remap utils
    regriddata = regrid(cplfield.data, regrid_space, axes(cplfield.data))

    return cplfield.data
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
    coupler_put!(coupler::CouplerState, fieldname::Symbol, fieldvalue, gridinfo, datetime::DateTime, units::Unitful.Units)

Updates coupler field `fieldname` with `fieldvalue`, the field's value at time `datetime`.

`gridinfo` and `units` inform the coupler of the format of the inputted data
allowing conversion to match the grid and units of the coupler field.
"""
function coupler_put!(
    coupler::CouplerState,
    fieldname::Symbol,
    fieldvalue,
)
    cplfield = coupler.CplStateDict[fieldname]

    # map new data to grid of coupler field; new data -> coupler grid
    regriddata = regrid(fieldvalue, axes(cplfield.data), axes(fieldvalue))

    cplfield.data .= regriddata

    return nothing
end

# TODO: maps cplfield data from `fromgrid` to `togrid`
# Would like a way to use the togrid and fromgrid to access the right remapping operator
# We only really (currently, and is cleanest) have space info via axes(fields) to use here.
# perhaps, a put regrids from :component to :coupler while a get regrids from :coupler to component.
# Operators need ordered keys to be accessed?
# Operators are auto-constructed during registration phase?
function regrid(data, togrid, fromgrid)
    return data
end

# convert from `fromunits` to `tounits`
# ustrip throws error if the units aren't dimensionally compatible
function uconvert(data, tounits::Unitful.Units, fromunits::Unitful.Units)
    tounits == fromunits && return data
    return ustrip.(tounits, data .* fromunits)
end

# set the time of the coupled field
function settime!(cplfield::CplFieldInfo, newtime::DateTime)
    cplfield.datetime = newtime
end

# display table of registered fields & their info
function Base.show(io::IO, coupler::CouplerState)
    #=
    ---------------------------------------------------------------------------
    |  Field Name  |  Units  |  Field Time (Y-m-d H:M:S)   |     Grid Type    |
    ---------------------------------------------------------------------------
    |  :OceanSST   |   °C    |  2021-01-01T00:00:00        |  DG Cubed Sphere |
    |  :Field2     |    m    |  2021-01-01T00:00:00        |  CG Cubed Sphere |
    ---------------------------------------------------------------------------
    =#
    fields = coupler.CplStateDict
    data = Array{Any}(undef, length(fields), 4)
    for (i, k) in enumerate(keys(coupler.CplStateDict))
        entry = coupler.CplStateDict[k]
        data[i, :] = [k entry.units entry.datetime entry.gridinfo.gridtype]
    end
    header = (["Field Name", "Units", "Field Time", "Grid Type"], ["", "", "Y-m-d H:M:S", ""])
    pretty_table(data, header = header)
end
