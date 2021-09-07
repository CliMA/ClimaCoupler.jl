using Unitful, Dates
using PrettyTables

export CouplerState
export coupler_push!, coupler_pull!, coupler_put!, coupler_get
export coupler_add_field!

# TODO: Build constructor that uses a model component's grid.
struct CplGridInfo{GT, GP, GH, GE}
    gridtype::GT
    gridparams::GP
    gridhandle::GH
    gridencoding::GE
end

mutable struct CplFieldInfo{AT, GT <: CplGridInfo, UT <: Unitful.Units}
    data::AT
    gridinfo::GT
    units::UT
    datetime::DateTime
end

mutable struct CouplerState{DT}
    CplStateBlob::DT
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
        grid,
        datetime::DateTime,
        units::Unitful.Units = Unitful.NoUnits, 
    )
    gridinfo = CplGridInfo(nothing, nothing, nothing, nothing)
    push!(coupler.CplStateBlob, fieldname => CplFieldInfo(fieldvalue, gridinfo, units, datetime))
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
function coupler_get(coupler::CouplerState, fieldname::Symbol, gridinfo, datetime::DateTime, units::Unitful.Units)
    cplfield = coupler.CplStateBlob[fieldname]

    # check that retrieving component and coupler are at same time
    datetime != cplfield.datetime && 
        throw(ErrorException("Retrieval time ($datetime) != coupler field time ($(cplfield.datetime))"))

    regriddata = regrid(cplfield.data, gridinfo, cplfield.gridinfo)

    return uconvert(regriddata, units, cplfield.units)
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
function coupler_put!(coupler::CouplerState, fieldname::Symbol, fieldvalue, gridinfo, datetime::DateTime, units::Unitful.Units)
    cplfield = coupler.CplStateBlob[fieldname]

    # map new data to grid of coupler field; new data -> coupler grid
    regriddata = regrid(fieldvalue, cplfield.gridinfo, gridinfo)

    # store new regridded data in correct units
    cplfield.data .= uconvert(regriddata, cplfield.units, units)

    # update timestamp of coupler field
    settime!(cplfield, datetime)

    return nothing
end

# TODO: maps cplfield data from `fromgrid` to `togrid`
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
    fields = coupler.CplStateBlob
    data = Array{Any}(undef, length(fields), 4)
    for (i,k) in enumerate( keys(coupler.CplStateBlob))
        entry = coupler.CplStateBlob[k]
        data[i,:] = [k entry.units entry.datetime entry.gridinfo.gridtype]
    end
    pretty_table(data, ["Field Name" "Units" "Field Time (Y-m-d H:M:S)" "Grid Type"])
end