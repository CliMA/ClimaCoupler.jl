"""
    TimeManager

This module facilitates calendar functions and temporal interpolations
of data.
"""
module TimeManager

using ..Interfacer
using Dates

using DocStringExtensions
const DSE = DocStringExtensions

export current_date,
    strdate_to_datetime,
    datetime_to_strdate,
    trigger_callback,
    AbstractFrequency,
    Monthly,
    EveryTimestep,
    trigger_callback!,
    HourlyCallback,
    MonthlyCallback,
    update_firstdayofmonth!


"""
    current_date(cs::CoupledSimulation, t::Int)

Return the model date at the current timestep.

# Arguments
- `cs`: [CoupledSimulation] containing info about the simulation
- `t`: [Real] number of seconds since simulation began
"""
current_date(cs::Interfacer.CoupledSimulation, t::Real) = cs.dates.date0[1] + Dates.Second(t)

"""
    strdate_to_datetime(strdate::String)

Convert from String ("YYYYMMDD") to Date format,
required by the official AMIP input files.

# Arguments
- `strdate`: [String] to be converted to Date type
"""
strdate_to_datetime(strdate::String) =
    Dates.DateTime(parse(Int, strdate[1:4]), parse(Int, strdate[5:6]), parse(Int, strdate[7:8]))

"""
    datetime_to_strdate(datetime::DateTime)

Convert from Date to String ("YYYYMMDD") format.

# Arguments
- `datetime`: [DateTime] object to be converted to string
"""
datetime_to_strdate(datetime::DateTime) =
    string(lpad(Dates.year(datetime), 4, "0")) *
    string(string(lpad(Dates.month(datetime), 2, "0"))) *
    string(lpad(Dates.day(datetime), 2, "0"))

abstract type AbstractFrequency end
struct Monthly <: AbstractFrequency end
struct EveryTimestep <: AbstractFrequency end

"""
    trigger_callback(cs, ::Monthly)

Returns `true` if the current date is equal to or exceeds the saved first of the month at time of 00:00:00.

# Arguments
- `cs`: [CoupledSimulation] containing info about the simulation
"""
trigger_callback(cs::CoupledSimulation, ::Monthly) = cs.dates.date[1] >= cs.dates.date1[1] ? true : false


abstract type CouplerCallback end

"""
    do_nothing(::Interfacer.CoupledSimulation, _)

This is a helper callback function that does nothing.
"""
do_nothing(::Interfacer.CoupledSimulation, _) = nothing

"""
    HourlyCallback{FT}

This is a callback type that triggers at intervals of 1h or multiple hours.

# Fields

$(DSE.FIELDS)
"""
@kwdef struct HourlyCallback{FT} <: CouplerCallback
    dt::FT = FT(1) # hours
    func::Function = do_nothing
    ref_date::Array = [Dates.DateTime(0)]
    active::Bool = false
    data::Array = []
end

"""
    MonthlyCallback{FT}

This is a callback type that triggers at intervals of 1 month or multiple months.

# Fields

$(DSE.FIELDS)
"""
@kwdef struct MonthlyCallback{FT} <: CouplerCallback
    dt::FT = FT(1) # months
    func::Function = do_nothing
    ref_date::Array = [Dates.DateTime(0)]
    active::Bool = false
    data::Array = []
end

dt_cb(cb::HourlyCallback) = Dates.Hour(cb.dt)
dt_cb(cb::MonthlyCallback) = Dates.Month(cb.dt)

"""
    trigger_callback!(cs::Interfacer.CoupledSimulation, cb::CouplerCallback)

This function triggers a callback function if the current date is equal to or exceeds the saved callback reference date.
As well as executing the functions `func`, it automatically updates the reference date, `ref_date`, for the next callback interval.
"""
function trigger_callback!(cs::Interfacer.CoupledSimulation, cb::CouplerCallback)
    if cb.active
        current_date = cs.dates.date[1]
        if current_date >= cb.ref_date[1]
            cb.func(cs, cb)
            cb.ref_date[1] = cb.ref_date[1] + dt_cb(cb)
        end
    end
end

"""
    update_firstdayofmonth!(cs::Interfacer.CoupledSimulation, _)

This function updates the first of the month reference date.
"""
function update_firstdayofmonth!(cs, _)
    cs.dates.date1[1] = cs.dates.date1[1] + Dates.Month(1)
    @info("update_firstdayofmonth! at $(cs.dates.date)")
end

end
