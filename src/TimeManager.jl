#=
    TimeManager

This module facilitates calendar functions and temporal interpolations
of data.
=#

module TimeManager

using ..Utilities
using Dates

export current_date, strdate_to_datetime, datetime_to_strdate


"""
    current_date(cs::CoupledSimulation, t::Int)

Return the model date at the current timestep

# Arguments
- `cs`: [CoupledSimulation] containing info about the simulation
- `t`: [Int] number of seconds since simulation began
"""
current_date(cs::CoupledSimulation, t::Int) = cs.dates.date0[1] + Dates.Second(t)

"""
    strdate_to_datetime(strdate::String)

Convert from String ("YYYYMMDD") to Date format,
required by the official AMIP input files

# Arguments
- `strdate`: [String] to be converted to Date type 
"""
strdate_to_datetime(strdate::String) =
    Dates.DateTime(parse(Int, strdate[1:4]), parse(Int, strdate[5:6]), parse(Int, strdate[7:8]))

"""
    datetime_to_strdate(datetime::DateTime)

Convert from Date to String ("YYYYMMDD") format

# Arguments
- `datetime`: [DateTime] object to be converted to string
"""
datetime_to_strdate(datetime::DateTime) =
    string(lpad(Dates.year(datetime), 4, "0")) *
    string(string(lpad(Dates.month(datetime), 2, "0"))) *
    string(lpad(Dates.day(datetime), 2, "0"))

end
