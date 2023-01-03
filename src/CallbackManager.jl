#=
    CallbackManager - TODO WIP

This module facilitates calendar functions and temporal interpolations
of data.
=#

module CallbackManager

using ..Utilities
using Dates

export current_date, strdate_to_datetime, datetime_to_strdate, calendar_callback

# TODO remove calls to calendar_callback - can be replaced with `<`

"""
    current_date(cs, t)

Return the model date
"""
current_date(cs, t) = cs.dates.date0[1] + Dates.Second(t)

"""
    strdate_to_datetime(strdate)

Convert from String ("YYYYMMDD") to Date format  
"""
strdate_to_datetime(strdate::String) =
    Dates.DateTime(parse(Int, strdate[1:4]), parse(Int, strdate[5:6]), parse(Int, strdate[7:8])) # required by the official AMIP input files

"""
    datetime_to_strdate(datetime)

Convert from Date to String ("YYYYMMDD") format  
"""
datetime_to_strdate(datetime::DateTime) =
    string(Dates.year(datetime)) *
    string(string(lpad(Dates.month(datetime), 2, "0"))) *
    string(lpad(Dates.day(datetime), 2, "0"))


end
