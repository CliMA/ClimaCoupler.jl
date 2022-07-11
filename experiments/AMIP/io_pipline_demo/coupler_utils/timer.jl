# time manager
#  - fascilitates calendar functions and temporal interpolations

"""
current_day(date0, t, FT)
- returns the date and simulation day index (NB: only dialy resolution with Dates)
"""

function current_date(date0, t, FT)
    day_secs = FT(60 * 60 * 24)
    sim_day = (t - mod(t, day_secs )) / day_secs

    date = date0 + Dates.Day(sim_day)
    return date
end

strdate_to_datetime(strdate) = Dates.Date(parse(Int,strdate[1:4]), parse(Int,strdate[5:6]), parse(Int,strdate[7:8]))

datetime_to_strdate(datetime) = string(Dates.year(datetime))*string(Dates.month(datetime))*string(Dates.day(datetime))

