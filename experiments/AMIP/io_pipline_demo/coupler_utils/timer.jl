# time manager
#  - fascilitates calendar functions and temporal interpolations

"""
    current_day(date0, t, FT)
- returns the date and simulation day index (NB: only dialy resolution with Dates)
"""

function current_date(date0, t, FT)
    day_secs = FT(86400)
    sim_day = (t - mod(t, day_secs )) / day_secs

    date = date0 + Dates.Day(sim_day)
    return date
end

strdate_to_datetime(strdate) = Dates.Date(parse(Int,strdate[1:4]), parse(Int,strdate[5:6]), parse(Int,strdate[7:8])) # required by the official AMIP input files

datetime_to_strdate(datetime) = string(Dates.year(datetime))*string(Dates.month(datetime))*string(Dates.day(datetime))

calendar_callback(exp, curr_date, callback_date) = Dates.days(curr_date - callback_date) < FT(0) ? nothing : eval(exp)


# # basic test 

# Δt_cpl = 0.5 * 86400
# tspan = (0,60* 86400)
# date0 = strdate_to_datetime("19000101")
# midmonth_dates = map(x -> date0 + Dates.Day(x), cumsum((15,31,29)))
# segment_idx = 1
# for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])
#     date = current_date(date0, t, FT)
#     @show date
#     next_month = midmonth_dates[segment_idx]
#     @show next_month
#     Dates.days(date - next_month) < FT(0) ? nothing : (println(Dates.days(date - next_month) ) , segment_idx +=1 )
# end

macro calendar_callback(ex, curr_date, callback_date)
    quote
        if Dates.days($curr_date - $callback_date) < FT(0)
            nothing
        else 
            eval($ex)
        end
    end
end

# # test for @calendar_callback
# Δt_cpl = 0.5 * 86400
# tspan = (0,60* 86400)
# date0 = strdate_to_datetime("19000101")
# midmonth_dates = map(x -> date0 + Dates.Day(x), cumsum((15,31,29)))
# segment_idx = 1
# for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])
#     date = current_date(date0, t, FT)
#     @show date
#     next_month = midmonth_dates[segment_idx]
#     @show next_month
#     #calendar_callback( :(println(Dates.days(esc(curr_date) - esc(callback_date)) ) , segment_idx +=1 ), date, next_month)
#     global expression = :(println(Dates.days(date - next_month) )  , segment_idx +=1 )
#     @calendar_callback( expression, date, next_month)
# end
