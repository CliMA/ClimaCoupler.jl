# calendar timer
#  - fascilitates calendar functions and temporal interpolations

"""
    current_date(t)
- returns the model date
"""
current_date(t) = date0 + Dates.Second(t)

"""
    strdate_to_datetime(strdate), datetime_to_strdate(datetime)
- convert between String ("YYYYMMDD") and Date formats  
"""
strdate_to_datetime(strdate::String) = Dates.Date(parse(Int,strdate[1:4]), parse(Int,strdate[5:6]), parse(Int,strdate[7:8])) # required by the official AMIP input files

datetime_to_strdate(datetime::Date) = string(Dates.year(datetime))*string(Dates.month(datetime))*string(Dates.day(datetime))

"""
    calendar_callback(ex, model_date, callback_date)
- evaluate `ex` at `callback_date` and do nothing otherwise
"""
macro calendar_callback(ex::Expr, model_date::Symbol, callback_date::Symbol)
    quote
        if Dates.days($model_date - $callback_date) < FT(0)
            nothing
        else 
            eval($ex)
        end
    end
end

# # test for @calendar_callback (TODO: modify and use when move to `src/`)
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
