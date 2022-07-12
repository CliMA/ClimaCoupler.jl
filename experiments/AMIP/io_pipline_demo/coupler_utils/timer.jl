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

# # test
# SST_info = bcfile_info_init(sst_data, "SST", boundary_space, segment_idx0 = [Int(1729)],interpolate_monthly = true, scaling_function = clean_sst)
# update_midmonth_data!(date0, SST_info)

# for t in ((tspan[1] + Δt_cpl):Δt_cpl:tspan[end])
#     #@show t
#     date = current_date(date0, t, FT)
#     next_monnth = SST_info.all_dates[SST_info.segment_idx[1] + Int(1)]
#     Dates.days(date - next_monnth) < FT(1) ? println(Dates.days(date - next_monnth) ) : println("yes:$(date) vs $(next_monnth)")  
# end
