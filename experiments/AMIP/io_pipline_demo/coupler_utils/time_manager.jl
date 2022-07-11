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

# IO - init
"""
BCDataInfo
- specific to each boundary condition from a file, per variable
"""
struct BCFileInf{S,D,M,C,L,O}
    datafile_cgll::S  
    weightfile::S
    all_dates::D
    segment_idx::Vector{Int}
    monthly_fields::C
    segment_length::Vector{Int}
    value_offset::O # TODO: turn this into a func for scaling + offsets etc
end

"""
bcfile_info_init
- regrids from lat-lon grid to cgll grid, saving the output in a new file, and returns info relevant for temporal interpolation 
"""
function bcfile_info_init(datafile_rll, varname, boundary_space, value_offset = FT(0))

    # regrid all times and save to file
    weightfile, datafile_cgll = ncreader_rll_to_cgll_from_space(FT, datafile_rll, varname, boundary_space, outfile = varname*"_cgll.g")
    ds = Dataset(datafile_cgll, "r")

    # init time tracking info
    current_fields = (Fields.zeros(FT, boundary_space), Fields.zeros(FT, boundary_space))
    segment_length = [Int(0)]
    segment_idx = [Int(1)]

    return BCFileInfo(datafile_cgll, weightfile, Dates.Date.(ds["time"][:]), segment_idx, current_fields, segment_length, value_offset) # TODO: generalize to DateTime (only dialy resol now)
end
# strdate_to_datetime.(string.( ds["time"][:])) # if string "date"


# IO - monthly
"""
update_midmonth_data!(bcf_info, midmonth_idx0 = Int(1))
- 
"""
function update_midmonth_data!(bcf_info, midmonth_idx0 = Int(1))
    
    all_dates = bcf_info.all_dates
    midmonth_idx = bcf_info.segment_idx
    datafile_cgll = bcf_info.datafile_cgll
    weightfile = bcf_info.weightfile
    monthly_fields = bcf_info.monthly_fields

    if (midmonth_idx == midmonth_idx0) && (Dates.days(date - all_dates[midmonth_idx])  < 0)
        @warn "no BC data for this time period - using file from $(all_dates[1])"
        bcf_info.monthly_fields[1] .= bcfields_from_file(datafile_cgll, weightfile, (Int(midmonth_idx),), axes(monthly_fields[1]))[1]
    elseif Dates.days(date - all_dates[end - 1]) > 0
        @warn "no BC data for this time period - using file from $(all_dates[end - 1])"
        bcf_info. monthly_fields[1] .= bcfields_from_file(datafile_cgll, weightfile, (Int(length(all_dates)),), axes(monthly_fields[1]))[1]  
    elseif Dates.days(date - all_dates[Int(midmonth_idx[1] + 1)]) > 20
        nearest_idx = argmin(abs.(parse(FT,datetime_to_strdate(date)) .- parse.(FT, datetime_to_strdate.(all_dates[:]))))
        @error "init data does not correspond to start date. Try initializing with `SIC_info.segment_idx = midmonth_idx = midmonth_idx0 = $nearest_idx` for this start date" # TODO: do this automatically w a warning
    elseif Dates.days(date - all_dates[Int(midmonth_idx[1])]) > 0
        bcf_info.segment_length .= Dates.days(all_dates[Int(midmonth_idx[1] + 1)] - all_dates[Int(midmonth_idx[1])])
        map(x -> bcf_info.monthly_fields[x] .= bcfields_from_file(datafile_cgll, weightfile, (Int(midmonth_idx[1]), Int(midmonth_idx[1] + 1)), axes(monthly_fields[1]))[x], (1,2)) 
    else
        nothing
    end
end

# IO - daily
"""
interpolate_midmonth_to_daily(date, bcf_info)
- 
"""
function interpolate_midmonth_to_daily(date, bcf_info)
    month_fraction = Dates.days(date - bcf_info.all_dates[Int(bcf_info.segment_idx[1])]) / bcf_info.segment_length[1]
    length(bcf_info.monthly_fields) > 1 ? intepol.(bcf_info.monthly_fields[1], bcf_info.monthly_fields[2], month_fraction, FT) .+  bcf_info.value_offset : bcf_info.monthly_fields .+  bcf_info.value_offset
end

intepol(ftuple1, ftuple2, month_fraction, FT) = ftuple1 * month_fraction + ftuple2 * (FT(1) - month_fraction)

# TODO:
# - Dates.epochdays2date
# - generalize for finer temp resol
