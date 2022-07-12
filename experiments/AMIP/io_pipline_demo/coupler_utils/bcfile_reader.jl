# bcfile_reader
# - coordinates reading of boundary conditions from NetCDF files, as well as regridding calls and temporal interpolations from monthly to daily intervals

# IO - init
"""
BCDataInfo
- specific to each boundary condition from a file, per variable
"""
struct BCFileInfo{S,V,D,C,O,I}
    datafile_cgll::S  
    varname::V
    weightfile::S
    all_dates::D
    segment_idx::Vector{Int}
    segment_idx0::Vector{Int}
    monthly_fields::C
    segment_length::Vector{Int}
    scaling_function::O # TODO: turn this into a func for scaling + offsets etc
    interpolate_monthly::I
end

"""
bcfile_info_init
- regrids from lat-lon grid to cgll grid, saving the output in a new file, and returns info relevant for temporal interpolation 
"""
function bcfile_info_init(datafile_rll, varname, boundary_space; interpolate_monthly = false, segment_idx0 = [Int(1)], scaling_function = false)

    # regrid all times and save to file
    weightfile, datafile_cgll = ncreader_rll_to_cgll_from_space(FT, datafile_rll, varname, boundary_space, outfile = varname*"_cgll.g")
    ds = Dataset(datafile_cgll, "r")

    # init time tracking info
    current_fields = interpolate_monthly ? (Fields.zeros(FT, boundary_space), Fields.zeros(FT, boundary_space)) : (Fields.zeros(FT, boundary_space),)
    segment_length = [Int(0)]

    if interpolate_monthly == true
        return BCFileInfo(datafile_cgll, varname, weightfile, Dates.Date.(ds["time"][:]), segment_idx0 .- Int(1), segment_idx0, current_fields, segment_length, scaling_function,[interpolate_monthly]) # TODO: generalize to DateTime (only dialy resol now)
    else
        return BCFileInfo(datafile_cgll, varname, weightfile, nothing, segment_idx0 .- Int(1), segment_idx0, current_fields, segment_length, scaling_function, [interpolate_monthly]) 
    end
end
# strdate_to_datetime.(string.( ds["time"][:])) # if string "date"


# IO - monthly
"""
update_midmonth_data!(bcf_info)
- 
"""
function update_midmonth_data!(date, bcf_info)

    # monthly count
    bcf_info.segment_idx[1] += Int(1)

    all_dates = bcf_info.all_dates
    midmonth_idx = bcf_info.segment_idx[1]
    midmonth_idx0 = bcf_info.segment_idx0[1]
    datafile_cgll = bcf_info.datafile_cgll
    weightfile = bcf_info.weightfile
    monthly_fields = bcf_info.monthly_fields
    scaling_function = bcf_info.scaling_function
    varname = bcf_info.varname

    if (all_dates == nothing)
        @warn "no temporally varying data, all months using the same field"
        bcf_info.monthly_fields[1] .= bcfields_from_file(datafile_cgll, varname, weightfile, (Int(midmonth_idx),), axes(monthly_fields[1]), scaling_function = scaling_function)[1]
        bcf_info.interpolate_monthly .= false
    elseif (midmonth_idx == midmonth_idx0) && (Dates.days(date - all_dates[midmonth_idx])  < 0)
        @warn "no BC data for this time period - using file from $(all_dates[1])"
        bcf_info.monthly_fields[1] .= bcfields_from_file(datafile_cgll, varname, weightfile, (Int(midmonth_idx),), axes(monthly_fields[1]), scaling_function = scaling_function)[1]
        bcf_info.interpolate_monthly .= false
    elseif Dates.days(date - all_dates[end - 1]) > 0
        @warn "no BC data for this time period - using file from $(all_dates[end - 1])"
        bcf_info. monthly_fields[1] .= bcfields_from_file(datafile_cgll, varname, weightfile, (Int(length(all_dates)),), axes(monthly_fields[1]), scaling_function = scaling_function)[1]  
        bcf_info.interpolate_monthly .= false
    elseif Dates.days(date - all_dates[Int(midmonth_idx + 1)]) > 20
        nearest_idx = argmin(abs.(parse(FT,datetime_to_strdate(date)) .- parse.(FT, datetime_to_strdate.(all_dates[:]))))
        @error "init data does not correspond to start date. Try initializing with `SIC_info.segment_idx = midmonth_idx = midmonth_idx0 = $nearest_idx` for this start date" # TODO: do this automatically w a warning
    elseif Dates.days(date - all_dates[Int(midmonth_idx)]) > 0
        bcf_info.segment_length .= Dates.days(all_dates[Int(midmonth_idx + 1)] - all_dates[Int(midmonth_idx)])
        map(x -> bcf_info.monthly_fields[x] .= bcfields_from_file(datafile_cgll, varname, weightfile, (Int(midmonth_idx), Int(midmonth_idx + 1)), axes(monthly_fields[1]), scaling_function = scaling_function)[x], (1,2)) 
        bcf_info.interpolate_monthly .= true
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
    if bcf_info.interpolate_monthly[1]
        month_fraction = FT(Dates.days(date - bcf_info.all_dates[Int(bcf_info.segment_idx[1])]) / bcf_info.segment_length[1])
        @assert abs(month_fraction) <= FT(1) "time interpolation weights must be <= 1, but month_fraction = $month_fraction"
        return intepol.(bcf_info.monthly_fields[1], bcf_info.monthly_fields[2], month_fraction, FT) 
    else
        return bcf_info.monthly_fields 
    end 
end

intepol(ftuple1, ftuple2, month_fraction, FT) = ftuple1 * month_fraction + ftuple2 * (FT(1) - month_fraction)

# TODO:
# - Dates.epochdays2date
# - generalize for finer temp resol
