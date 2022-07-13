# bcfile_reader
# - coordinates reading of boundary conditions from NetCDF files, as well as regridding calls and temporal interpolations from monthly to daily intervals

# IO - init
"""
    BCDataInfo

Stores information specific to each boundary condition from a file and each variable.
"""
struct BCFileInfo{S, V, R, D, C, O, I}
    datafile_cgll::S
    varname::V
    weightfile::S
    regrid_space::R
    all_dates::D
    segment_idx::Vector{Int}
    segment_idx0::Vector{Int}
    monthly_fields::C
    segment_length::Vector{Int}
    scaling_function::O # TODO: turn this into a func for scaling + offsets etc
    interpolate_monthly::I
end

"""
    Bcfile_info_init(datafile_rll, varname, boundary_space; interpolate_monthly = false, segment_idx0 = [Int(1)], scaling_function = false)

Regrids from lat-lon grid to cgll grid, saving the output in a new file, and returns the info packaged in a single struct
"""
function bcfile_info_init(
    datafile_rll,
    varname,
    boundary_space;
    interpolate_monthly = false,
    segment_idx0 = [Int(1)],
    scaling_function = false,
)

    # regrid all times and save to file
    weightfile, datafile_cgll, regrid_space =
        ncreader_rll_to_cgll_from_space(datafile_rll, varname, boundary_space, outfile = varname * "_cgll.g")
    ds = Dataset(datafile_cgll, "r")

    # init time tracking info
    current_fields =
        interpolate_monthly ? (Fields.zeros(FT, boundary_space), Fields.zeros(FT, boundary_space)) :
        (Fields.zeros(FT, boundary_space),)
    segment_length = [Int(0)]

    if "time" in ds
        data_dates = Dates.DateTime.(ds["time"][:])
    elseif "date" in ds
        data_dates = strdate_to_datetime.(string.(ds["date"][:]))
    else
        @warn "No dates availabe in file $datafile_rll"
        data_dates = nothing
    end

    return BCFileInfo(
        datafile_cgll,
        varname,
        weightfile,
        regrid_space,
        data_dates,
        segment_idx0 .- Int(1),
        segment_idx0,
        current_fields,
        segment_length,
        scaling_function,
        [interpolate_monthly],
    ) # TODO: generalize to DateTime (only dialy resol now)

end

# IO - monthly
"""
    update_midmonth_data!(bcf_info)

Extracts boundary condition data from regridded (to model grid) NetCDF files (which times, depends on the specifications in the `bcf_info` struct).
"""
function update_midmonth_data!(date, bcf_info)

    # monthly count
    bcf_info.segment_idx[1] += Int(1)

    all_dates = bcf_info.all_dates
    midmonth_idx = bcf_info.segment_idx[1]
    midmonth_idx0 = bcf_info.segment_idx0[1]
    monthly_fields = bcf_info.monthly_fields
    datafile_cgll = bcf_info.datafile_cgll
    weightfile = bcf_info.weightfile
    scaling_function = bcf_info.scaling_function
    varname = bcf_info.varname
    regrid_space = bcf_info.regrid_space

    if (all_dates == nothing) # temporally invariant BCs
        @warn "no temporally varying data, all months using the same field"
        map(
            x ->
                bcf_info.monthly_fields[x] .= ncreader_cgll_sparse_to_field(
                    datafile_cgll,
                    varname,
                    weightfile,
                    (Int(midmonth_idx),),
                    regrid_space,
                    scaling_function = scaling_function,
                )[1],
            Tuple(1:length(monthly_fields)),
        )
    elseif (midmonth_idx == midmonth_idx0) && (Dates.days(date - all_dates[midmonth_idx]) < 0) # for init
        @warn "this time period is before BC data - using file from $(all_dates[midmonth_idx0])"
        map(
            x ->
                bcf_info.monthly_fields[x] .= ncreader_cgll_sparse_to_field(
                    datafile_cgll,
                    varname,
                    weightfile,
                    (Int(midmonth_idx0),),
                    regrid_space,
                    scaling_function = scaling_function,
                )[1],
            Tuple(1:length(monthly_fields)),
        )
    elseif Dates.days(date - all_dates[end - 1]) > 0 # for fini
        @warn "this time period is after BC data - using file from $(all_dates[end - 1])"
        map(
            x ->
                bcf_info.monthly_fields[x] .= ncreader_cgll_sparse_to_field(
                    datafile_cgll,
                    varname,
                    weightfile,
                    (Int(length(all_dates)),),
                    regrid_space,
                    scaling_function = scaling_function,
                )[1],
            Tuple(1:length(monthly_fields)),
        )
    elseif Dates.days(date - all_dates[Int(midmonth_idx)]) > 2 # throw error when there are closer initial indices for the bc file data that matches this date0
        nearest_idx =
            argmin(abs.(parse(FT, datetime_to_strdate(date)) .- parse.(FT, datetime_to_strdate.(all_dates[:]))))
        @error "init data does not correspond to start date. Try initializing with `SIC_info.segment_idx = midmonth_idx = midmonth_idx0 = $nearest_idx` for this start date" # TODO: do this automatically w a warning
    elseif Dates.days(date - all_dates[Int(midmonth_idx - 1)]) > 0 # date crosses to the next month
        bcf_info.segment_length .= Dates.days(all_dates[Int(midmonth_idx + 1)] - all_dates[Int(midmonth_idx)])
        map(
            x ->
                bcf_info.monthly_fields[x] .= ncreader_cgll_sparse_to_field(
                    datafile_cgll,
                    varname,
                    weightfile,
                    (Int(midmonth_idx), Int(midmonth_idx + 1)),
                    regrid_space,
                    scaling_function = scaling_function,
                )[x],
            Tuple(1:length(monthly_fields)),
        )
    else
        @error "Check boundary file specification"
    end
end

next_month_date(bcfile_info) = bcfile_info.all_dates[bcfile_info.segment_idx[1] + Int(1)]

# IO - daily
"""
    interpolate_midmonth_to_daily(date, bcf_info)

Interpolates linearly between two `Fields` in the `bcf_info` struct, or returns the first Field if interpolation is switched off. 
"""
function interpolate_midmonth_to_daily(date, bcf_info)
    if bcf_info.interpolate_monthly[1]
        month_fraction = FT(Dates.days(date - bcf_info.all_dates[Int(bcf_info.segment_idx[1])]) / bcf_info.segment_length[1])
        @assert abs(month_fraction) <= FT(1) "time interpolation weights must be <= 1, but month_fraction = $month_fraction"
        return intepol.(bcf_info.monthly_fields[1], bcf_info.monthly_fields[2], month_fraction, FT) 
    else
        return bcf_info.monthly_fields[1] 
    end 
end

intepol(ftuple1, ftuple2, month_fraction, FT) = ftuple1 * month_fraction + ftuple2 * (FT(1) - month_fraction)

# TODO:
# - design unit test for `update_midmonth_data!` and `interpolate_midmonth_to_daily`
# - replace if statements with dipatches, write better abstractions
