# bcfile_reader
# - coordinates reading of boundary conditions from NetCDF files, as well as regridding calls and temporal interpolations from monthly to daily intervals

# IO - init
"""
    BCDataInfo

Stores information specific to each boundary condition from a file and each variable.
The inputs are:
    FT::F                           # float type
    datafile_cgll::S                # file containing all regridded fields
    varname::V                      # name of the variable
    weightfile::S                   # file containing regridding weights
    all_dates::D                    # all dates contained in the original data file
    segment_idx::Vector{Int}        # index of the monthly data in the file
    segment_idx0::Vector{Int}       # `segment_idx` of the file data that is closest to date0
    monthly_fields::C               # Tuple of the two monthly fields, that will be used for the daily interpolation
    segment_length::Vector{Int}     # length of each month segment (used in the daily interpolation)
    scaling_function::O             # function that scales, offsets or transforms the raw variable
    interpolate_daily::I            # switch to trigger daily interpolation
    land_mask::M                     # mask with 1 = land, 0 = ocean / sea-ice
"""
struct BCFileInfo{F, S, V, D, C, O, I, M}
    FT::F
    datafile_cgll::S
    varname::V
    weightfile::S
    all_dates::D
    segment_idx::Vector{Int}
    segment_idx0::Vector{Int}
    monthly_fields::C
    segment_length::Vector{Int}
    scaling_function::O
    interpolate_daily::I
    land_mask::M
end

"""
    Bcfile_info_init(datafile_rll, varname, boundary_space; interpolate_daily = false, segment_idx0 = [Int(1)], scaling_function = false)

Regrids from lat-lon grid to cgll grid, saving the output in a new file, and returns the info packaged in a single struct
"""
function bcfile_info_init(
    FT,
    datafile_rll,
    varname,
    boundary_space;
    interpolate_daily = false,
    segment_idx0 = nothing,
    scaling_function = no_scaling,
    land_mask = nothing,
    date0 = nothing,
)

    # regrid all times and save to file
    weightfile, datafile_cgll =
        ncreader_rll_to_cgll_from_space(datafile_rll, varname, boundary_space, outfile = varname * "_cgll.g")
    ds = Dataset(datafile_cgll, "r")

    # init time tracking info
    current_fields =
        interpolate_daily ? (Fields.zeros(FT, boundary_space), Fields.zeros(FT, boundary_space)) :
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

    # unless the start file date is specified, find the closest one to the start date
    segment_idx0 =
        segment_idx0 != nothing ? segment_idx0 :
        [argmin(abs.(parse(FT, datetime_to_strdate(date0)) .- parse.(FT, datetime_to_strdate.(data_dates[:]))))]

    return BCFileInfo(
        FT,
        datafile_cgll,
        varname,
        weightfile,
        data_dates,
        segment_idx0 .- Int(1),
        segment_idx0,
        current_fields,
        segment_length,
        scaling_function,
        [interpolate_daily],
        land_mask,
    )

end

no_scaling(x, _info) = swap_space!(x, axes(_info.land_mask))

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
    FT = bcf_info.FT

    if (all_dates == nothing) # temporally invariant BCs
        @warn "no temporally varying data, all months using the same field"
        map(
            x ->
                bcf_info.monthly_fields[x] .= scaling_function(
                    ncreader_cgll_sparse_to_field(
                        datafile_cgll,
                        varname,
                        weightfile,
                        (Int(midmonth_idx),),
                        boundary_space,
                    )[1],
                    bcf_info,
                ),
            Tuple(1:length(monthly_fields)),
        )
    elseif (midmonth_idx == midmonth_idx0) && (Dates.days(date - all_dates[midmonth_idx]) < 0) # for init
        @warn "this time period is before BC data - using file from $(all_dates[midmonth_idx0])"
        map(
            x ->
                bcf_info.monthly_fields[x] .= scaling_function(
                    ncreader_cgll_sparse_to_field(
                        datafile_cgll,
                        varname,
                        weightfile,
                        (Int(midmonth_idx0),),
                        boundary_space,
                    )[1],
                    bcf_info,
                ),
            Tuple(1:length(monthly_fields)),
        )
    elseif Dates.days(date - all_dates[end - 1]) > 0 # for fini
        @warn "this time period is after BC data - using file from $(all_dates[end - 1])"
        map(
            x ->
                bcf_info.monthly_fields[x] .= scaling_function(
                    ncreader_cgll_sparse_to_field(
                        datafile_cgll,
                        varname,
                        weightfile,
                        (Int(length(all_dates)),),
                        boundary_space,
                    )[1],
                    bcf_info,
                ),
            Tuple(1:length(monthly_fields)),
        )
    elseif Dates.days(date - all_dates[Int(midmonth_idx)]) > 2 # throw error when there are closer initial indices for the bc file data that matches this date0
        nearest_idx =
            argmin(abs.(parse(FT, datetime_to_strdate(date)) .- parse.(FT, datetime_to_strdate.(all_dates[:]))))
        @error "init data does not correspond to start date. Try initializing with `SIC_info.segment_idx = midmonth_idx = midmonth_idx0 = $nearest_idx` for this start date" # TODO: do this automatically w a warning
    elseif Dates.days(date - all_dates[Int(midmonth_idx - 1)]) > 0 # date crosses to the next month
        @warn "updating monthly data file"
        bcf_info.segment_length .= Dates.days(all_dates[Int(midmonth_idx + 1)] - all_dates[Int(midmonth_idx)])
        map(
            x ->
                bcf_info.monthly_fields[x] .= scaling_function(
                    ncreader_cgll_sparse_to_field(
                        datafile_cgll,
                        varname,
                        weightfile,
                        (Int(midmonth_idx), Int(midmonth_idx + 1)),
                        boundary_space,
                    )[x],
                    bcf_info,
                ),
            Tuple(1:length(monthly_fields)),
        )
    else
        @error "Check boundary file specification"
    end
end

next_date_in_file(bcfile_info) = bcfile_info.all_dates[bcfile_info.segment_idx[1] + Int(1)]
