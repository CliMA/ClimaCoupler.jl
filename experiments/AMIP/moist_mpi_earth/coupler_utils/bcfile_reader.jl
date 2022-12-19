# bcfile_reader
# - coordinates reading of boundary conditions from NetCDF files, as well as regridding calls and temporal interpolations from monthly to daily intervals

# IO - init
"""
    BCDataInfo

Stores information specific to each boundary condition from a file and each variable.
The inputs are:
    comms_ctx::X                    # communications context used for MPI
    datafile_cgll::S                # file containing all regridded fields
    varname::V                      # name of the variable
    all_dates::D                    # all dates contained in the original data file
    monthly_fields::C               # Tuple of the two monthly fields, that will be used for the daily interpolation
    scaling_function::O             # function that scales, offsets or transforms the raw variable
    land_mask::M                    # mask with 1 = land, 0 = ocean / sea-ice
    segment_idx::Vector{Int}        # index of the monthly data in the file
    segment_idx0::Vector{Int}       # `segment_idx` of the file data that is closest to date0
    segment_length::Vector{Int}     # length of each month segment (used in the daily interpolation)
    interpolate_daily::Bool         # switch to trigger daily interpolation

"""
struct BCFileInfo{FT, X, S, V, D, C, O, M}
    comms_ctx::X
    hd_outfile_root::S
    varname::V
    all_dates::D
    monthly_fields::C
    scaling_function::O
    land_mask::M
    segment_idx::Vector{Int}
    segment_idx0::Vector{Int}
    segment_length::Vector{Int}
    interpolate_daily::Bool
end

BCFileInfo{FT}(args...) = BCFileInfo{FT, typeof.(args[1:7])...}(args...)

float_type(::BCFileInfo{FT}) where {FT} = FT

"""
    bcfile_info_init(FT, comms_ctx, datafile_rll, varname, boundary_space; interpolate_daily = false, segment_idx0 = [Int(1)], scaling_function = false)

Regrids from lat-lon grid to cgll grid, saving the output in a new file, and returns the info packaged in a single struct
"""
function bcfile_info_init(
    FT,
    comms_ctx,
    datafile_rll,
    varname,
    boundary_space;
    interpolate_daily = false,
    segment_idx0 = nothing,
    scaling_function = no_scaling,
    land_mask = nothing,
    date0 = nothing,
    mono = true,
)

    # regrid all times and save to hdf5 files
    hd_outfile_root = varname * "_cgll"
    if ClimaComms.iamroot(comms_ctx)
        hdwrite_regridfile_rll_to_cgll(
            comms_ctx,
            datafile_rll,
            varname,
            boundary_space;
            hd_outfile_root = hd_outfile_root,
            mono = mono,
        )
    end
    ClimaComms.barrier(comms_ctx)
    data_dates = load(joinpath(REGRID_DIR, hd_outfile_root * "_times.jld2"), "times")

    # init time tracking info
    current_fields = Fields.zeros(FT, boundary_space), Fields.zeros(FT, boundary_space)
    segment_length = [Int(0)]

    # unless the start file date is specified, find the closest one to the start date
    segment_idx0 =
        segment_idx0 != nothing ? segment_idx0 :
        [argmin(abs.(parse(FT, datetime_to_strdate(date0)) .- parse.(FT, datetime_to_strdate.(data_dates[:]))))]

    return BCFileInfo{FT}(
        comms_ctx,
        hd_outfile_root,
        varname,
        data_dates,
        current_fields,
        scaling_function,
        land_mask,
        deepcopy(segment_idx0),
        segment_idx0,
        segment_length,
        interpolate_daily,
    )

end

no_scaling(x, _info) = swap_space!(x, axes(_info.land_mask))

# IO - monthly
"""
    update_midmonth_data!(date, bcf_info)

Extracts boundary condition data from regridded (to model grid) NetCDF files (which times, depends on the specifications in the `bcf_info` struct).
"""
function update_midmonth_data!(date, bcf_info)
    # monthly count
    FT = float_type(bcf_info)
    all_dates = bcf_info.all_dates
    midmonth_idx = bcf_info.segment_idx[1]
    midmonth_idx0 = bcf_info.segment_idx0[1]
    monthly_fields = bcf_info.monthly_fields
    outfile_root = bcf_info.hd_outfile_root
    scaling_function = bcf_info.scaling_function
    varname = bcf_info.varname
    interpolate_daily = bcf_info.interpolate_daily
    comms_ctx = bcf_info.comms_ctx

    ClimaComms.barrier(comms_ctx)


    if (midmonth_idx == midmonth_idx0) && (Dates.days(date - all_dates[midmonth_idx]) < 0) # for init
        midmonth_idx = bcf_info.segment_idx[1] -= Int(1)
        midmonth_idx = midmonth_idx < Int(1) ? midmonth_idx + Int(1) : midmonth_idx
        @warn "this time period is before BC data - using file from $(all_dates[midmonth_idx0])"
        bcf_info.monthly_fields[1] .= scaling_function(
            hdread_regridfile(comms_ctx, outfile_root, all_dates[Int(midmonth_idx0)], varname),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= deepcopy(bcf_info.monthly_fields[1])
        bcf_info.segment_length .= FT(0)
    elseif Dates.days(date - all_dates[end - 1]) > 0 # for fini
        @warn "this time period is after BC data - using file from $(all_dates[end - 1])"
        bcf_info.monthly_fields[1] .= scaling_function(
            hdread_regridfile(comms_ctx, outfile_root, all_dates[Int(length(all_dates))], varname),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= deepcopy(bcf_info.monthly_fields[1])
        bcf_info.segment_length .= FT(0)
    elseif Dates.days(date - all_dates[Int(midmonth_idx + 1)]) > 2 # throw error when there are closer initial indices for the bc file data that matches this date0
        nearest_idx =
            argmin(abs.(parse(FT, datetime_to_strdate(date)) .- parse.(FT, datetime_to_strdate.(all_dates[:]))))
        @error "init data does not correspond to start date. Try initializing with `SIC_info.segment_idx = midmonth_idx = midmonth_idx0 = $nearest_idx` for this start date" # TODO: do this automatically w a warning
    elseif Dates.days(date - all_dates[Int(midmonth_idx)]) > 0 # date crosses to the next month
        midmonth_idx = bcf_info.segment_idx[1] += Int(1)
        @warn "On $date updating monthly data files: mid-month dates = [ $(all_dates[Int(midmonth_idx)]) , $(all_dates[Int(midmonth_idx+1)]) ]"
        bcf_info.segment_length .= (all_dates[Int(midmonth_idx + 1)] - all_dates[Int(midmonth_idx)]).value
        bcf_info.monthly_fields[1] .= scaling_function(
            hdread_regridfile(comms_ctx, outfile_root, all_dates[Int(midmonth_idx)], varname),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= scaling_function(
            hdread_regridfile(comms_ctx, outfile_root, all_dates[Int(midmonth_idx + 1)], varname),
            bcf_info,
        )
    else
        @error "Check boundary file specification"
    end
end

next_date_in_file(bcfile_info) = bcfile_info.all_dates[bcfile_info.segment_idx[1] + Int(1)]

# IO - daily
"""
    interpolate_midmonth_to_daily(date, bcf_info)

Interpolates linearly between two `Fields` in the `bcf_info` struct, or returns the first Field if interpolation is switched off. 
"""
function interpolate_midmonth_to_daily(date, bcf_info)

    if bcf_info.interpolate_daily && bcf_info.segment_length[1] > FT(0)
        segment_length = bcf_info.segment_length
        segment_idx = bcf_info.segment_idx
        all_dates = bcf_info.all_dates
        monthly_fields = bcf_info.monthly_fields

        return interpol.(
            monthly_fields[1],
            monthly_fields[2],
            FT((date - all_dates[Int(segment_idx[1])]).value),
            FT(segment_length[1]),
        )
    else
        return bcf_info.monthly_fields[1]
    end
end

"""
    interpol(f1::FT, f2::FT, Δt_tt1::FT, Δt_t2t1::FT) where {FT}

Performs linear interpolation of `f` at time `t` within a segment `Δt_t2t1 = (t2 - t1)`, of fields `f1` and `f2`, with `t2 > t1`. 

`Δt_tt1 = (t - t1)`
`f(t1) = f1`

"""
function interpol(f1::FT, f2::FT, Δt_tt1::FT, Δt_t2t1::FT) where {FT}
    interp_fraction = Δt_tt1 / Δt_t2t1
    @assert abs(interp_fraction) <= FT(1) "time interpolation weights must be <= 1, but `interp_fraction` = $interp_fraction"
    return f1 * interp_fraction + f2 * (FT(1) - interp_fraction)
end
