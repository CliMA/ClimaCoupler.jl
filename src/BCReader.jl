"""
    BCReader

This module coordinates reading of boundary conditions from NetCDF files,
as well as regridding calls and temporal interpolations from
monthly to daily intervals.
"""
module BCReader

import JLD2
import ClimaComms
import ClimaCore as CC
import ..Utilities, ..Regridder, ..TimeManager

export BCFileInfo,
    float_type_bcf, bcfile_info_init, update_midmonth_data!, next_date_in_file, interpolate_midmonth_to_daily


"""
    BCFileInfo

Stores information specific to each boundary condition from a file and each variable.

# Inputs:
- bcfile_dir::b                   # directory of the BC file
- comms_ctx::X                    # communication context used for MPI
- hd_outfile_root::S              # filename root for regridded data
- varname::V                      # name of the variable
- all_dates::D                    # vector of all dates contained in the original data file
- monthly_fields::C               # tuple of the two monthly fields, that will be used for the daily interpolation
- scaling_function::O             # function that scales, offsets or transforms the raw variable
- land_fraction::M                # fraction with 1 = 100% land, 0 = 100% ocean and/or sea-ice
- segment_idx::Vector{Int}        # index of the monthly data in the file
- segment_idx0::Vector{Int}       # `segment_idx` of the file data that is closest to date0
- segment_length::Vector{Int}     # length of each month segment (used in the daily interpolation)
- interpolate_daily::Bool         # switch to trigger daily interpolation
- mono::Bool                      # flag for monotone remapping of input data
"""
struct BCFileInfo{FT <: Real, B, X, S, V, D, C, O, M, VI}
    bcfile_dir::B
    comms_ctx::X
    hd_outfile_root::S
    varname::V
    all_dates::D
    monthly_fields::C
    scaling_function::O
    land_fraction::M
    segment_idx::VI
    segment_idx0::VI
    segment_length::VI
    interpolate_daily::Bool
    mono::Bool
end

BCFileInfo{FT}(args...) where {FT} = BCFileInfo{FT, typeof.(args[1:9])...}(args...)

float_type_bcf(::BCFileInfo{FT}) where {FT} = FT

"""
    no_scaling(field, bcf_info)

Remap the values of a `field` onto the space of the
`bcf_info`'s land_fraction without scaling.

# Arguments
- `field`: [CC.Fields.Field] contains the values to be remapped.
- `bcf_info`: [BCFileInfo] contains a land_fraction to remap onto the space of.
"""
no_scaling(field::CC.Fields.Field, bcf_info::BCFileInfo{FT}) where {FT} =
    Utilities.swap_space!(zeros(axes(bcf_info.land_fraction)), field)

"""
    bcfile_info_init(
        FT,
        bcfile_dir,
        datafile_rll,
        varname,
        boundary_space,
        comms_ctx;
        interpolate_daily = false,
        segment_idx0 = nothing,
        scaling_function = no_scaling,
        land_fraction = nothing,
        date0 = nothing,
        mono = true,
    )

Regrids from lat-lon grid to cgll grid, saving the output in a new file,
and returns the info packaged in a single struct.

# Arguments
- `FT`: [DataType] Float type.
- `bcfile_dir`: [String] directory the BC file is stored in.
- `datafile_rll`: [String] file containing data to regrid.
- `varname`: [String] name of the variable to be regridded.
- `boundary_space`: [Spaces.AbstractSpace] the space to which we are mapping.
- `comms_ctx`: [ClimaComms.AbstractCommsContext] context used for this operation.
- `interpolate_daily`: [Bool] switch to trigger daily interpolation.
- `segment_idx0`: [Vector{Int}] reference date which, after initialization, refers to the the first file date index used minus 1 (segment_idx[1] - 1)
- `scaling function`: [Function] scales, offsets or transforms `varname`.
- `land_fraction`: [CC.Fields.field] fraction with 1 = land, 0 = ocean / sea-ice.
- `date0`: [Dates.DateTime] start date of the file data.
- `mono`: [Bool] flag for monotone remapping of `datafile_rll`.

# Returns
- `BCFileInfo`
"""
function bcfile_info_init(
    FT,
    bcfile_dir,
    datafile_rll,
    varname,
    boundary_space,
    comms_ctx;
    interpolate_daily = false,
    segment_idx0 = nothing,
    scaling_function = no_scaling,
    land_fraction = nothing,
    date0 = nothing,
    mono = true,
)

    # regrid all times and save to hdf5 files
    hd_outfile_root = varname * "_cgll"
    if ClimaComms.iamroot(comms_ctx)
        Regridder.hdwrite_regridfile_rll_to_cgll(
            FT,
            bcfile_dir,
            datafile_rll,
            varname,
            boundary_space;
            hd_outfile_root = hd_outfile_root,
            mono = mono,
        )
    end
    ClimaComms.barrier(comms_ctx)
    data_dates = JLD2.load(joinpath(bcfile_dir, hd_outfile_root * "_times.jld2"), "times")

    # init time tracking info
    current_fields = CC.Fields.zeros(FT, boundary_space), CC.Fields.zeros(FT, boundary_space)
    segment_length = [Int(0)]

    # unless the start file date is specified, find the closest one to the start date
    segment_idx0 =
        segment_idx0 != nothing ? segment_idx0 :
        [
            argmin(
                abs.(
                    parse(FT, TimeManager.datetime_to_strdate(date0)) .-
                    parse.(FT, TimeManager.datetime_to_strdate.(data_dates[:]))
                ),
            ),
        ]

    return BCFileInfo{FT}(
        bcfile_dir,
        comms_ctx,
        hd_outfile_root,
        varname,
        data_dates,
        current_fields,
        scaling_function,
        land_fraction,
        deepcopy(segment_idx0),
        segment_idx0,
        segment_length,
        interpolate_daily,
        mono,
    )
end

# IO - monthly
"""
    update_midmonth_data!(date, bcf_info::BCFileInfo{FT}) where {FT}

Extracts boundary condition data from regridded (to model grid) NetCDF files.
The times for which data is extracted depends on the specifications in the
`bcf_info` struct).

# Arguments
- `date`: [Dates.DateTime] start date for data.
- `bcf_info`: [BCFileInfo] containing boundary condition data.
"""
function update_midmonth_data!(date, bcf_info::BCFileInfo{FT}) where {FT}
    # monthly count
    (; bcfile_dir, comms_ctx, hd_outfile_root, varname, all_dates, scaling_function) = bcf_info
    segment_idx = bcf_info.segment_idx[1] # index of the current date in the file. [segment_idx, segment_idx+1] indexes the current segment between which we interpolate
    segment_idx0 = bcf_info.segment_idx0[1] # reference index (first segment_idx - 1)

    # upon initialization (segment_idx == segment_idx0) with model date before the final file date
    if (segment_idx == segment_idx0) && !((date - all_dates[end]).value >= 0)
        #  case 1: model date is before the first segment read from file
        if (date - all_dates[segment_idx0]).value < 0
            @warn "This time period is before BC data - using file from $(all_dates[segment_idx0])"
            bcf_info.monthly_fields[1] .= scaling_function(
                Regridder.read_from_hdf5(bcfile_dir, hd_outfile_root, all_dates[Int(segment_idx0)], varname, comms_ctx),
                bcf_info,
            )
            bcf_info.monthly_fields[2] .= deepcopy(bcf_info.monthly_fields[1])
            bcf_info.segment_length .= Int(0)
            bcf_info.segment_idx[1] -= Int(1)
            bcf_info.segment_idx0[1] -= Int(2)

            # case 2: model date is after the first segment read from file
        elseif (date - all_dates[Int(segment_idx0) + 1]).value >= 0
            nearest_idx = argmin(
                abs.(
                    parse(FT, TimeManager.datetime_to_strdate(date)) .-
                    parse.(FT, TimeManager.datetime_to_strdate.(all_dates[:]))
                ),
            )
            @warn "Initializing with `segment_idx = $nearest_idx"
            bcf_info.segment_idx[1] = nearest_idx
            bcf_info.segment_idx0[1] = nearest_idx
            update_midmonth_data!(date, bcf_info)

            # case 3: model date is within the first segment read from file
        elseif (date - all_dates[segment_idx0]).value >= 0
            @warn "On $date updating file data reads: file dates = [ $(all_dates[segment_idx]) , $(all_dates[segment_idx+1]) ]"
            bcf_info.segment_length .= (all_dates[segment_idx + 1] - all_dates[segment_idx]).value
            bcf_info.monthly_fields[1] .= scaling_function(
                Regridder.read_from_hdf5(bcfile_dir, hd_outfile_root, all_dates[segment_idx], varname, comms_ctx),
                bcf_info,
            )
            bcf_info.monthly_fields[2] .= scaling_function(
                Regridder.read_from_hdf5(bcfile_dir, hd_outfile_root, all_dates[segment_idx + 1], varname, comms_ctx),
                bcf_info,
            )
            bcf_info.segment_idx0[1] -= Int(1)
        end

        # case 4: date is at or after the last date in file
    elseif (date - all_dates[end]).value >= 0
        @warn "This time period is after BC data - using file from $(all_dates[end])"
        bcf_info.monthly_fields[1] .= scaling_function(
            Regridder.read_from_hdf5(
                bcfile_dir,
                hd_outfile_root,
                all_dates[Int(length(all_dates))],
                varname,
                comms_ctx,
            ),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= deepcopy(bcf_info.monthly_fields[1])
        bcf_info.segment_length .= Int(0)

        # case 5: model date crosses to the next segment
    elseif (date - all_dates[Int(segment_idx) + 1]).value >= 0
        segment_idx = bcf_info.segment_idx[1] += Int(1)

        bcf_info.segment_length .= (all_dates[Int(segment_idx + 1)] - all_dates[Int(segment_idx)]).value

        bcf_info.monthly_fields[1] .= scaling_function(
            Regridder.read_from_hdf5(bcfile_dir, hd_outfile_root, all_dates[Int(segment_idx)], varname, comms_ctx),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= scaling_function(
            Regridder.read_from_hdf5(bcfile_dir, hd_outfile_root, all_dates[Int(segment_idx + 1)], varname, comms_ctx),
            bcf_info,
        )

        # case 6: undefined condition
    else
        throw(
            ErrorException(
                "Check boundary file specification: segment: $(all_dates[segment_idx]) - $(all_dates[segment_idx+1]), date: $date",
            ),
        )
    end
end

"""
    next_date_in_file(bcf_info)

Returns the next date stored in the file `bcfile_info` struct after the
current date index given by `segment_idx`.
Note: this function does not update `segment_idx`, so repeated calls will
return the same value unless `segment_idx` is modified elsewhere in between.

# Arguments
- `bcf_info`: [BCFileInfo] containing the date information.

# Returns
- Dates.DateTime
"""
next_date_in_file(bcf_info::BCFileInfo{FT}) where {FT} = bcf_info.all_dates[bcf_info.segment_idx[1] + Int(1)]

# IO - daily
"""
    interpolate_midmonth_to_daily(date, bcf_info::BCFileInfo{FT}) where {FT}

Interpolates linearly between two `Fields` in the `bcf_info` struct,
or returns the first Field if interpolation is switched off.

# Arguments
- `date`: [Dates.DateTime] start date for data.
- `bcf_info`: [BCFileInfo] contains fields to be interpolated.

# Returns
- CC.Fields.field
"""
function interpolate_midmonth_to_daily(date, bcf_info::BCFileInfo{FT}) where {FT}
    (; segment_length, segment_idx, all_dates, monthly_fields, interpolate_daily) = bcf_info
    if interpolate_daily && segment_length[1] > FT(0) && date > all_dates[Int(segment_idx[1])]
        Δt_tt1 = FT((date - all_dates[Int(segment_idx[1])]).value)
        interp_fraction = Δt_tt1 / FT(segment_length[1])
        @assert abs(interp_fraction) <= FT(1) "time interpolation weights must be <= 1, but `interp_fraction` = $interp_fraction"
        return interpol.(monthly_fields[1], monthly_fields[2], Δt_tt1, FT(segment_length[1]))
    else
        return monthly_fields[1]
    end
end

"""
    interpol(f1::FT, f2::FT, Δt_tt1::FT, Δt_t2t1::FT)

Performs linear interpolation of `f` at time `t` within
a segment `Δt_t2t1 = (t2 - t1)`, of fields `f1` and `f2`, with `t2 > t1`.

# Arguments
- `f1`: [FT] first value to be interpolated (`f(t1) = f1`).
- `f2`: [FT] second value to be interpolated.
- `Δt_tt1`: [FT] time between `t1` and some `t` (`Δt_tt1 = (t - t1)`).
- `Δt_t2t1`: [FT] time between `t1` and `t2`.

# Returns
- FT
"""
function interpol(f1::FT, f2::FT, Δt_tt1::FT, Δt_t2t1::FT) where {FT}
    interp_fraction = Δt_tt1 / Δt_t2t1
    return f1 * (FT(1) - interp_fraction) + f2 * interp_fraction
end

end
