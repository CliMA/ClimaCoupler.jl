"""
    Regridder

This module contains functions to regrid information between CC.Spaces.
Many of the functions used in this module call TempestRemap functions
via ClimaCoreTempestRemap wrappers.
"""
module Regridder

import Dates
import JLD2
import NCDatasets
import ClimaComms
import ClimaCore as CC
import ClimaCoreTempestRemap as CCTR
import ..Interfacer, ..Utilities, ..TimeManager
import ClimaUtilities.Regridders

export read_available_dates,
    read_from_hdf5,
    dummmy_remap!,
    remap_field_cgll_to_rll,
    land_fraction,
    update_surface_fractions!,
    combine_surfaces!,
    combine_surfaces_from_sol!,
    binary_mask,
    nans_to_zero,
    truncate_dataset #= Converts NaNs to zeros of the same type. =#



nans_to_zero(v) = isnan(v) ? typeof(v)(0) : v

"""
    read_from_hdf5(REGIRD_DIR, hd_outfile_root, time, varname, comms_ctx)

Read in a variable `varname` from an HDF5 file.
If a CommsContext other than `SingletonCommsContext` is used for `comms_ctx`,
the input HDF5 file must be readable by multiple MPI processes.

# Arguments
- `REGRID_DIR`: [String] directory to save output files in.
- `hd_outfile_root`: [String] root of the output file name.
- `time`: [Dates.DateTime] the timestamp of the data being written.
- `varname`: [String] variable name of data.
- `comms_ctx`: [ClimaComms.AbstractCommsContext] context used for this operation.

# Returns
- Field or FieldVector
"""
function read_from_hdf5(REGRID_DIR, hd_outfile_root, time, varname, comms_ctx)
    hdfreader = CC.InputOutput.HDF5Reader(
        joinpath(REGRID_DIR, hd_outfile_root * "_" * varname * "_" * string(time) * ".hdf5"),
        comms_ctx,
    )

    field = CC.InputOutput.read_field(hdfreader, varname)
    Base.close(hdfreader)
    return field
end

"""
    dummmy_remap!(target, source)

Simple stand-in function for remapping.
For AMIP we don't need regridding of surface model CC.Fields.
When we do, we re-introduce the ClimaCoreTempestRemap remapping functions.

# Arguments
- `target`: [CC.Fields.Field] destination of remapping.
- `source`: [CC.Fields.Field] source of remapping.
"""
function dummmy_remap!(target, source)
    parent(target) .= parent(source)
end

"""
    write_datafile_cc(datafile_cc, field, name)

Write the data stored in `field` to an NCDataset file `datafile_cc`.

# Arguments
- `datafile_cc`: [String] filename of output file.
- `field`: [CC.Fields.Field] to be written to file.
- `name`: [Symbol] variable name.
"""
function write_datafile_cc(datafile_cc, field, name)
    space = axes(field)
    # write data
    NCDatasets.NCDataset(datafile_cc, "c") do nc
        CCTR.def_space_coord(nc, space; type = "cgll")
        nc_field = NCDatasets.defVar(nc, name, Float64, space)
        nc_field[:, 1] = field

        nothing
    end
end

"""
    remap_field_cgll_to_rll(
        name,
        field::CC.Fields.Field,
        remap_tmpdir,
        datafile_rll;
        nlat = 90,
        nlon = 180
    )

Remap an individual FT-valued Field from model (CGLL) nodes to a lat-lon (RLL)
grid using TempestRemap.

# Arguments
- `name`: [Symbol] variable name.
- `field`: [CC.Fields.Field] data to be remapped.
- `remap_tmpdir`: [String] directory used for remapping.
- `datafile_rll`: [String] filename of remapped data output.
"""
function remap_field_cgll_to_rll(name, field::CC.Fields.Field, remap_tmpdir, datafile_rll; nlat = 90, nlon = 180)
    space = axes(field)
    hspace = :topology in propertynames(space) ? space : CC.Spaces.horizontal_space(space)
    Nq = CC.Spaces.Quadratures.polynomial_degree(CC.Spaces.quadrature_style(hspace)) + 1

    # write out our cubed sphere mesh
    meshfile_cc = remap_tmpdir * "/mesh_cubedsphere.g"
    CCTR.write_exodus(meshfile_cc, CC.Spaces.topology(hspace))

    meshfile_rll = remap_tmpdir * "/mesh_rll.g"
    CCTR.rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

    meshfile_overlap = remap_tmpdir * "/mesh_overlap.g"
    CCTR.overlap_mesh(meshfile_overlap, meshfile_cc, meshfile_rll)

    weightfile = remap_tmpdir * "/remap_weights.nc"
    CCTR.remap_weights(weightfile, meshfile_cc, meshfile_rll, meshfile_overlap; in_type = "cgll", in_np = Nq)

    datafile_cc = remap_tmpdir * "/datafile_cc.nc"
    write_datafile_cc(datafile_cc, field, name)

    CCTR.apply_remap( # TODO: this can be done online
        datafile_rll,
        datafile_cc,
        weightfile,
        [string(name)],
    )
end

"""
    read_available_dates(ds::NCDatasets.NCDataset)

Return all the dates in the given NCDataset. The dates are read from the "time"
or "date" datasets. If none is available, return an empty vector.

Code taken from ClimaUtilities
"""
function read_available_dates(ds)
    if "time" in keys(ds.dim)
        return Dates.DateTime.(reinterpret.(Ref(NCDatasets.DateTimeStandard), ds["time"][:]))
    elseif "date" in keys(ds.dim)
        return yyyymmdd_to_datetime.(string.(ds["date"][:]))
    else
        return Dates.DateTime[]
    end
end

"""
    strdate_to_datetime(strdate::String)

Convert from String ("YYYYMMDD") to Date format.

# Arguments
- `yyyymmdd`: [String] to be converted to Date type

Code taken from ClimaUtilities
"""
function yyyymmdd_to_datetime(strdate::String)
    length(strdate) == 8 || error("$strdate does not have the YYYYMMDD format")
    return Dates.DateTime(parse(Int, strdate[1:4]), parse(Int, strdate[5:6]), parse(Int, strdate[7:8]))
end
"""
    function land_fraction(
        FT,
        REGRID_DIR,
        comms_ctx::ClimaComms.AbstractCommsContext,
        infile,
        varname,
        boundary_space;
        outfile_root = "land_sea_cgll",
        mono = false,
        threshold = 0.7,
    )

Initialize a fraction for land/sea classification of grid squares over the space.
With `mono` = `true`, remappings are monotone and conservative, (slower).
With `mono` = `false`, values outside of `threshold` are cutoff (faster).

See https://github.com/CliMA/ClimaCoupler.jl/wiki/ClimaCoupler-Lessons-Learned
    for a detailed comparison of remapping approaches.

# Arguments
- `FT`: [DataType] Float type
- `REGRID_DIR`: [String] directory to save output files in.
- `comms_ctx`: [ClimaComms.AbstractCommsContext] context used for this operation.
- `infile`: [String] filename containing input data.
- `varname`: [Symbol] variable name.
- `boundary_space`: [CC.Spaces.AbstractSpace] over which we are mapping data.
- `outfile_root`: [String] root for output file name.
- `mono`: [Bool] flag for monotone remapping.
- `threshold`: [FT] cutoff value for `binary_mask` when non-monotone remapping.

# Returns
- CC.Fields.Field
"""
function land_fraction(
    FT,
    REGRID_DIR,
    comms_ctx::ClimaComms.AbstractCommsContext,
    infile,
    varname,
    boundary_space;
    mono = false,
    threshold = 0.7,
)

    if ClimaComms.iamroot(comms_ctx)
        Regridders.TempestRegridder(boundary_space, varname, infile; regrid_dir = REGRID_DIR, mono)
    end
    ClimaComms.barrier(comms_ctx)
    # dates are already read in when using Regridders.TempestRegridder, but they are not
    # returned, so we need to read them again
    dates = NCDatasets.NCDataset(infile, "r") do ds
        read_available_dates(ds)
    end
    outfile_root = varname
    fraction = read_from_hdf5(REGRID_DIR, outfile_root, dates[1], varname, comms_ctx)
    fraction = Utilities.swap_space!(boundary_space, fraction) # needed if we are reading from previous run
    return mono ? fraction : binary_mask.(fraction, threshold)
end

"""
    update_surface_fractions!(cs::Interfacer.CoupledSimulation)

Updates dynamically changing area fractions.
Maintains the invariant that the sum of area fractions is 1 at all points.

# Arguments
- `cs`: [Interfacer.CoupledSimulation] containing area fraction information.
"""
function update_surface_fractions!(cs::Interfacer.CoupledSimulation)
    FT = Interfacer.float_type(cs)

    ice_fraction_before = Interfacer.get_field(cs.model_sims.ice_sim, Val(:area_fraction))

    # static fraction
    land_fraction = Interfacer.get_field(cs.model_sims.land_sim, Val(:area_fraction))

    # update component models with dynamic area fractions
    # max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
    Interfacer.update_field!(
        cs.model_sims.ice_sim,
        Val(:area_fraction),
        max.(min.(ice_fraction_before, FT(1) .- land_fraction), FT(0)),
    )
    ice_fraction = Interfacer.get_field(cs.model_sims.ice_sim, Val(:area_fraction))

    Interfacer.update_field!(
        cs.model_sims.ocean_sim,
        Val(:area_fraction),
        max.(FT(1) .- (ice_fraction .+ land_fraction), FT(0)),
    )
    ocean_fraction = Interfacer.get_field(cs.model_sims.ocean_sim, Val(:area_fraction))

    # check that the sum of area fractions is 1
    @assert minimum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ FT(1)
    @assert maximum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ FT(1)
end

"""
    binary_mask(var, threshold)

Converts a number `var` to 1, if `var` is greater or equal than a given `threshold` value,
or 0 otherwise, keeping the same type.

# Arguments
- `var`: [FT] value to be converted.
- `threshold`: [FT] cutoff value for conversions.
"""
binary_mask(var, threshold) = var >= threshold ? one(var) : zero(var)

"""
    binary_mask(var)

Converts a number `var` to 1, if `var` is greater or equal than `eps(FT)`,
or 0 otherwise, keeping the same type.

# Arguments
- `var`: [FT] value to be converted.
"""
binary_mask(var) = binary_mask(var, eps(eltype(var)))

"""
    combine_surfaces!(combined_field::CC.Fields.Field, sims, field_name::Val)

Sums the fields, specified by `field_name`, weighted by the respective area fractions of all
surface simulations. THe result is saved in `combined_field`.

# Arguments
- `combined_field`: [CC.Fields.Field] output object containing weighted values.
- `sims`: [NamedTuple] containing simulations .
- `field_name`: [Val] containing the name Symbol of the field t be extracted by the `Interfacer.get_field` functions.

# Example
- `combine_surfaces!(temp_field, cs.model_sims, Val(:surface_temperature))`
"""
function combine_surfaces!(combined_field::CC.Fields.Field, sims::NamedTuple, field_name::Val)
    combined_field .= eltype(combined_field)(0)
    for sim in sims
        if sim isa Interfacer.SurfaceModelSimulation
            combined_field .+=
                Interfacer.get_field(sim, Val(:area_fraction)) .* nans_to_zero.(Interfacer.get_field(sim, field_name)) # this ensures that unitialized (masked) areas do not affect (TODO: move to mask / remove)
        end
    end
end

"""
    combine_surfaces_from_sol!(combined_field::CC.Fields.Field, fractions::NamedTuple, fields::NamedTuple)

Sums Field objects in `fields` weighted by the respective area fractions, and updates
these values in `combined_field`.
NamedTuples `fields` and `fractions` must have matching field names.
This method can be used to combine fields that were saved in the solution history.

# Arguments
- `combined_field`: [CC.Fields.Field] output object containing weighted values.
- `fractions`: [NamedTuple] containing weights used on values in `fields`.
- `fields`: [NamedTuple] containing values to be weighted by `fractions`.
"""
function combine_surfaces_from_sol!(combined_field::CC.Fields.Field, fractions::NamedTuple, fields::NamedTuple)
    combined_field .= 0
    warn_nans = false
    for surface_name in propertynames(fields) # could use dot here?
        if any(x -> isnan(x), getproperty(fields, surface_name))
            warn_nans = true
        end
        combined_field .+= getproperty(fractions, surface_name) .* nans_to_zero.(getproperty(fields, surface_name))
    end
    warn_nans && @warn "NaNs were detected and converted to zeros."
end

"""
    read_remapped_field(name::Symbol, datafile_latlon::String, lev_name = "z")

Extract data and coordinates from `datafile_latlon`.
"""
function read_remapped_field(name::Symbol, datafile_latlon::String, lev_name = "z")
    out = NCDatasets.NCDataset(datafile_latlon, "r") do nc
        lon = Array(nc["lon"])
        lat = Array(nc["lat"])
        lev = lev_name in keys(nc) ? Array(nc[lev_name]) : Float64(-999)
        var = Array(nc[name])
        coords = (; lon = lon, lat = lat, lev = lev)

        (var, coords)
    end

    return out
end

"""
    truncate_dataset(datafile, filename, varname, datapath_trunc, date0, t_start, t_end, comms_ctx)

Truncates given data set, and constructs a new dataset containing only
the dates that are used in the simulation
"""
function truncate_dataset(
    datafile,
    filename,
    varname,
    datapath_trunc,
    date0,
    t_start,
    t_end,
    comms_ctx::ClimaComms.AbstractCommsContext,
)
    date_start = date0 + Dates.Second(t_start)
    date_end = date0 + Dates.Second(t_start + t_end)

    filename_truncated = replace(
        string(lowercase(filename), "_truncated_data_", string(date_start), string(date_end), ".nc"),
        r":" => "",
    )
    datafile_truncated = joinpath(datapath_trunc, filename_truncated)

    if ClimaComms.iamroot(comms_ctx)
        ds = NCDatasets.NCDataset(datafile, "r")
        dates = ds["time"][:]

        # Find the bounding indices of the dates we need
        (start_id, end_id) = find_idx_bounding_dates(dates, date_start, date_end)

        var_truncated = NCDatasets.nomissing(NCDatasets.view(ds, time = start_id:end_id)[varname])

        # Create new dataset to fill with truncated data
        ds_truncated = NCDatasets.NCDataset(datafile_truncated, "c")

        # Keep all dimensions of original dataset (except for time, which we truncate)
        ds_dim_names = NCDatasets.dimnames(ds[varname])
        for dim_name in ds_dim_names
            dim_name != "time" && NCDatasets.defDim(ds_truncated, dim_name, ds.dim[dim_name])
        end
        dates_truncated = dates[start_id:end_id]
        NCDatasets.defDim(ds_truncated, "time", length(dates_truncated))
        ds_truncated.attrib["title"] = ds.attrib["title"] * " (dates truncated)"

        # Define dimension variables
        for dim_name in ds_dim_names
            if dim_name == "time"
                var = NCDatasets.defVar(ds_truncated, dim_name, dates_truncated, (dim_name,))
            else
                var = NCDatasets.defVar(ds_truncated, dim_name, ds[dim_name][:], (dim_name,))
            end
        end

        # Create variable of interest in new dataset, and fill with input dataset values
        var = NCDatasets.defVar(ds_truncated, varname, var_truncated, ds_dim_names)

        close(ds)
        close(ds_truncated)

        return datafile_truncated
    end
end

"""
    find_idx_bounding_dates(dates, date_start, date_end)

Returns the index range from dates that contains date_start to date_end
"""
function find_idx_bounding_dates(dates, date_start, date_end)
    # if the simulation start date is before our first date in the dataset
    # leave the beginning of the truncated dataset to be first date available
    if date_start < dates[1]
        start_id = 1
        # if the simulation start date is after the last date in the dataset
        # start the truncated dataset at its last possible date
    elseif date_start > last(dates)
        start_id = length(dates)
        # if the simulation start date falls within the range of the dataset
        # find the closest date to the start date and truncate there
    else
        (~, start_id) = findmin(x -> abs(x - date_start), dates)
        # if the closest date is after the start date, add one more date before
        if dates[start_id] > date_start
            start_id = start_id - 1
        end
    end

    # if the simulation end date is before our first date in the dataset
    # truncate the end of the dataset to be the first date
    if date_end < dates[1]
        end_id = 1
        # if the simulation end date is after the last date in the dataset
        # leave the end of the dataset as is
    elseif date_end > last(dates)
        end_id = length(dates)
        # if the simulation end date falls within the range of the dataset
        # find the closest date to the end date and truncate there
    else
        (~, end_id) = findmin(x -> abs(x - date_end), dates)
        # if the closest date is before the end date, add one more date after
        if dates[end_id] < date_end
            end_id = end_id + 1
        end
    end

    return (; start_id, end_id)
end

end # Module
