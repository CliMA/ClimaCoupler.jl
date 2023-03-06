"""
    Regridder

This module contains functions to regrid information between spaces.
Many of the functions used in this module call TempestRemap functions
via ClimaCoreTempestRemap wrappers.
"""
module Regridder

using ..Utilities
using ..TimeManager
using ClimaCore: Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaComms
using NCDatasets
using ClimaCoreTempestRemap
using Dates
using JLD2

export write_to_hdf5,
    read_from_hdf5,
    dummmy_remap!,
    remap_field_cgll_to_rll,
    land_sea_mask,
    update_masks!,
    combine_surfaces!,
    binary_mask,
    nans_to_zero


#= Converts NaNs to zeros of the same type. =#
nans_to_zero(v) = isnan(v) ? typeof(v)(0) : v

"""
    reshape_cgll_sparse_to_field!(field::Fields.Field, in_array::Array, R)

Reshapes a sparse vector array `in_array` (CGLL, raw output of the TempestRemap),
and uses its data to populate the input Field object `field`.
Redundant nodes are populated using `dss` operations.

# Arguments
- `field`: [Fields.Field] object populated with the input array.
- `in_array`: [Array] input used to fill `field`.
- `R`: [NamedTuple] containing `target_idxs` and `row_indices` used for indexing.
"""
function reshape_cgll_sparse_to_field!(field::Fields.Field, in_array::Array, R)
    field_array = parent(field)

    fill!(field_array, zero(eltype(field_array)))
    Nf = size(field_array, 3)

    for (n, row) in enumerate(R.row_indices)
        it, jt, et = (view(R.target_idxs[1], n), view(R.target_idxs[2], n), view(R.target_idxs[3], n))
        for f in 1:Nf
            field_array[it, jt, f, et] .= in_array[row]
        end
    end

    # broadcast to the redundant nodes using unweighted dss
    space = axes(field)
    topology = Spaces.topology(space)
    hspace = Spaces.horizontal_space(space)
    Spaces.dss2!(Fields.field_values(field), topology, hspace.quadrature_style)
end

"""
    hdwrite_regridfile_rll_to_cgll(
        FT,
        REGRID_DIR,
        datafile_rll,
        varname,
        space;
        hd_outfile_root = "data_cgll",
        mono = false,
    )

Reads and regrids data of the `varname` variable from an input NetCDF file and
saves it as another NetCDF file using Tempest Remap.
The input NetCDF fileneeds to be `Exodus` formatted, and can contain
time-dependent data. The output NetCDF file is then read back, the output
arrays converted into Fields and saved as HDF5 files (one per time slice).
This function should be called by the root process.
The saved regridded HDF5 output is readable by multiple MPI processes.

# Arguments
- `FT`: [DataType] Float type.
- `REGRID_DIR`: [String] directory to save output files in.
- `datafile_rll`: [String] filename of RLL dataset to be mapped to CGLL.
- `varname`: [String] the name of the variable to be remapped.
- `space`: [Spaces.AbstractSpace] the space to which we are mapping.
- `hd_outfile_root`: [String] root of the output file name.
- `mono`: [Bool] flag to specify monotone remapping.
"""
function hdwrite_regridfile_rll_to_cgll(
    FT,
    REGRID_DIR,
    datafile_rll,
    varname,
    space;
    hd_outfile_root = "data_cgll",
    mono = false,
)
    out_type = "cgll"

    outfile = hd_outfile_root * ".nc"
    outfile_root = mono ? outfile[1:(end - 3)] * "_mono" : outfile[1:(end - 3)]
    datafile_cgll = joinpath(REGRID_DIR, outfile_root * ".g")

    meshfile_rll = joinpath(REGRID_DIR, outfile_root * "_mesh_rll.g")
    meshfile_cgll = joinpath(REGRID_DIR, outfile_root * "_mesh_cgll.g")
    meshfile_overlap = joinpath(REGRID_DIR, outfile_root * "_mesh_overlap.g")
    weightfile = joinpath(REGRID_DIR, outfile_root * "_remap_weights.nc")

    topology = Topologies.Topology2D(space.topology.mesh, Topologies.spacefillingcurve(space.topology.mesh))
    Nq = Spaces.Quadratures.polynomial_degree(space.quadrature_style) + 1
    space_undistributed = Spaces.SpectralElementSpace2D(topology, Spaces.Quadratures.GLL{Nq}())

    if isfile(datafile_cgll) == false
        isdir(REGRID_DIR) ? nothing : mkpath(REGRID_DIR)

        nlat, nlon = NCDataset(datafile_rll) do ds
            (ds.dim["lat"], ds.dim["lon"])
        end
        # write lat-lon mesh
        rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

        # write cgll mesh, overlap mesh and weight file
        write_exodus(meshfile_cgll, topology)
        overlap_mesh(meshfile_overlap, meshfile_rll, meshfile_cgll)

        # 'in_np = 1' and 'mono = true' arguments ensure mapping is conservative and monotone
        # Note: for a kwarg not followed by a value, set it to true here (i.e. pass 'mono = true' to produce '--mono')
        # Note: out_np = degrees of freedom = polynomial degree + 1
        kwargs = (; out_type = out_type, out_np = Nq)
        kwargs = mono ? (; (kwargs)..., in_np = 1, mono = mono) : kwargs
        remap_weights(weightfile, meshfile_rll, meshfile_cgll, meshfile_overlap; kwargs...)
        apply_remap(datafile_cgll, datafile_rll, weightfile, [varname])
    else
        @warn "Using the existing $datafile_cgll : check topology is consistent"
    end

    function get_time(ds)
        if "time" in ds
            data_dates = Dates.DateTime.(ds["time"][:])
        elseif "date" in ds
            data_dates = TimeManager.strdate_to_datetime.(string.(ds["date"][:]))
        else
            @warn "No dates available in file $datafile_rll"
            data_dates = [Dates.DateTime(0)]
        end
    end

    # read the remapped file with sparse matrices
    offline_outvector, times = NCDataset(datafile_cgll, "r") do ds_wt
        (
            offline_outvector = ds_wt[varname][:][:, :], # ncol, times
            times = get_time(ds_wt),
        )
    end

    # weightfile info needed to populate all nodes and save into fields with
    #  sparse matrices
    _, _, row_indices = NCDataset(weightfile, "r") do ds_wt
        (Array(ds_wt["S"]), Array(ds_wt["col"]), Array(ds_wt["row"]))
    end

    target_unique_idxs =
        out_type == "cgll" ? collect(Spaces.unique_nodes(space_undistributed)) :
        collect(Spaces.all_nodes(space_undistributed))
    target_unique_idxs_i = map(row -> target_unique_idxs[row][1][1], row_indices)
    target_unique_idxs_j = map(row -> target_unique_idxs[row][1][2], row_indices)
    target_unique_idxs_e = map(row -> target_unique_idxs[row][2], row_indices)
    target_unique_idxs = (target_unique_idxs_i, target_unique_idxs_j, target_unique_idxs_e)

    R = (; target_idxs = target_unique_idxs, row_indices = row_indices)

    offline_field = Fields.zeros(FT, space_undistributed)

    offline_fields = ntuple(x -> similar(offline_field), length(times))

    ntuple(x -> reshape_cgll_sparse_to_field!(offline_fields[x], offline_outvector[:, x], R), length(times))

    # TODO: extend write! to handle time-dependent fields
    map(
        x -> write_to_hdf5(
            REGRID_DIR,
            hd_outfile_root,
            times[x],
            offline_fields[x],
            varname,
            ClimaComms.SingletonCommsContext(),
        ),
        1:length(times),
    )
    jldsave(joinpath(REGRID_DIR, hd_outfile_root * "_times.jld2"); times = times)
end

"""
    write_to_hdf5(REGRID_DIR, hd_outfile_root, time, field, varname, comms_ctx)

Function to save individual HDF5 files after remapping.
If a CommsContext other than `SingletonCommsContext` is used for `comms_ctx`,
the HDF5 output is readable by multiple MPI processes.

# Arguments
- `REGRID_DIR`: [String] directory to save output files in.
- `hd_outfile_root`: [String] root of the output file name.
- `time`: [Dates.DateTime] the timestamp of the data being written.
- `field`: [Fields.Field] object to be written.
- `varname`: [String] variable name of data.
- `comms_ctx`: [ClimaComms.AbstractCommsContext] context used for this operation.
"""
function write_to_hdf5(REGRID_DIR, hd_outfile_root, time, field, varname, comms_ctx)
    t = Dates.datetime2unix.(time)
    hdfwriter = InputOutput.HDF5Writer(joinpath(REGRID_DIR, hd_outfile_root * "_" * string(time) * ".hdf5"), comms_ctx)

    InputOutput.HDF5.write_attribute(hdfwriter.file, "unix time", t) # TODO: a better way to write metadata, CMIP convention
    InputOutput.write!(hdfwriter, field, string(varname))
    Base.close(hdfwriter)
end

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
    hdfreader = InputOutput.HDF5Reader(joinpath(REGRID_DIR, hd_outfile_root * "_" * string(time) * ".hdf5"), comms_ctx)

    field = InputOutput.read_field(hdfreader, varname)
    Base.close(hdfreader)
    return field
end

"""
    dummmy_remap!(target, source)

Simple stand-in function for remapping.
For AMIP we don't need regridding of surface model fields.
When we do, we re-introduce the ClimaCoreTempestRemap remapping functions.

# Arguments
- `target`: [Fields.Field] destination of remapping.
- `source`: [Fields.Field] source of remapping.
"""
function dummmy_remap!(target, source)
    parent(target) .= parent(source)
end

"""
    write_datafile_cc(datafile_cc, field, name)

Write the data stored in `field` to an NCDataset file `datafile_cc`.

# Arguments
- `datafile_cc`: [String] filename of output file.
- `field`: [Fields.Field] to be written to file.
- `name`: [Symbol] variable name.
"""
function write_datafile_cc(datafile_cc, field, name)
    space = axes(field)
    # write data
    NCDataset(datafile_cc, "c") do nc
        def_space_coord(nc, space; type = "cgll")
        nc_field = defVar(nc, name, Float64, space)
        nc_field[:, 1] = field

        nothing
    end
end

"""
    remap_field_cgll_to_rll(
        name,
        field::Fields.Field,
        remap_tmpdir,
        datafile_rll;
        nlat = 90,
        nlon = 180
    )

Remap an individual FT-valued Field from model (CGLL) nodes to a lat-lon (RLL)
grid using TempestRemap.

# Arguments
- `name`: [Symbol] variable name.
- `field`: [Fields.Field] data to be remapped.
- `remap_tmpdir`: [String] directory used for remapping.
- `datafile_rll`: [String] filename of remapped data output.
"""
function remap_field_cgll_to_rll(name, field::Fields.Field, remap_tmpdir, datafile_rll; nlat = 90, nlon = 180)
    space = axes(field)
    hspace = :topology in propertynames(space) ? space : space.horizontal_space
    Nq = Spaces.Quadratures.polynomial_degree(hspace.quadrature_style) + 1

    # write out our cubed sphere mesh
    meshfile_cc = remap_tmpdir * "/mesh_cubedsphere.g"
    write_exodus(meshfile_cc, hspace.topology)

    meshfile_rll = remap_tmpdir * "/mesh_rll.g"
    rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

    meshfile_overlap = remap_tmpdir * "/mesh_overlap.g"
    overlap_mesh(meshfile_overlap, meshfile_cc, meshfile_rll)

    weightfile = remap_tmpdir * "/remap_weights.nc"
    remap_weights(weightfile, meshfile_cc, meshfile_rll, meshfile_overlap; in_type = "cgll", in_np = Nq)

    datafile_cc = remap_tmpdir * "/datafile_cc.nc"
    write_datafile_cc(datafile_cc, field, name)

    apply_remap( # TODO: this can be done online
        datafile_rll,
        datafile_cc,
        weightfile,
        [string(name)],
    )
end

"""
    function land_sea_mask(
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

Initialize a mask for land/sea classification of grid squares over the space.
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
- `boundary_space`: [Spaces.AbstractSpace] over which we are mapping data.
- `outfile_root`: [String] root for output file name.
- `mono`: [Bool] flag for monotone remapping.
- `threshold`: [FT] cutoff value for `binary_mask` when non-monotone remapping.

# Returns
- Fields.Field
"""
function land_sea_mask(
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

    if ClimaComms.iamroot(comms_ctx)
        hdwrite_regridfile_rll_to_cgll(
            FT,
            REGRID_DIR,
            infile,
            varname,
            boundary_space;
            hd_outfile_root = outfile_root,
            mono = mono,
        )
    end
    ClimaComms.barrier(comms_ctx)
    file_dates = JLD2.load(joinpath(REGRID_DIR, outfile_root * "_times.jld2"), "times")
    mask = read_from_hdf5(REGRID_DIR, outfile_root, file_dates[1], varname, comms_ctx)
    mask = swap_space!(zeros(boundary_space), mask) # needed if we are reading from previous run
    return mono ? mask : binary_mask.(mask, threshold = threshold)
end

"""
    update_masks!(cs::CoupledSimulation)

Updates dynamically changing masks.
Maintains the invariant that the sum of masks is 1 at all points.

# Arguments
- `cs`: [CoupledSimulation] containing mask information.
"""
function update_masks!(cs::CoupledSimulation)
    # dynamic masks
    ice_d = cs.model_sims.ice_sim.integrator.p.ice_mask
    FT = eltype(ice_d)

    # static mask
    land_s = cs.surface_masks.land

    # max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
    cs.surface_masks.ice .= max.(min.(ice_d, FT(1) .- land_s), FT(0))
    cs.surface_masks.ocean .= max.(FT(1) .- (cs.surface_masks.ice .+ land_s), FT(0))

    @assert minimum(cs.surface_masks.ice .+ cs.surface_masks.land .+ cs.surface_masks.ocean) ≈ FT(1)
    @assert maximum(cs.surface_masks.ice .+ cs.surface_masks.land .+ cs.surface_masks.ocean) ≈ FT(1)

end


"""
    binary_mask(var::FT; threshold = 0.5)

Converts a number `var` to 1, if `var` is greater than a given `threshold` value, or 0 otherwise, keeping the same type.

# Arguments
- `var`: [FT] value to be converted.
- `threshold`: [Float] cutoff value for conversions.
"""
binary_mask(var::FT; threshold = 0.5) where {FT} = var > FT(threshold) ? FT(1) : FT(0)

"""
    combine_surfaces!(combined_field::Fields.Field, masks::NamedTuple, fields::NamedTuple)

Sums Field objects in `fields` weighted by the respective masks, and updates
these values in `combined_field`.
NamedTuples `fields` and `masks` must have matching field names.

# Arguments
- `combined_field`: [Fields.Field] output object containing weighted values.
- `masks`: [NamedTuple] containing weights used on values in `fields`.
- `fields`: [NamedTuple] containing values to be weighted by `masks`.
"""
function combine_surfaces!(combined_field::Fields.Field, masks::NamedTuple, fields::NamedTuple)
    combined_field .= eltype(combined_field)(0)
    warn_nans = false
    for surface_name in propertynames(fields) # could use dot here?
        if any(x -> isnan(x), getproperty(fields, surface_name))
            warn_nans = true
        end
        combined_field .+= getproperty(masks, surface_name) .* nans_to_zero.(getproperty(fields, surface_name))
    end
    warn_nans && @warn "NaNs were detected and converted to zeros."
end

end # Module
