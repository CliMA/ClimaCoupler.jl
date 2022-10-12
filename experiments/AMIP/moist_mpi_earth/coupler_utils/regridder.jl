
using ClimaCore
using ClimaCore: Domains, Topologies, Meshes, Spaces, Fields, InputOutput
using NCDatasets
using ClimaCoreTempestRemap
using Dates

REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath(".", "regrid_tmp/")
if ClimaComms.iamroot(comms_ctx)
    rm(REGRID_DIR; recursive = true, force = true)
end
ClimaComms.barrier(comms_ctx)

"""
    reshape_cgll_sparse_to_field!(field::Fields.Field, in_array::Array, R)

Reshapes a sparse vector array (raw output of the TempestRemap) to a Field object. 
Redundant nodes are populated using `dss` operations. 
"""
function reshape_cgll_sparse_to_field!(field::Fields.Field, in_array::Array, R)
    field_array = parent(field)

    fill!(field_array, zero(eltype(field_array)))
    Nf = size(field_array, 3)

    f = 1
    for (n, row) in enumerate(R.row_indices)
        it, jt, et = (view(R.target_idxs[1], n), view(R.target_idxs[2], n), view(R.target_idxs[3], n))
        for f in 1:Nf
            field_array[it, jt, f, et] .= in_array[row]
        end
    end
    # broadcast to the redundant nodes using unweighted dss
    topology = Spaces.topology(axes(field))
    Spaces.dss_interior_faces!(topology, Fields.field_values(field))
    Spaces.dss_local_vertices!(topology, Fields.field_values(field))
    return field
end

"""
hdwrite_regridfile_rll_to_cgll(comms_ctx, datafile_rll, varname, space; hd_outfile_root = "data_cgll", mono = false)

Reads and regrids data of the `varname` variable from an input NetCDF file and saves it as another NetCDF file using Tempest Remap. 
The input NetCDF file needs to be `Exodus` formatted, and can contain time-dependent data. The output NetCDF file
is then read back, the output arrays converted into Fields and saved as HDF5 files (one per time slice). This function should 
be called by the root process. The regridded HDF5 output is readable by multiple MPI processes. 

"""
function hdwrite_regridfile_rll_to_cgll(
    comms_ctx,
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
    space_undistributed = ClimaCore.Spaces.SpectralElementSpace2D(topology, ClimaCore.Spaces.Quadratures.GLL{Nq}())

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
        kwargs = mono ? (; (kwargs)..., in_np = mono ? 1 : false, mono = mono) : kwargs
        remap_weights(weightfile, meshfile_rll, meshfile_cgll, meshfile_overlap; kwargs...)
        # remap
        apply_remap(datafile_cgll, datafile_rll, weightfile, [varname])
    else
        @warn "Using the existing $datafile_cgll : check topology is consistent"
    end

    function get_time(ds)
        if "time" in ds
            data_dates = Dates.DateTime.(ds["time"][:])
        elseif "date" in ds
            data_dates = strdate_to_datetime.(string.(ds["date"][:]))
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

    # weightfile info needed to populate all nodes and save into fields with sparse matrices
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

    # save save_hdf5 # TODO: extend write! to handle time-dependent fields
    map(x -> save_remap_hdf5(hd_outfile_root, times[x], offline_fields[x], varname), 1:length(times))
    jldsave(joinpath(REGRID_DIR, hd_outfile_root * "_times.jld2"); times = times)

    return times
end

"""
    save_remap_hdf5(hd_outfile_root, tx, field, varname)

Helper to save individual hdf5 files after remapping.
"""
function save_remap_hdf5(hd_outfile_root, tx, field, varname)
    t = Dates.datetime2unix.(tx)
    hdfwriter = InputOutput.HDF5Writer(joinpath(REGRID_DIR, hd_outfile_root * "_" * string(tx) * ".hdf5"))

    InputOutput.HDF5.write_attribute(hdfwriter.file, "unix time", t) # TODO: a better way to write metadata, CMIP convention
    InputOutput.write!(hdfwriter, field, string(varname))
    Base.close(hdfwriter)
end

function hdread_regridfile(comms_ctx, hd_outfile_root, time, varname)
    hdfreader = InputOutput.HDF5Reader(joinpath(REGRID_DIR, hd_outfile_root * "_" * string(time) * ".hdf5"), comms_ctx)
    field = InputOutput.read_field(hdfreader, varname)
    Base.close(hdfreader)
    return field
end

function dummmy_remap!(target, source)
    parent(target) .= parent(source)
end

"""
    remap_field_cgll2rll(name::Symbol, field::Fields.Field, remap_tmpdir, datafile_latlon)
Remap individual Field from model (CGLL) nodes to a lat-lon (RLL) grid using TempestRemap
"""
function remap_field_cgll2rll(name::Symbol, field::Fields.Field, remap_tmpdir, datafile_latlon; nlat = 90, nlon = 180)

    space = axes(field)
    hspace = :topology in propertynames(space) ? space : space.horizontal_space
    topology = hspace.topology
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
        datafile_latlon,
        datafile_cc,
        weightfile,
        [string(name)],
    )
end
function write_datafile_cc(datafile_cc, field, name)
    space = axes(field)
    # write data
    NCDataset(datafile_cc, "c") do nc
        def_space_coord(nc, space; type = "cgll")
        # nc_time = def_time_coord(nc)
        nc_field = defVar(nc, name, Float64, space)
        nc_field[:, 1] = field

        nothing
    end
end


function create_space(; R = FT(6371e3), ne = 4, polynomial_degree = 3)
    domain = Domains.SphereDomain(R)
    mesh = Meshes.EquiangularCubedSphere(domain, ne)
    topology = Topologies.Topology2D(mesh)
    Nq = polynomial_degree + 1
    quad = Spaces.Quadratures.GLL{Nq}()
    space = Spaces.SpectralElementSpace2D(topology, quad)
end # debugging tool, also to be used in tests
