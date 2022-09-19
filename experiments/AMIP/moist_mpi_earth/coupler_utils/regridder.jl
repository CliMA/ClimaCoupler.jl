import ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields
using NCDatasets

using TempestRemap_jll
using Test
using ClimaCoreTempestRemap

REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath(".", "regrid_tmp/")
rm(REGRID_DIR; recursive = true, force = true)

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
    ncreader_rll_to_cgll(datafile_rll, varname [; outfile = "data_cgll.nc", ne = 4, R = 5.0, Nq = 5])

Reads and regrids data of the `varname` variable from an input NetCDF file and saves it as another NetCDF file. 
The input NetCDF file needs to be `Exodus` formatted, and can contain time-dependent data. 

"""
function ncreader_rll_to_cgll_from_space(datafile_rll, varname, space; outfile = "data_cgll.nc", mono = false)

    outfile_root = mono ? outfile[1:(end - 3)] * "_mono" : outfile[1:(end - 3)]
    datafile_cgll = joinpath(REGRID_DIR, outfile_root * ".g")

    meshfile_rll = joinpath(REGRID_DIR, outfile_root * "_mesh_rll.g")
    meshfile_cgll = joinpath(REGRID_DIR, outfile_root * "_mesh_cgll.g")
    meshfile_overlap = joinpath(REGRID_DIR, outfile_root * "_mesh_overlap.g")
    weightfile = joinpath(REGRID_DIR, outfile_root * "_remap_weights.nc")

    topology = space.topology
    Nq = Spaces.Quadratures.polynomial_degree(space.quadrature_style) + 1

    if isfile(datafile_cgll) == false
        isdir(REGRID_DIR) ? nothing : mkpath(REGRID_DIR)

        ds = NCDataset(datafile_rll)
        nlat = ds.dim["lat"]
        nlon = ds.dim["lon"]

        # write lat-lon mesh
        rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

        # write cgll mesh, overlap mesh and weight file 
        write_exodus(meshfile_cgll, topology)
        overlap_mesh(meshfile_overlap, meshfile_rll, meshfile_cgll)

        # 'in_np = 1' and 'mono = true' arguments ensure mapping is conservative and monotone
        # Note: for a kwarg not followed by a value, set it to true here (i.e. pass 'mono = true' to produce '--mono')
        # Note: out_np = degrees of freedom = polynomial degree + 1

        kwargs = (; out_type = "cgll", out_np = Nq)
        kwargs = mono ? (; (kwargs)..., in_np = mono ? 1 : false, mono = mono) : kwargs
        remap_weights(weightfile, meshfile_rll, meshfile_cgll, meshfile_overlap; kwargs...)
        # remap
        apply_remap(datafile_cgll, datafile_rll, weightfile, [varname])
    else
        @warn "Using the existing $datafile_cgll : check topology is consistent"
    end

    return weightfile, datafile_cgll
end

"""
    ncreader_cgll_sparse_to_field(datafile_cgll, varname, weightfile, t_i_tuple, space; scaling_function = FT_dot, clean_exodus = false)

Given time `t_i_tuple` indices of the NetCDF file data, this reads in the required data of the specified `varname` variable and converts the sparse vector to a `Field` object
The NetCDF file needs to be of the Exodus format and have a time dimension.
"""
function ncreader_cgll_sparse_to_field(
    datafile_cgll,
    varname,
    weightfile,
    t_i_tuple,
    space;
    scaling_function = FT_dot,
    clean_exodus = false,
)
    # read the remapped file
    offline_outvector = NCDataset(datafile_cgll, "r") do ds_wt
        ds_wt[varname][:][:, [t_i_tuple...]] # ncol, times
    end

    # weightfile info needed to populate all nodes and save into fields
    _, _, row_indices = NCDataset(weightfile, "r") do ds_wt
        (Array(ds_wt["S"]), Array(ds_wt["col"]), Array(ds_wt["row"]))
    end

    out_type = "cgll"

    target_unique_idxs = out_type == "cgll" ? collect(Spaces.unique_nodes(space)) : collect(Spaces.all_nodes(space))

    target_unique_idxs_i = map(row -> target_unique_idxs[row][1][1], row_indices)
    target_unique_idxs_j = map(row -> target_unique_idxs[row][1][2], row_indices)
    target_unique_idxs_e = map(row -> target_unique_idxs[row][2], row_indices)

    target_unique_idxs = (target_unique_idxs_i, target_unique_idxs_j, target_unique_idxs_e)

    R = (; target_idxs = target_unique_idxs, row_indices = row_indices)

    # TODO: this could be taken out for fewer allocations? 
    offline_field = Fields.zeros(FT, space)

    offline_fields = ntuple(x -> similar(offline_field), length(t_i_tuple))

    clean_exodus ? run(`mkdir -p $REGRID_DIR`) : nothing

    ntuple(
        x -> scaling_function(reshape_cgll_sparse_to_field!(offline_fields[x], offline_outvector[:, x], R)),
        length(t_i_tuple),
    )
end

FT_dot(x) = FT.(x)

# for AMIP we don't need regridding of surface model fields. When we do, we re-introduce the ClimaCoreTempestRemap 
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

# TODO:
# - add unit tests 
