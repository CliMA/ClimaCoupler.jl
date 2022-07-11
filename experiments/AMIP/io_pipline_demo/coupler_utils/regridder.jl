import ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields
using NCDatasets

using TempestRemap_jll
using Test
using ClimaCoreTempestRemap

REGRID_DIR = "regrid_tmp/"

function reshape_sparse_to_field!(field::Fields.Field, in_array::Array, R)
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
    ncreader_rll_to_cgll()
- reads and regrids a variable from input nc file and saves it as another nc file 
- nc file needs to be in the exodus format and it can handle files with multiple times
"""
function ncreader_rll_to_cgll(
    FT,
    datafile_rll,
    varname;
    outfile = "data_cgll.nc",
    ne = 4,
    R = 5.0,
    Nq = 5,
    clean_exodus = true,
)

    outfile_root = outfile[1:(end - 3)]
    datafile_cgll = joinpath(REGRID_DIR, outfile)

    meshfile_rll = joinpath(REGRID_DIR, outfile_root * "_mesh_rll.g")
    meshfile_cgll = joinpath(REGRID_DIR, outfile_root * "_mesh_cgll.g")
    meshfile_overlap = joinpath(REGRID_DIR, outfile_root * "_mesh_overlap.g")
    weightfile = joinpath(REGRID_DIR, outfile_root * "_remap_weights.nc")

    domain = Domains.SphereDomain(R)
    mesh = Meshes.EquiangularCubedSphere(domain, ne)
    topology = Topologies.Topology2D(mesh)
    quad = Spaces.Quadratures.GLL{Nq}()
    space = Spaces.SpectralElementSpace2D(topology, quad)
    coords = Fields.coordinate_field(space)


    if isfile(datafile_cgll) == false
        run(`mkdir -p $REGRID_DIR`)
        ds = NCDataset(datafile_rll)
        nlat = ds.dim["lat"]
        nlon = ds.dim["lon"]

        # write lat-lon mesh
        rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

        # write cgll mesh, overlap mesh and weight file 
        write_exodus(meshfile_cgll, topology)
        overlap_mesh(meshfile_overlap, meshfile_rll, meshfile_cgll)
        remap_weights(weightfile, meshfile_rll, meshfile_cgll, meshfile_overlap; out_type = "cgll", out_np = Nq)

        # remap
        apply_remap(datafile_cgll, datafile_rll, weightfile, [varname])
    else
        @warn "Using the existing $datafile_cgll"
    end

    return weightfile, datafile_cgll
end
"""
    bcfields_from_file
- given time indices, 
- nc file needs to be of the exodus format and havea time dimension
"""
function bcfields_from_file(datafile_cgll, weightfile, t_i_tple, boundary_space, clean_exodus = false)
    # read the remapped file
    offline_outvector = NCDataset(datafile_cgll, "r") do ds_wt
        ds_wt[varname][:][:, [t_i_tple...]] # ncol, times
    end

    offline_outvector = FT.(offline_outvector)

    # weightfile info needed to populate all nodes and save into fields
    weights, col_indices, row_indices = NCDataset(weightfile, "r") do ds_wt
        (Array(ds_wt["S"]), Array(ds_wt["col"]), Array(ds_wt["row"]))
    end

    out_type = "cgll"

    target_unique_idxs = out_type == "cgll" ? collect(Spaces.unique_nodes(boundary_space)) : collect(boundary_Spaces.all_nodes(boundary_space))

    target_unique_idxs_i = map(row -> target_unique_idxs[row][1][1], row_indices)
    target_unique_idxs_j = map(row -> target_unique_idxs[row][1][2], row_indices)
    target_unique_idxs_e = map(row -> target_unique_idxs[row][2], row_indices)

    target_unique_idxs = (target_unique_idxs_i, target_unique_idxs_j, target_unique_idxs_e)

    R = (; target_idxs = target_unique_idxs, row_indices = row_indices)

    # this could be taken out for fewer dynamic allocations? 
    offline_field = Fields.zeros(FT, boundary_space)

    offline_fields = ntuple(x -> similar(offline_field), length(t_i_tple))

    clean_exodus ? run(`mkdir -p $REGRID_DIR`) : nothing

    ntuple( x -> reshape_sparse_to_field!(offline_fields[x], offline_outvector[:,x] , R), length(t_i_tple))
end


function ncreader_rll_to_cgll_from_space(FT, infile, varname, h_space; outfile = "outfile_cgll.nc")
    R = h_space.topology.mesh.domain.radius
    ne = h_space.topology.mesh.ne
    Nq = Spaces.Quadratures.polynomial_degree(h_space.quadrature_style) + 1

    weightfile, datafile_cgll = ncreader_rll_to_cgll(FT, infile, varname, ne = ne, R = R, Nq = Nq, outfile = outfile)
end

# for AMIP we don't need regridding. WHen we do we re-introduce the ClimaCoreTempestRemap 
function dummmy_remap!(target, source)  # TODO: bring back Tempest regrid
    parent(target) .= parent(source)
end



