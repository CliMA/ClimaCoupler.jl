import ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields
using NCDatasets
using TempestRemap_jll
using Test
using ClimaCoreTempestRemap

OUTPUT_DIR = "."

function reshape_sparse_to_field!(field::Fields.Field, in_array::Array, R)
    field_array = parent(field)

    fill!(field_array, zero(eltype(field_array)))
    Nf = size(field_array, 3)

    f = 1
    for (n, row) in enumerate(R.row_indices)
        it, jt, et = (
            view(R.target_idxs[1], n),
            view(R.target_idxs[2], n),
            view(R.target_idxs[3], n),
        )
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
- nc file needs to be in the exodus format
"""
function ncreader_rll_to_cgll(FT, infile, varname, outfile =  "data_cc.nc", ne = 4, R = 5.0, Nq = 5)
    ds = NCDataset(infile)
    nlat = ds.dim["lat"]
    nlon = ds.dim["lon"]

    meshfile_rll = joinpath(OUTPUT_DIR, "mesh_rll.g")
    rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

    meshfile_cgll = joinpath(OUTPUT_DIR, "mesh_cgll.g")

    domain = Domains.SphereDomain(R)
    mesh = Meshes.EquiangularCubedSphere(domain, ne)
    topology = Topologies.Topology2D(mesh)
    quad = Spaces.Quadratures.GLL{Nq}()
    space = Spaces.SpectralElementSpace2D(topology, quad)
    coords = Fields.coordinate_field(space)

    # write mesh
    meshfile_cc = joinpath(OUTPUT_DIR, meshfile_cgll)
    write_exodus(meshfile_cc, topology)

    meshfile_overlap = joinpath(OUTPUT_DIR, "mesh_overlap.g")
    overlap_mesh(meshfile_overlap, meshfile_rll, meshfile_cc)

    weightfile = joinpath(OUTPUT_DIR, "remap_weights.nc")
    remap_weights(
        weightfile,
        meshfile_rll,
        meshfile_cc,
        meshfile_overlap;
        out_type = "cgll",
        out_np = Nq,
    )

    datafile_cc = joinpath(OUTPUT_DIR, outfile)

    apply_remap(datafile_cc, datafile_rll, weightfile, [varname])

    # read the remapped file
    offline_outarray = NCDataset(datafile_cc, "r") do ds_wt
        ds_wt[varname][:][:, 1]
    end

    offline_outarray = FT.(offline_outarray)
    
    field_o = Fields.zeros(FT, space)
    
    # need to populate all nodes
    weights, col_indices, row_indices = NCDataset(weightfile, "r") do ds_wt
        (Array(ds_wt["S"]), Array(ds_wt["col"]), Array(ds_wt["row"]))
    end
    
    out_type = "cgll"
    
    target_unique_idxs =
        out_type == "cgll" ? collect(Spaces.unique_nodes(space)) :
        collect(Spaces.all_nodes(space))
    
    target_unique_idxs_i =
        map(row -> target_unique_idxs[row][1][1], row_indices)
    target_unique_idxs_j =
        map(row -> target_unique_idxs[row][1][2], row_indices)
    target_unique_idxs_e = map(row -> target_unique_idxs[row][2], row_indices)
    
    target_unique_idxs =
        (target_unique_idxs_i, target_unique_idxs_j, target_unique_idxs_e)
        
    R = (;target_idxs = target_unique_idxs, row_indices = row_indices )
    
    offline_field = similar(field_o)
    reshape_sparse_to_field!(offline_field, offline_outarray, R)

end
