import ClimaCore
using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields
using NCDatasets
import Pkg; Pkg.add("TempestRemap_jll")
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
function ncreader_rll_to_cgll(FT, datafile_rll, varname; outfile =  "data_cgll.nc", ne = 4, R = 5.0, Nq = 5, clean_exodus = true)
    
    run(`mkdir -p $REGRID_DIR`)
    ds = NCDataset(datafile_rll)
    nlat = ds.dim["lat"]
    nlon = ds.dim["lon"]

    meshfile_rll = joinpath(REGRID_DIR, "mesh_rll.g")
    rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

    meshfile_cgll = joinpath(REGRID_DIR, "mesh_cgll.g")

    domain = Domains.SphereDomain(R)
    mesh = Meshes.EquiangularCubedSphere(domain, ne)
    topology = Topologies.Topology2D(mesh)
    quad = Spaces.Quadratures.GLL{Nq}()
    space = Spaces.SpectralElementSpace2D(topology, quad)
    coords = Fields.coordinate_field(space)

    # write mesh
    write_exodus(meshfile_cgll, topology)

    meshfile_overlap = joinpath(REGRID_DIR, "mesh_overlap.g")
    overlap_mesh(meshfile_overlap, meshfile_rll, meshfile_cgll)

    weightfile = joinpath(REGRID_DIR, "remap_weights.nc")
    remap_weights(
        weightfile,
        meshfile_rll,
        meshfile_cgll,
        meshfile_overlap;
        out_type = "cgll",
        out_np = Nq,
    )

    datafile_cgll = joinpath(REGRID_DIR, outfile)

    apply_remap(datafile_cgll, datafile_rll, weightfile, [varname])

    # read the remapped file
    offline_outarray = NCDataset(datafile_cgll, "r") do ds_wt
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

    clean_exodus ? run(`mkdir -p $REGRID_DIR`) : nothing

    reshape_sparse_to_field!(offline_field, offline_outarray, R)

end

function ncreader_rll_to_cgll_from_space(FT, infile,  varname, h_space)
    R = h_space.topology.mesh.domain.radius
    ne = h_space.topology.mesh.ne
    Nq = Spaces.Quadratures.polynomial_degree(h_space.quadrature_style) + 1

    mask = ncreader_rll_to_cgll(FT, infile,  varname, ne = ne, R = R, Nq = Nq)    
end 

# for AMIP we don't need regridding. WHen we do we re-introduce the ClimaCoreTempestRemap 
function dummmy_remap!(target, source)  # TODO: bring back Tempest regrid
    parent(target) .= parent(source)
end

