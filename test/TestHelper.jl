#=
    TestHelper

This module defines helper functions, objects, and constants to be used by
various files in the test folder.
=#
module TestHelper

using ClimaCore: Geometry, Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaComms
using NCDatasets

export create_space, gen_ncdata

"""
    create_space(FT; comms_ctx = ClimaComms.SingletonCommsContext(),
        R = FT(6371e3), ne = 4, polynomial_degree = 3, nz = 1)

Initialize a space on a sphere with the given parameters.
Used for debugging and testing.

# Arguments
- FT: [DataType] Float type
- comms_ctx: [ClimaComms.AbstractCommsContext] context used for this operation.
- `R`: [FT] radius of the sphere underlying space.
- `ne`: [Integer] number of elements used in the space's mesh.
- `polynomial_degree`: [Integer] degree of the polynomial used to represent
    the space (number of GLL nodes - 1).
- `nz`: [Integer] number of vertical elements
"""
function create_space(
    FT;
    comms_ctx = ClimaComms.SingletonCommsContext(),
    R = FT(6371e3),
    ne = 4,
    polynomial_degree = 3,
    nz = 1,
)
    domain = Domains.SphereDomain(R)
    mesh = Meshes.EquiangularCubedSphere(domain, ne)

    if comms_ctx isa ClimaComms.SingletonCommsContext
        topology = Topologies.Topology2D(comms_ctx, mesh, Topologies.spacefillingcurve(mesh))
    else
        topology = Topologies.DistributedTopology2D(comms_ctx, mesh, Topologies.spacefillingcurve(mesh))
    end

    Nq = polynomial_degree + 1
    quad = Spaces.Quadratures.GLL{Nq}()
    sphere_space = Spaces.SpectralElementSpace2D(topology, quad)

    if nz > 1
        vertdomain =
            Domains.IntervalDomain(Geometry.ZPoint{FT}(0), Geometry.ZPoint{FT}(100); boundary_tags = (:bottom, :top))
        vertmesh = Meshes.IntervalMesh(vertdomain, nelems = nz)
        vert_center_space = Spaces.CenterFiniteDifferenceSpace(vertmesh)
        return Spaces.ExtrudedFiniteDifferenceSpace(sphere_space, vert_center_space)
    else
        return sphere_space
    end
end

"""
    gen_ncdata(FT, path, varname, val)

Create an NCDataset with lat/lon dimensions containing the value `val` for the
variable `varname`, and store it at `path`.

# Arguments
- `FT`: [DataType] Float type.
- `path`: [String] location to store output datafile.
- `varname`: [Symbol] variable name.
- `val`: [FT] value to store as `varname` at all indices.
"""
function gen_ncdata(FT, path, varname, val)
    isfile(path) ? rm(path) : nothing

    # Create dataset of all ones
    nc = NCDataset(path, "c")

    # Define dataset information
    defDim(nc, "lon", 64)
    defDim(nc, "lat", 32)
    nc.attrib["title"] = "this is an NCDataset containing all 1s on a lat/lon grid"

    # Define variables
    lon = defVar(nc, "lon", FT, ("lon",))
    lat = defVar(nc, "lat", FT, ("lat",))
    v = defVar(nc, varname, FT, ("lon", "lat"))

    # Populate lon and lat
    lon[:] = [i for i in 0.0:(360 / 64):(360 - (360 / 64))]
    lat[:] = [i for i in (-90 + (180 / 64)):(180 / 32):(90 - (180 / 64))]

    # Generate some example data and write it to v
    v[:, :] = [val for i in 1:nc.dim["lon"], j in 1:nc.dim["lat"]]

    close(nc)
end

end
