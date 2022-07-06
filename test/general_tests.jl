using Test
using ClimaCoupler
import ClimaCore: Fields, Domains, Topologies, Meshes, Spaces, Geometry
using IntervalSets

function spectral_space_2D(; n1 = 1, n2 = 1, Nij = 4)
    domain = Domains.RectangleDomain(
        Geometry.XPoint(-1.0) .. Geometry.XPoint(1.0),
        Geometry.YPoint(-1.0) .. Geometry.YPoint(1.0),
        x1periodic = false,
        x2periodic = false,
        x1boundary = (:east, :west),
        x2boundary = (:south, :north),
    )
    mesh = Meshes.RectilinearMesh(domain, n1, n2)
    grid_topology = Topologies.Topology2D(mesh)

    quad = Spaces.Quadratures.GLL{Nij}()
    space = Spaces.SpectralElementSpace2D(grid_topology, quad)
    return space
end

@testset "Enable Mismatched Space Broadcasting" begin
    space1 = spectral_space_2D()
    space2 = spectral_space_2D()
    field1 = ones(space1)
    field2 = 2 .* ones(space2)
    @test Fields.is_diagonalized_spaces(typeof(space1), typeof(space2))
    @test_nowarn field1 .= field2
    @test parent(field1) == parent(field2)
end
