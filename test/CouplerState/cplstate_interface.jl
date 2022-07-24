using Test
using Random
using ClimaCoupler, Dates, Unitful
using IntervalSets
using ClimaCore: Domains, Meshes, Geometry, Topologies, Spaces, Fields, Operators

struct SimulationA{T} <: ClimaCoupler.AbstractSimulation
    data::T
end
ClimaCoupler.name(::SimulationA) = :simA

struct SimulationB{T} <: ClimaCoupler.AbstractSimulation
    data::T
end
ClimaCoupler.name(::SimulationB) = :simB

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

@testset "Coupler Interface" begin

    spaceA = spectral_space_2D()
    spaceB = spectral_space_2D(n1 = 2, n2 = 2)

    simA = SimulationA(ones(spaceA))
    simB = SimulationB(zeros(spaceB))
    coupler = CouplerState(1.0)

    coupler_add_field!(coupler, :test1, simA.data; write_sim = simA)

    map = Operators.LinearRemap(spaceB, spaceA)
    coupler_add_map!(coupler, map)

    @show coupler

    @testset "coupler_get" begin
        @test simA.data === coupler_get(coupler, :test1)

        # test remapping
        @test map === ClimaCoupler.coupler_get_map(coupler, axes(simB.data), axes(simA.data))
        @test ones(spaceB) ≈ coupler_get(coupler, :test1, axes(simB.data))
        target_field = zeros(spaceB)
        coupler_get!(target_field, coupler, :test1)
        @test ones(spaceB) ≈ target_field

        # key not in coupler dict
        @test_throws KeyError coupler_get(coupler, :idontexist)
        @test_throws KeyError coupler_get(coupler, :idontexist, axes(simB.data))
        @test_throws KeyError coupler_get!(target_field, coupler, :idontexist)
    end

    @testset "coupler_put!" begin
        newdata = zeros(spaceA)
        # simA can write to :test1
        coupler_put!(coupler, :test1, newdata, simA)
        @test newdata == coupler_get(coupler, :test1)
        # coupler_put! is in-place; original data array has been modified
        @test simA.data === coupler_get(coupler, :test1)

        # simB cannot write to :test1
        @test_throws AssertionError coupler_put!(coupler, :test1, newdata, simB)

        # coupler_put! must be to a previously added field
        @test_throws KeyError coupler_put!(coupler, :idontexist, newdata, simA)
    end
end
