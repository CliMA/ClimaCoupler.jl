using Test
using Random
using ClimaCoupler, Dates, Unitful
using IntervalSets
using ClimaCore: Domains, Meshes, Geometry, Topologies, Spaces, Fields, Operators

# Load file to test
include("../../experiments/ClimaCore/CoupledSims/coupled_sim.jl")

@testset "Clock" begin
    time_info = (start = 0.0, dt = 0.5, stop = 2.0)
    clock = Clock(time_info...)

    tick!(clock)
    @test clock.time == time_info.dt
    while !stop_time_exceeded(clock)
        tick!(clock)
    end
    @test clock.time == time_info.stop
end


struct SimA{T} <: AbstractSim
    data::T
end
name(::SimA) = :simA

struct SimB{T} <: AbstractSim
    data::T
end
name(::SimB) = :simB

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

    simA = SimA(ones(spaceA))
    simB = SimB(zeros(spaceB))
    coupler = CouplerState(1.0)

    coupler_add_field!(coupler, :test1, simA.data; write_sim = simA)

    map = Operators.LinearRemap(spaceB, spaceA)
    coupler_add_map!(coupler, :simA_to_simB, map)

    @testset "coupler_get" begin
        @test simA.data === coupler_get(coupler, :test1)

        # test remapping
        @test map === get_remap_operator(coupler, name(simB), name(simA))
        @test ones(spaceB) ≈ coupler_get(coupler, :test1, simB)
        target_field = zeros(spaceB)
        coupler_get!(target_field, coupler, :test1, simB)
        @test ones(spaceB) ≈ target_field

        # key not in coupler dict
        @test_throws KeyError coupler_get(coupler, :idontexist)
        @test_throws KeyError coupler_get(coupler, :idontexist, simB)
        @test_throws KeyError coupler_get!(target_field, coupler, :idontexist, simB)
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
        # incoming data must match dimensions/space of added field
        @test_throws ErrorException coupler_put!(coupler, :test1, ones(spaceB), simA)
    end
end
