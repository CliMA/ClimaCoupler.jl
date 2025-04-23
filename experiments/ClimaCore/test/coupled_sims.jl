import Test: @test, @testset, @test_throws
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC

# Load file to test
include("../CoupledSims/coupled_sim.jl")

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
    xdomain = CC.Domains.IntervalDomain(
        CC.Geometry.XPoint(-1.0),
        CC.Geometry.XPoint(1.0),
        periodic = false,
        boundary_names = (:east, :west),
    )
    ydomain = CC.Domains.IntervalDomain(
        CC.Geometry.YPoint(-1.0),
        CC.Geometry.YPoint(1.0),
        periodic = false,
        boundary_names = (:south, :north),
    )
    domain = CC.Domains.RectangleDomain(xdomain, ydomain)
    mesh = CC.Meshes.RectilinearMesh(domain, n1, n2)
    comms_ctx = ClimaComms.SingletonCommsContext()
    grid_topology = CC.Topologies.Topology2D(comms_ctx, mesh)

    quad = CC.Spaces.Quadratures.GLL{Nij}()
    space = CC.Spaces.SpectralElementSpace2D(grid_topology, quad)
    return space
end

@testset "Coupler Interface" begin

    spaceA = spectral_space_2D()
    spaceB = spectral_space_2D(n1 = 2, n2 = 2)

    simA = SimA(ones(spaceA))
    simB = SimB(zeros(spaceB))
    coupler = CouplerState(1.0)

    coupler_add_field!(coupler, :test1, simA.data; write_sim = simA)

    map = CC.Operators.LinearRemap(spaceB, spaceA)
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
