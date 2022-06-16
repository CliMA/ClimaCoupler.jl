using Test
using Random
using ClimaCoupler, Dates, Unitful
using IntervalSets
using ClimaCore: Domains, Meshes, Topologies, Spaces, Fields

struct SimulationA <: ClimaCoupler.AbstractSimulation end
ClimaCoupler.name(::SimulationA) = :simA

struct SimulationB <: ClimaCoupler.AbstractSimulation end
ClimaCoupler.name(::SimulationB) = :simB

@testset "Coupler Interface" begin
    Random.seed!(26)

    simA = SimulationA()
    simB = SimulationB()
    coupler = CouplerState()

    data = rand(10, 10)
    coupler_add_field!(coupler, :test1, data; write_sim = simA)
    # coupler_add_field!(coupler, :test2, data; write_sim = simB)

    @show coupler

    @testset "coupler_get" begin
        @test data === coupler_get(coupler, :test1)

        # key not in coupler dict
        @test_throws KeyError coupler_get(coupler, :idontexist)
    end

    @testset "coupler_put!" begin
        newdata = rand(10, 10)
        # simA can write to :test1
        coupler_put!(coupler, :test1, newdata, simA)
        @test newdata == coupler_get(coupler, :test1)
        # coupler_put! is in-place; original data array has been modified
        @test data === coupler_get(coupler, :test1)

        # simB cannot write to :test1
        @test_throws AssertionError coupler_put!(coupler, :test1, newdata, simB)

        # coupler_put! must be to a previously added field
        @test_throws KeyError coupler_put!(coupler, :idontexist, newdata, simA)
        # incoming data must match dimensions of added field
        @test_throws DimensionMismatch coupler_put!(coupler, :test1, rand(10, 5), simA)
    end
end