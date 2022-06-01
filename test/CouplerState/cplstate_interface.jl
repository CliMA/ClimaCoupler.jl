using Test
using Random
using ClimaCoupler, Dates, Unitful

@testset "Coupler Interface" begin
    Random.seed!(26)

    coupler = CouplerState()

    data = rand(10, 10)
    coupler_add_field!(coupler, :test1, data)
    coupler_add_field!(coupler, :test2, data)
    coupler_add_field!(coupler, :test3, data)

    @show coupler

    @testset "coupler_get" begin
        @test data === coupler_get(coupler, :test1)

        # key not in coupler dict
        @test_throws KeyError coupler_get(coupler, :idontexist)
    end

    @testset "coupler_put!" begin
        newdata = rand(10, 10)
        coupler_put!(coupler, :test1, newdata)

        @test newdata == coupler_get(coupler, :test1)
        # coupler_put! is in-place; original data array has been modified
        @test data === coupler_get(coupler, :test1)

        # coupler_put! must be to a previously add_fielded field
        @test_throws KeyError coupler_put!(coupler, :idontexist, newdata)
        # incoming data must match dimensions of add_fielded field
        @test_throws DimensionMismatch coupler_put!(coupler, :test1, rand(10, 5))
    end
end
