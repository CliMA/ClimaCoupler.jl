using Test
using Random
using CouplerMachine, Dates, Unitful

@testset "Coupler Interface" begin
    Random.seed!(26)

    coupler = CouplerState()

    data = rand(10, 10)
    date = DateTime(2021)
    coupler_add_field!(coupler, :test1, data, nothing, date, u"kg")
    coupler_add_field!(coupler, :test2, data, nothing, date, u"km/hr")
    coupler_add_field!(coupler, :test3, data, nothing, date, u"°C")

    @show coupler

    @testset "coupler_get" begin
        @test data === coupler_get(coupler, :test1, nothing, date, u"kg")
        # unit conversion
        @test data .* 1000 == coupler_get(coupler, :test1, nothing, date, u"g")
        @test data .* (5 / 18) == coupler_get(coupler, :test2, nothing, date, u"m/s")
        @test ustrip.(u"°F", data * u"°C") == coupler_get(coupler, :test3, nothing, date, u"°F")

        # key not in coupler dict
        @test_throws KeyError coupler_get(coupler, :idontexist, nothing, date, u"kg")
        # retreival at inconsistent datetime
        @test_throws ErrorException coupler_get(coupler, :test1, nothing, DateTime(2022), u"kg")
        # incompatible units
        @test_throws Unitful.DimensionError coupler_get(coupler, :test1, nothing, date, u"°C")
    end

    @testset "coupler_put!" begin
        newdata = rand(10, 10)
        newdate = DateTime(2022)
        coupler_put!(coupler, :test1, newdata, nothing, newdate, u"kg")

        @test newdata == coupler_get(coupler, :test1, nothing, newdate, u"kg")
        # coupler_put! is in-place; original data array has been modified
        @test data === coupler_get(coupler, :test1, nothing, newdate, u"kg")
        # unit conversion; data is converted from g -> kg when put in coupler
        coupler_put!(coupler, :test1, newdata, nothing, newdate, u"g")
        @test newdata ≈ 1000 * data

        # coupler_put! must be to a previously add_fielded field
        @test_throws KeyError coupler_put!(coupler, :idontexist, newdata, nothing, newdate, u"kg")
        # incoming data must match dimensions of add_fielded field
        @test_throws DimensionMismatch coupler_put!(coupler, :test1, rand(10, 5), nothing, newdate, u"kg")
        # incompatible units
        @test_throws Unitful.DimensionError coupler_put!(coupler, :test1, newdata, nothing, newdate, u"J/m")
        # coupler_put! updates coupler timestamp
        @test_throws ErrorException coupler_get(coupler, :test1, nothing, date, u"kg")
    end

    @testset "grid interpolation" begin
        # placeholder
    end
end
