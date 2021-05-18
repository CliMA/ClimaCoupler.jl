using Test
using Random
using CouplerMachine, Dates, Unitful

@testset "Coupler Interface" begin
    Random.seed!(26)

    coupler = CplState()

    data = rand(10,10)
    date = DateTime(2021)
    register_cpl_field!(coupler, :test1, data, nothing, date, u"kg")
    register_cpl_field!(coupler, :test2, data, nothing, date, u"km/hr")
    register_cpl_field!(coupler, :test3, data, nothing, date, u"°C")

    @testset "coupler_get" begin
        @test data === CouplerMachine.coupler_get(coupler, :test1, nothing, date, u"kg")
        # unit conversion
        @test data .* 1000 == CouplerMachine.coupler_get(coupler, :test1, nothing, date, u"g")
        @test data .* (5 / 18) == CouplerMachine.coupler_get(coupler, :test2, nothing, date, u"m/s")
        @test ustrip.(u"°F", data*u"°C") == CouplerMachine.coupler_get(coupler, :test3, nothing, date, u"°F")

        # key not in coupler dict
        @test_throws KeyError CouplerMachine.coupler_get(coupler, :idontexist, nothing, date, u"kg")
        # retreival at inconsistent datetime
        @test_throws ErrorException CouplerMachine.coupler_get(coupler, :test1, nothing, DateTime(2022), u"kg")
        # incompatible units
        @test_throws Unitful.DimensionError CouplerMachine.coupler_get(coupler, :test1, nothing, date, u"°C")
    end

    @testset "coupler_put!" begin
        newdata = rand(10,10)
        newdate = DateTime(2022)
        CouplerMachine.coupler_put!(coupler, :test1, newdata, nothing, newdate, u"kg")

        @test newdata == CouplerMachine.coupler_get(coupler, :test1, nothing, newdate, u"kg")
        # coupler_put! is in-place
        @test newdata !== CouplerMachine.coupler_get(coupler, :test1, nothing, newdate, u"g")
        # unit conversion
        CouplerMachine.coupler_put!(coupler, :test1, newdata, nothing, newdate, u"g")
        @test newdata ≈ CouplerMachine.coupler_get(coupler, :test1, nothing, newdate, u"g")

        # coupler_put! must be to a previously registered field
        @test_throws KeyError CouplerMachine.coupler_put!(coupler, :idontexist, newdata, nothing, newdate, u"kg")
        # incoming data must match dimensions of registered field
        @test_throws DimensionMismatch CouplerMachine.coupler_put!(coupler, :test1, rand(10,5), nothing, newdate, u"kg")
        # incompatible units
        @test_throws Unitful.DimensionError CouplerMachine.coupler_put!(coupler, :test1, newdata, nothing, newdate, u"J/m")
        # coupler_put! updates coupler timestamp
        @test_throws ErrorException CouplerMachine.coupler_get(coupler, :test1, nothing, date, u"kg")
    end

    @testset "grid interpolation" begin
        # placeholder
    end
end
