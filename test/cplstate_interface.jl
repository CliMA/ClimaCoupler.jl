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

    @testset "get" begin
        @test data === CouplerMachine.get(coupler, :test1, nothing, date, u"kg")
        # unit conversion
        @test data .* 1000 == CouplerMachine.get(coupler, :test1, nothing, date, u"g")
        @test data .* (5 / 18) == CouplerMachine.get(coupler, :test2, nothing, date, u"m/s")
        @test ustrip.(u"°F", data*u"°C") == CouplerMachine.get(coupler, :test3, nothing, date, u"°F")

        # key not in coupler dict
        @test_throws KeyError CouplerMachine.get(coupler, :idontexist, nothing, date, u"kg")
        # retreival at inconsistent datetime
        @test_throws ErrorException CouplerMachine.get(coupler, :test1, nothing, DateTime(2022), u"kg")
        # incompatible units
        @test_throws Unitful.DimensionError CouplerMachine.get(coupler, :test1, nothing, date, u"°C")
    end

    @testset "put!" begin
        newdata = rand(10,10)
        newdate = DateTime(2022)
        CouplerMachine.put!(coupler, :test1, newdata, nothing, newdate, u"kg")

        @test newdata == CouplerMachine.get(coupler, :test1, nothing, newdate, u"kg")
        # put! is in-place
        @test newdata !== CouplerMachine.get(coupler, :test1, nothing, newdate, u"g")
        # unit conversion
        CouplerMachine.put!(coupler, :test1, newdata, nothing, newdate, u"g")
        @test newdata ≈ CouplerMachine.get(coupler, :test1, nothing, newdate, u"g")

        # put! must be to a previously registered field
        @test_throws KeyError CouplerMachine.put!(coupler, :idontexist, newdata, nothing, newdate, u"kg")
        # incoming data must match dimensions of registered field
        @test_throws DimensionMismatch CouplerMachine.put!(coupler, :test1, rand(10,5), nothing, newdate, u"kg")
        # incompatible units
        @test_throws Unitful.DimensionError CouplerMachine.put!(coupler, :test1, newdata, nothing, newdate, u"J/m")
        # put! updates coupler timestamp
        @test_throws ErrorException CouplerMachine.get(coupler, :test1, nothing, date, u"kg")
    end

    @testset "grid interpolation" begin
        # placeholder
    end
end
