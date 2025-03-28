#=
    Unit tests for ClimaCoupler TimeManager module
=#
import Test: @testset, @test
import Dates
import ClimaCoupler: Interfacer, TimeManager

@testset "test strdate_to_datetime" begin
    @test TimeManager.strdate_to_datetime("19000101") == Dates.DateTime(1900, 1, 1)
    @test TimeManager.strdate_to_datetime("00000101") == Dates.DateTime(0, 1, 1)
end

@testset "test datetime_to_strdate" begin
    @test TimeManager.datetime_to_strdate(Dates.DateTime(1900, 1, 1)) == "19000101"
    @test TimeManager.datetime_to_strdate(Dates.DateTime(0, 1, 1)) == "00000101"
end
