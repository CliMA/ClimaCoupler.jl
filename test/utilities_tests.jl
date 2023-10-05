#=
    Unit tests for ClimaCoupler Utilities module
=#

using Test
using ClimaCoupler: Utilities, TestHelper
using ClimaCore: Fields

for FT in (Float32, Float64)
    @testset "test swap_space!" begin
        space1 = TestHelper.create_space(FT, R = FT(6371e3))
        space2 = TestHelper.create_space(FT, R = FT(6371e3))

        field1 = ones(space1)
        field2 = ones(space2)

        field2 = Utilities.swap_space!(field2, field1)

        @test parent(field1) == parent(field2)
        @test axes(field2) == space2
    end
end
