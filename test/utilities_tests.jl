#=
    Unit tests for ClimaCoupler Utilities module
=#

using Test
using ClimaCoupler: Utilities, TestHelper
using ClimaCore: Fields

for FT in (Float32, Float64)
    @testset "test float_type for FT=$FT" begin
        cs = Utilities.CoupledSimulation{FT}(
            nothing, # comms_ctx
            nothing, # dates
            nothing, # boundary_space
            nothing, # fields
            nothing, # parsed_args
            nothing, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # t
            Int(200), # Δt_cpl
            (;), # surface_masks
            (;), # model_sims
            (;), # mode
            (), # diagnostics
        )

        @test Utilities.float_type(cs) == FT
    end

    @testset "test swap_space!" begin
        space1 = TestHelper.create_space(FT, R = FT(6371e3))
        space2 = TestHelper.create_space(FT, R = FT(6371e3) / 2)

        field1 = ones(space1)
        field2 = Utilities.swap_space!(field1, space2)

        @test parent(field1) == parent(field2)
        @test axes(field2) == space2
    end
end
