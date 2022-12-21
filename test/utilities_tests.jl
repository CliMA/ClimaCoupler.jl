#= 
    Unit tests for ClimaCoupler Utilities module
=#

using Test
using ClimaCoupler: Utilities, TestHelper
using ClimaCore: Fields

for FT in (Float32, Float64)
    @testset "test float_type_cs for FT=$FT" begin
        cs = Utilities.CoupledSimulation(
            nothing, # comms_ctx
            nothing, # tspan
            nothing, # dates
            nothing, # boundary_space
            nothing, # fields
            nothing, # parsed_args
            nothing, # conservation_checks
            FT(200), # t
            FT(1 * 24 * 3600), # Δt_cpl
            (;), # surface_masks
            (;), # model_sims
            (;), # mode
            (;), # monthly_3d_diags
            (;), # monthly_2d_diags
        )
        @test Utilities.float_type_cs(cs) == FT
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
