# this folder contains temporary tests for component code that has not been moved to any src code (location TBD)
using Test
using ClimaCore
using ClimaCore: Fields
using ClimaCoupler.Regridder
import ClimaCoupler.Regridder: binary_mask

include("TestHelper.jl")

# sea ice
include("../experiments/AMIP/modular/components/flux_calculator.jl")
include("../experiments/AMIP/modular/components/ocean/slab_seaice_init.jl")

for FT in (Float32, Float64)
    @testset "test sea-ice energy slab for FT=$FT" begin
        function test_sea_ice_rhs(; F_rad = 0.0, T_base = 271.2, global_mask = 1.0)
            space = TestHelper.create_space(FT)
            params = IceSlabParameters(
                FT(2),  # ice thickness
                FT(900.0), # density of sea ice
                FT(2100.0), # specific heat of sea ice
                FT(T_base), # temperature of sea water at the ice base
                FT(1e-3), # roughness length for momentum
                FT(1e-5), # roughness length for tracers
                FT(271.2), # freezing point of sea water
                FT(2.0),# thermal condictivity of ice
                FT(0.8), # sea ice albedo
            )

            Y = slab_ice_space_init(FT, space, params)
            dY = slab_ice_space_init(FT, space, params) .* FT(0.0)

            ice_fraction = Fields.ones(space) .* FT(global_mask)
            dt = FT(1.0)

            additional_cache = (;
                F_aero = ClimaCore.Fields.zeros(space),
                F_rad = ClimaCore.Fields.zeros(space) .+ FT(F_rad),
                ice_fraction = ice_fraction,
                dt = dt,
            )

            p = (; additional_cache..., params = params)

            ice_rhs!(dY, Y, p, 0)

            return dY, Y, p
        end

        # check that nothing changes with no fluxes
        dY, Y, p = test_sea_ice_rhs()
        @test sum([i for i in extrema(dY)] .≈ [FT(0.0), FT(0.0)]) == 2

        # check that extracting expected T due to input atmopsheric fluxes
        dY, Y, p = test_sea_ice_rhs(F_rad = 1.0)
        dT_expected = -1.0 / (p.params.h * p.params.ρ * p.params.c)
        @test sum([i for i in extrema(dY)] .≈ [FT(dT_expected), FT(dT_expected)]) == 2

        # check that tendency not added if T of ice would have done above freezing
        dY, Y, p = test_sea_ice_rhs(F_rad = 0.0, T_base = 330.0) # Float32 requires a large number here!
        @test sum([i for i in extrema(dY)] .≈ [FT(0.0), FT(0.0)]) == 2

        # check that the correct tendency was added due to basal flux
        dY, Y, p = test_sea_ice_rhs(F_rad = 0.0, T_base = 269.2, global_mask = 1.0)
        dT_expected = -2.0 * p.params.k_ice / (p.params.h * p.params.h * p.params.ρ * p.params.c)
        @test sum([i for i in extrema(dY)] .≈ [FT(dT_expected), FT(dT_expected)]) == 2

        # check that no tendency is applied in a masked case
        dY, Y, p = test_sea_ice_rhs(F_rad = 0.0, T_base = 269.2, global_mask = 0.0)
        @test sum([i for i in extrema(dY)] .≈ [FT(0.0), FT(0.0)]) == 2
    end
end
