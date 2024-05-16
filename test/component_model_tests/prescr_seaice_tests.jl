import Test: @test, @testset
import ClimaCore as CC
import Thermodynamics.Parameters as TDP
import ClimaCoupler
import ClimaCoupler: TestHelper

include(pkgdir(ClimaCoupler, "experiments/ClimaEarth/components/ocean/prescr_seaice.jl"))

for FT in (Float32, Float64)
    @testset "test sea-ice energy slab for FT=$FT" begin
        function test_sea_ice_rhs(; F_radiative = 0.0, T_base = 271.2, global_mask = 1.0)
            space = TestHelper.create_space(FT)
            params = IceSlabParameters{FT}(T_base = T_base)

            Y = slab_ice_space_init(FT, space, params)
            dY = slab_ice_space_init(FT, space, params) .* FT(0.0)

            ice_fraction = CC.Fields.ones(space) .* FT(global_mask)
            dt = FT(1.0)

            thermo_params = TDP.ThermodynamicsParameters(FT)

            additional_cache = (;
                F_turb_energy = CC.Fields.zeros(space),
                F_radiative = CC.Fields.zeros(space) .+ FT(F_radiative),
                area_fraction = ice_fraction,
                q_sfc = CC.Fields.zeros(space),
                ρ_sfc = CC.Fields.ones(space),
                thermo_params = thermo_params,
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
        dY, Y, p = test_sea_ice_rhs(F_radiative = 1.0)
        dT_expected = -1.0 / (p.params.h * p.params.ρ * p.params.c)
        @test sum([i for i in extrema(dY)] .≈ [FT(dT_expected), FT(dT_expected)]) == 2

        # check that tendency not added if T of ice would have done above freezing
        dY, Y, p = test_sea_ice_rhs(F_radiative = 0.0, T_base = 330.0) # Float32 requires a large number here!
        @test sum([i for i in extrema(dY)] .≈ [FT(0.0), FT(0.0)]) == 2

        # check that the correct tendency was added due to basal flux
        dY, Y, p = test_sea_ice_rhs(F_radiative = 0.0, T_base = 269.2, global_mask = 1.0)
        dT_expected = -2.0 * p.params.k_ice / (p.params.h * p.params.h * p.params.ρ * p.params.c)
        @test sum([i for i in extrema(dY)] .≈ [FT(dT_expected), FT(dT_expected)]) == 2

        # check that no tendency is applied in a masked case
        dY, Y, p = test_sea_ice_rhs(F_radiative = 0.0, T_base = 269.2, global_mask = 0.0)
        @test sum([i for i in extrema(dY)] .≈ [FT(0.0), FT(0.0)]) == 2
    end

    @testset "dss_state! SeaIceModelSimulation for FT=$FT" begin
        # use TestHelper to create space
        boundary_space = TestHelper.create_space(FT)

        # construct dss buffer to put in cache
        dss_buffer = CC.Spaces.create_dss_buffer(CC.Fields.zeros(boundary_space))

        # set up objects for test
        integrator = (;
            u = (;
                state_field1 = FT.(CC.Fields.ones(boundary_space)),
                state_field2 = FT.(CC.Fields.zeros(boundary_space)),
            ),
            p = (; cache_field = FT.(CC.Fields.zeros(boundary_space)), dss_buffer = dss_buffer),
        )
        sim = PrescribedIceSimulation(nothing, nothing, nothing, integrator)

        # make field non-constant to check the impact of the dss step
        coords_lat = CC.Fields.coordinate_field(sim.integrator.u.state_field2).lat
        @. sim.integrator.u.state_field2 = sin(coords_lat)

        # apply DSS
        integrator_copy = deepcopy(integrator)
        dss_state!(sim)

        # test that uniform field and cache are unchanged, non-constant is changed
        # note: uniform field is changed slightly by dss
        @test sim.integrator.u.state_field1 ≈ integrator_copy.u.state_field1
        @test sim.integrator.u.state_field2 != integrator_copy.u.state_field2
        @test sim.integrator.p.cache_field == integrator_copy.p.cache_field
    end
end
