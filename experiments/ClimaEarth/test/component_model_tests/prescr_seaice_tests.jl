import Test: @test, @testset
import Dates
import NCDatasets
import ClimaCore as CC
import Thermodynamics.Parameters as TDP
import ClimaParams # required for TDP
import ClimaCoupler

include(joinpath("..", "..", "components", "ocean", "prescr_seaice.jl"))

for FT in (Float32, Float64)
    @testset "test sea-ice energy slab for FT=$FT" begin
        function test_sea_ice_rhs(; F_radiative = 0.0, T_base = 271.2)
            space = CC.CommonSpaces.CubedSphereSpace(
                FT;
                radius = FT(6371e3),
                n_quad_points = 4,
                h_elem = 4,
            )
            params = IceSlabParameters{FT}(T_base = T_base)

            Y = slab_ice_space_init(FT, space, params)
            dY = slab_ice_space_init(FT, space, params) .* FT(0.0)

            dt = FT(1.0)
            t_start = 0.0

            thermo_params = TDP.ThermodynamicsParameters(FT)

            # Set up prescribed sea ice concentration object
            sic_data = try
                joinpath(
                    @clima_artifact("historical_sst_sic"),
                    "MODEL.ICE.HAD187001-198110.OI198111-202206.nc",
                )
            catch error
                @warn "Using lowres SIC. If you want the higher resolution version, you have to obtain it from ClimaArtifacts"
                joinpath(
                    @clima_artifact("historical_sst_sic_lowres"),
                    "MODEL.ICE.HAD187001-198110.OI198111-202206_lowres.nc",
                )
            end
            SIC_timevaryinginput = TimeVaryingInput(
                sic_data,
                "SEAICE",
                space,
                reference_date = Dates.DateTime("20100101", Dates.dateformat"yyyymmdd"),
                file_reader_kwargs = (; preprocess_func = (data) -> data / 100,), ## convert to fraction
            )
            # Get initial SIC values and use them to calculate ice fraction
            SIC_init = CC.Fields.zeros(space)
            evaluate!(SIC_init, SIC_timevaryinginput, t_start)

            cache = (;
                F_turb_energy = CC.Fields.zeros(space),
                F_radiative = CC.Fields.zeros(space) .+ FT(F_radiative),
                area_fraction = SIC_init,
                SIC_timevaryinginput = SIC_timevaryinginput,
                land_fraction = CC.Fields.zeros(space),
                thermo_params = thermo_params,
                dt = dt,
            )

            p = (; cache..., params = params)

            ice_rhs!(dY, Y, p, t_start)

            return dY, Y, p
        end

        # check that nothing changes with no fluxes (zero conduction when T_base == initial T_sfc)
        T_base_eq = IceSlabParameters{FT}().T_freeze - FT(5.0)
        dY, Y, p = test_sea_ice_rhs(T_base = T_base_eq)
        @test all([i for i in extrema(dY)] .≈ [FT(0.0), FT(0.0)])

        # check expected dT due to radiative flux only (again set T_base == initial T_sfc)
        dY, Y, p = test_sea_ice_rhs(F_radiative = 1.0, T_base = T_base_eq)
        dT_expected = -1.0 / (p.params.h * p.params.ρ * p.params.c)
        @test all(extrema(dY) .≈ FT(dT_expected))

        # check that tendency will not result in above freezing T
        dY, Y, p = test_sea_ice_rhs(F_radiative = 0.0, T_base = 330.0) # Float32 requires a large number here!
        dT_maximum = @. (p.params.T_freeze - Y.T_sfc) / p.dt
        @test minimum(dT_maximum .- dY.T_sfc) >= FT(0.0)

        # check that the correct tendency was added due to basal conductive flux
        dY, Y, p = test_sea_ice_rhs(F_radiative = 0.0, T_base = 269.2)
        dT_expected =
            (p.params.k_ice / (p.params.h * p.params.h * p.params.ρ * p.params.c)) *
            (p.params.T_base - minimum(Y.T_sfc))
        @test minimum(dY) ≈ FT(0)
        @test maximum(dY) ≈ FT(dT_expected)
    end

    @testset "dss_state! SeaIceModelSimulation for FT=$FT" begin
        boundary_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6371e3),
            n_quad_points = 4,
            h_elem = 4,
        )

        # construct dss buffer to put in cache
        dss_buffer = CC.Spaces.create_dss_buffer(CC.Fields.zeros(boundary_space))

        # set up objects for test
        u = CC.Fields.FieldVector(;
            state_field1 = CC.Fields.ones(boundary_space),
            state_field2 = CC.Fields.zeros(boundary_space),
        )
        p = (;
            cache_field = CC.Fields.zeros(boundary_space),
            dss_buffer = CC.Spaces.create_dss_buffer(u),
        )
        integrator = (; u, p)
        sim = PrescribedIceSimulation(nothing, integrator)

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
