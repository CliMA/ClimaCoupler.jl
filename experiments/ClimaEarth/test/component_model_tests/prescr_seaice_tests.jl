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
        function test_sea_ice_rhs(; SW_d = 0.0, LW_d = 0.0, T_base = 271.2)
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
                SW_d = CC.Fields.zeros(space) .+ FT(SW_d),
                LW_d = CC.Fields.zeros(space) .+ FT(LW_d),
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

        # TODO: get sigma from parameters
        σ = FT(5.67e-8)

        # check expected dT due to upwelling longwave flux only
        # (zero conduction when T_base == initial T_bulk)
        T_base_eq = IceSlabParameters{FT}().T_freeze - FT(5.0)
        dY, Y, p = test_sea_ice_rhs(T_base = T_base_eq)
        dT_expected =
            (-p.params.ϵ * σ * T_base_eq^4) / (p.params.h * p.params.ρ * p.params.c)
        @test all(T -> T == FT(0) || T ≈ dT_expected, Array(parent(dY.T_bulk)))

        # check expected dT due to downwelling shortwave flux and upwelling longwave flux
        # (again set T_base == initial T_bulk)
        dY, Y, p = test_sea_ice_rhs(SW_d = 1.0, LW_d = 0.0, T_base = T_base_eq)
        dT_expected =
            ((1 - p.params.α) * 1.0 - p.params.ϵ * σ * T_base_eq^4) /
            (p.params.h * p.params.ρ * p.params.c)
        @test all(T -> T == FT(0) || T ≈ dT_expected, Array(parent(dY.T_bulk)))

        # check expected dT due to downwelling and upwelling longwave flux
        # (again set T_base == initial T_bulk)
        dY, Y, p = test_sea_ice_rhs(SW_d = 0.0, LW_d = 2.0, T_base = T_base_eq)
        dT_expected =
            (p.params.ϵ * (2.0 - σ * T_base_eq^4)) / (p.params.h * p.params.ρ * p.params.c)
        @test all(T -> T == FT(0) || T ≈ dT_expected, Array(parent(dY.T_bulk)))

        # check that tendency will not result in above freezing T
        dY, Y, p = test_sea_ice_rhs(SW_d = 0.0, LW_d = 0.0, T_base = 330.0) # Float32 requires a large number here!
        dT_maximum = @. (p.params.T_freeze - Y.T_bulk) / p.dt
        @test minimum(dT_maximum .- dY.T_bulk) >= FT(0.0)

        # check that the correct tendency was added due to basal conductive flux and upwelling longwave flux
        T_base = 269.2
        dY, Y, p = test_sea_ice_rhs(SW_d = 0.0, LW_d = 0.0, T_base = T_base)
        T_bulk = minimum(Y.T_bulk) # get the non-zero temperature value
        (; k_ice, h, ρ, c, T_base, ϵ) = p.params
        dT_expected =
            (k_ice / (h * h * ρ * c)) * (T_base - T_bulk) -
            (ϵ * σ * ice_surface_temperature(T_bulk, T_base)^4) / (h * ρ * c)
        @test minimum(dY) ≈ FT(dT_expected)
        @test maximum(dY) ≈ FT(0)
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
