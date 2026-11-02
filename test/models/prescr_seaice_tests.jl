import Test: @test, @testset
import Dates
import NCDatasets
import ClimaCore as CC
import Thermodynamics.Parameters as TDP
import ClimaParams as CP # required for TDP
import ClimaCoupler: Interfacer, Models, Utilities
import ClimaTimeSteppers as CTS
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.ClimaArtifacts: @clima_artifact

for FT in (Float32, Float64)
    @testset "test sea-ice energy slab for FT=$FT" begin
        coupled_param_dict = CP.create_toml_dict(FT)
        thermo_params = TDP.ThermodynamicsParameters(coupled_param_dict)

        space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = coupled_param_dict["planet_radius"], # in meters
            n_quad_points = 4,
            h_elem = 4,
        )

        function test_sea_ice_rhs(
            coupled_param_dict,
            thermo_params,
            space;
            SW_d = 0.0,
            LW_d = 0.0,
            T_base = 271.2,
        )
            params = Models.IceSlabParameters{FT}(coupled_param_dict; T_base)

            Y = Models.slab_ice_space_init(FT, space, params, nothing)
            dY = Models.slab_ice_space_init(FT, space, params, nothing) .* FT(0.0)

            dt = FT(1.0)
            t_start = 0.0

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
                start_date = Dates.DateTime("20100101", Dates.dateformat"yyyymmdd"),
                file_reader_kwargs = (; preprocess_func = (data) -> data / 100,), ## convert to fraction
            )
            # Get initial SIC values and use them to calculate ice fraction
            SIC_init = CC.Fields.zeros(space)
            evaluate!(SIC_init, SIC_timevaryinginput, t_start)

            α_ice = CC.Fields.zeros(space)
            @. α_ice = Models.ice_surface_albedo(
                Models.ice_surface_temperature(
                    Y.T_bulk,
                    params.T_base,
                    params.T_sfc_min,
                    params.T_freeze,
                ),
                params.T_freeze,
                params.ΔT_melt,
                params.α_cold,
                params.α_melt,
            )

            cache = (;
                F_turb_energy = CC.Fields.zeros(space),
                SW_d = CC.Fields.zeros(space) .+ FT(SW_d),
                LW_d = CC.Fields.zeros(space) .+ FT(LW_d),
                α_ice = α_ice,
                area_fraction = SIC_init,
                SIC_timevaryinginput = SIC_timevaryinginput,
                land_fraction = CC.Fields.zeros(space),
                thermo_params = thermo_params,
                dt = dt,
                binary_area_fraction = true,
                domain_type = "global",
            )

            p = (; cache..., params = params)

            Models.ice_rhs!(dY, Y, p, t_start)

            return dY, Y, p
        end

        # check expected dT due to upwelling longwave flux only
        # (zero conduction when T_base == initial T_bulk)
        T_base_eq = Models.IceSlabParameters{FT}(coupled_param_dict).T_freeze - FT(5.0)
        dY, Y, p =
            test_sea_ice_rhs(coupled_param_dict, thermo_params, space, T_base = T_base_eq)
        σ = coupled_param_dict["stefan_boltzmann_constant"]
        dT_expected =
            (-p.params.ϵ * σ * T_base_eq^4) / (p.params.h * p.params.ρ * p.params.c)
        @test all(T -> T == FT(0) || T ≈ dT_expected, Array(parent(dY.T_bulk)))

        # check expected dT due to downwelling shortwave flux and upwelling longwave flux
        # (again set T_base == initial T_bulk)
        dY, Y, p = test_sea_ice_rhs(
            coupled_param_dict,
            thermo_params,
            space,
            SW_d = 1.0,
            LW_d = 0.0,
            T_base = T_base_eq,
        )
        dT_expected =
            ((1 - p.params.α_cold) * 1.0 - p.params.ϵ * σ * T_base_eq^4) /
            (p.params.h * p.params.ρ * p.params.c)
        @test all(T -> T == FT(0) || T ≈ dT_expected, Array(parent(dY.T_bulk)))

        # check expected dT due to downwelling and upwelling longwave flux
        # (again set T_base == initial T_bulk)
        dY, Y, p = test_sea_ice_rhs(
            coupled_param_dict,
            thermo_params,
            space,
            SW_d = 0.0,
            LW_d = 2.0,
            T_base = T_base_eq,
        )
        dT_expected =
            (p.params.ϵ * (2.0 - σ * T_base_eq^4)) / (p.params.h * p.params.ρ * p.params.c)
        @test all(T -> T == FT(0) || T ≈ dT_expected, Array(parent(dY.T_bulk)))

        # check that tendency will not result in above freezing T
        dY, Y, p = test_sea_ice_rhs(
            coupled_param_dict,
            thermo_params,
            space,
            SW_d = 0.0,
            LW_d = 0.0,
            T_base = 330.0,
        ) # Float32 requires a large number here!
        dT_maximum = @. (p.params.T_freeze - Y.T_bulk) / p.dt
        @test minimum(dT_maximum .- dY.T_bulk) >= FT(0.0)

        # check that the correct tendency was added due to basal conductive flux and upwelling longwave flux
        T_base = 269.2
        dY, Y, p = test_sea_ice_rhs(
            coupled_param_dict,
            thermo_params,
            space,
            SW_d = 0.0,
            LW_d = 0.0,
            T_base = T_base,
        )
        T_bulk = minimum(Y.T_bulk) # get the non-zero temperature value
        (; k_ice, h, ρ, c, T_base, ϵ, T_sfc_min, T_freeze) = p.params
        T_sfc = Models.ice_surface_temperature(T_bulk, T_base, T_sfc_min, T_freeze)
        dT_expected =
            (k_ice / (h * h * ρ * c)) * (T_base - T_sfc) - (ϵ * σ * T_sfc^4) / (h * ρ * c)
        @test minimum(dY) ≈ FT(dT_expected)
        @test maximum(dY) ≈ FT(0)

        # Zero-allocation check for ice_rhs! via @noinline wrapper
        @noinline call_ice_rhs!(dY, Y, p, t) = (Models.ice_rhs!(dY, Y, p, t); nothing)
        call_ice_rhs!(dY, Y, p, 0.0)
        @test (@allocated call_ice_rhs!(dY, Y, p, 0.0)) == 0
    end

    @testset "ice_surface_albedo ramp and coupling-step constancy for FT=$FT" begin
        coupled_param_dict = CP.create_toml_dict(FT)
        thermo_params = TDP.ThermodynamicsParameters(coupled_param_dict)
        params = Models.IceSlabParameters{FT}(coupled_param_dict)
        (; T_freeze, ΔT_melt, α_cold, α_melt) = params

        @test α_cold == FT(0.78)
        @test α_melt == FT(0.55)
        @test ΔT_melt == FT(3)

        # Exactness at both ramp endpoints (including Float32)
        @test Models.ice_surface_albedo(T_freeze, T_freeze, ΔT_melt, α_cold, α_melt) ===
              α_melt
        @test Models.ice_surface_albedo(
            T_freeze + FT(2),
            T_freeze,
            ΔT_melt,
            α_cold,
            α_melt,
        ) === α_melt
        @test Models.ice_surface_albedo(
            T_freeze - ΔT_melt,
            T_freeze,
            ΔT_melt,
            α_cold,
            α_melt,
        ) === α_cold
        @test Models.ice_surface_albedo(
            T_freeze - ΔT_melt - FT(15),
            T_freeze,
            ΔT_melt,
            α_cold,
            α_melt,
        ) === α_cold
        @test Models.ice_surface_albedo(
            T_freeze - ΔT_melt / FT(2),
            T_freeze,
            ΔT_melt,
            α_cold,
            α_melt,
        ) ≈ FT(0.5) * (α_cold + α_melt)

        space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = coupled_param_dict["planet_radius"],
            n_quad_points = 4,
            h_elem = 4,
        )
        sim = Models.PrescribedIceSimulation(
            FT;
            tspan = (0.0, 200.0),
            dt = 50.0,
            saveat = [100.0],
            boundary_space = space,
            coupled_param_dict = coupled_param_dict,
            thermo_params = thermo_params,
            start_date = Dates.DateTime(2010, 1, 1),
            land_fraction = CC.Fields.zeros(space),
            comms_ctx = nothing,
        )

        # Direct and diffuse albedos both alias cached p.α_ice
        @test Interfacer.get_field(sim, Val(:surface_direct_albedo)) ===
              sim.integrator.p.α_ice
        @test Interfacer.get_field(sim, Val(:surface_diffuse_albedo)) ===
              sim.integrator.p.α_ice

        # Place T_bulk in the interior of the melt ramp and apply SW heating
        sim.integrator.u.T_bulk .= T_freeze - FT(1.5)
        Models.update_sea_ice_albedo!(sim)
        sim.integrator.p.SW_d .= FT(400)
        α_before_step = copy(Array(parent(sim.integrator.p.α_ice)))
        T_before_step = copy(Array(parent(sim.integrator.u.T_bulk)))

        # Stepping the ODE integrator directly runs all RK4 stages; p.α_ice must remain invariant
        CTS.step!(sim.integrator, 50.0, true)
        @test Array(parent(sim.integrator.u.T_bulk)) != T_before_step
        @test Array(parent(sim.integrator.p.α_ice)) == α_before_step

        # Stepping through Interfacer.step! refreshes p.α_ice at the coupling-step boundary
        Interfacer.step!(sim, 100.0)
        @test Array(parent(sim.integrator.p.α_ice)) != α_before_step
        expected_α = Models.ice_surface_albedo.(
            Interfacer.get_field(sim, Val(:surface_temperature)),
            T_freeze,
            ΔT_melt,
            α_cold,
            α_melt,
        )
        @test Array(parent(sim.integrator.p.α_ice)) ≈ Array(parent(expected_α))
    end

    @testset "dss_state! AbstractSeaIceSimulation for FT=$FT" begin
        coupled_param_dict = CP.create_toml_dict(FT)
        boundary_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = coupled_param_dict["planet_radius"], # in meters
            n_quad_points = 4,
            h_elem = 4,
        )

        # set up objects for test
        u = CC.Fields.FieldVector(;
            state_field1 = CC.Fields.ones(boundary_space),
            state_field2 = CC.Fields.zeros(boundary_space),
        )
        p = (;
            cache_field = CC.Fields.zeros(boundary_space),
            dss_buffer = Utilities.init_dss_buffer(u),
            domain_type = "global",
        )
        integrator = (; u, p)
        sim = Models.PrescribedIceSimulation(nothing, integrator)

        # make field non-constant to check the impact of the dss step
        coords_lat = CC.Fields.coordinate_field(sim.integrator.u.state_field2).lat
        @. sim.integrator.u.state_field2 = sin(coords_lat)

        # apply DSS
        integrator_copy = deepcopy(integrator)
        Models.dss_state!(sim)

        # test that uniform field and cache are unchanged, non-constant is changed
        # note: uniform field is changed slightly by dss
        @test sim.integrator.u.state_field1 ≈ integrator_copy.u.state_field1
        @test sim.integrator.u.state_field2 != integrator_copy.u.state_field2
        @test sim.integrator.p.cache_field == integrator_copy.p.cache_field
    end
end
