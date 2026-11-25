import Test: @test, @testset
import Dates
import NCDatasets
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC
import Thermodynamics.Parameters as TDP
import ClimaParams as CP # required for TDP
import ClimaCoupler: Interfacer, Models, Utilities
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

            cache = (;
                F_turb_energy = CC.Fields.zeros(space),
                SW_d = CC.Fields.zeros(space) .+ FT(SW_d),
                LW_d = CC.Fields.zeros(space) .+ FT(LW_d),
                α_sfc = CC.Fields.zeros(space),
                h_ice = Models.ice_thickness.(
                    CC.Fields.coordinate_field(space).lat,
                    params,
                ),
                area_fraction = SIC_init,
                SIC_timevaryinginput = SIC_timevaryinginput,
                land_fraction = CC.Fields.zeros(space),
                thermo_params = thermo_params,
                dt = dt,
                binary_area_fraction = true,
                domain_type = "global",
            )

            p = (; cache..., params = params)

            Models.update_ice_albedo!(Y, p)
            Models.ice_rhs!(dY, Y, p, t_start)

            return dY, Y, p
        end

        # check expected dT due to upwelling longwave flux only
        # (zero conduction when T_base == initial T_bulk)
        T_base_eq = Models.IceSlabParameters{FT}(coupled_param_dict).T_freeze - FT(5.0)
        dY, Y, p =
            test_sea_ice_rhs(coupled_param_dict, thermo_params, space, T_base = T_base_eq)
        σ = coupled_param_dict["stefan_boltzmann_constant"]
        dT_nh = (-p.params.ϵ * σ * T_base_eq^4) / (p.params.h_nh * p.params.ρ * p.params.c)
        dT_sh = (-p.params.ϵ * σ * T_base_eq^4) / (p.params.h_sh * p.params.ρ * p.params.c)
        @test all(T -> T == FT(0) || T ≈ dT_nh || T ≈ dT_sh, Array(parent(dY.T_bulk)))

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
        dT_nh =
            ((1 - p.params.α_cold) * 1.0 - p.params.ϵ * σ * T_base_eq^4) /
            (p.params.h_nh * p.params.ρ * p.params.c)
        dT_sh =
            ((1 - p.params.α_cold) * 1.0 - p.params.ϵ * σ * T_base_eq^4) /
            (p.params.h_sh * p.params.ρ * p.params.c)
        @test all(T -> T == FT(0) || T ≈ dT_nh || T ≈ dT_sh, Array(parent(dY.T_bulk)))

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
        dT_nh =
            (p.params.ϵ * (2.0 - σ * T_base_eq^4)) /
            (p.params.h_nh * p.params.ρ * p.params.c)
        dT_sh =
            (p.params.ϵ * (2.0 - σ * T_base_eq^4)) /
            (p.params.h_sh * p.params.ρ * p.params.c)
        @test all(T -> T == FT(0) || T ≈ dT_nh || T ≈ dT_sh, Array(parent(dY.T_bulk)))

        # check that the cap keeps T_sfc from exceeding T_freeze: start T_sfc 0.1 K below
        # T_freeze with strong warming over a long step (uncapped, T_sfc would rise ~1 K),
        # so T_sfc after one step lands exactly on T_freeze only if the cap is active
        params = Models.IceSlabParameters{FT}(coupled_param_dict)
        dt = FT(3600)
        T_bulk = (params.T_base + params.T_freeze) / 2 - FT(0.05)
        dT = Models.ice_energy_balance_tendency(
            T_bulk,
            one(FT),
            zero(FT),
            FT(1000),
            FT(400),
            params.α_melt_nh,
            params.h_nh,
            params,
            dt,
        )
        T_sfc_next = 2 * (T_bulk + dT * dt) - params.T_base # unclamped extrapolation
        @test isapprox(T_sfc_next, params.T_freeze; atol = FT(1e-3))

        # check that the correct tendency was added due to basal conductive flux and upwelling longwave flux
        T_base = 271.2
        dY, Y, p = test_sea_ice_rhs(
            coupled_param_dict,
            thermo_params,
            space,
            SW_d = 0.0,
            LW_d = 0.0,
            T_base = T_base,
        )
        T_bulk = minimum(Y.T_bulk) # get the non-zero temperature value
        (;
            k_ice,
            h_nh,
            h_sh,
            ρ,
            c,
            T_base,
            ϵ,
            T_sfc_min,
            T_freeze,
            ΔT_melt,
            α_cold,
            α_melt_nh,
            α_melt_sh,
        ) = p.params
        T_sfc = Models.ice_surface_temperature(T_bulk, T_base, T_sfc_min, T_freeze)
        dT_expected_nh =
            (k_ice / (h_nh * h_nh * ρ * c)) * (T_base - T_sfc) -
            (ϵ * σ * T_sfc^4) / (h_nh * ρ * c)
        dT_expected_sh =
            (k_ice / (h_sh * h_sh * ρ * c)) * (T_base - T_sfc) -
            (ϵ * σ * T_sfc^4) / (h_sh * ρ * c)
        @test minimum(dY) ≈ FT(min(dT_expected_nh, dT_expected_sh))
        @test maximum(dY) ≈ FT(0)

        # check temperature-dependent albedo ramp and hemispheric helpers
        @test Models.ice_surface_albedo(T_freeze, T_freeze, ΔT_melt, α_cold, α_melt_nh) ==
              α_melt_nh
        @test Models.ice_surface_albedo(
            T_freeze - ΔT_melt,
            T_freeze,
            ΔT_melt,
            α_cold,
            α_melt_nh,
        ) == α_cold
        @test Models.ice_surface_albedo(
            T_freeze - ΔT_melt / FT(2),
            T_freeze,
            ΔT_melt,
            α_cold,
            α_melt_nh,
        ) ≈ FT(0.5) * (α_cold + α_melt_nh)
        @test Models.ice_thickness(FT(70), p.params) == h_nh
        @test Models.ice_thickness(FT(-65), p.params) == h_sh
        @test Models.ice_melt_albedo(FT(70), p.params) == α_melt_nh
        @test Models.ice_melt_albedo(FT(-65), p.params) == α_melt_sh

        # check that `update_ice_albedo!` applies the hemispheric melt albedo pointwise:
        # with T_sfc capped at T_freeze the albedo is exactly α_melt_nh (lat ≥ 0) or
        # α_melt_sh (lat < 0)
        Y_melt = copy(Y)
        Y_melt.T_bulk .= T_freeze + FT(10)
        Models.update_ice_albedo!(Y_melt, p)
        lats = Array(parent(CC.Fields.coordinate_field(space).lat))
        α_sfc = Array(parent(p.α_sfc))
        @test any(lats .< 0) && any(lats .>= 0)
        @test all(α_sfc[lats .>= 0] .== α_melt_nh)
        @test all(α_sfc[lats .< 0] .== α_melt_sh)
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

@testset "PrescribedIceSimulation cached coupling-step albedo" begin
    FT = Float32
    coupled_param_dict = CP.create_toml_dict(FT)
    comms_ctx = ClimaComms.context()
    space = Utilities._column_boundary_space(FT, (FT(75), FT(0)), comms_ctx)
    sim = Models.PrescribedIceSimulation(
        FT;
        tspan = (0.0, 1200.0),
        dt = 600.0,
        saveat = [0.0, 1200.0],
        boundary_space = space,
        coupled_param_dict,
        thermo_params = TDP.ThermodynamicsParameters(coupled_param_dict),
        comms_ctx,
        start_date = Dates.DateTime(2010, 1, 1),
        land_fraction = CC.Fields.zeros(space),
        domain_type = "column",
    )

    # Mid-ramp T_sfc with warming forcing: RK4 stages keep `p.α_sfc` fixed,
    # while `Interfacer.step!` refreshes `p.α_sfc` at the coupling step.
    sim.integrator.u.T_bulk .=
        (sim.params.T_freeze - sim.params.ΔT_melt / FT(2) + sim.params.T_base) / FT(2)
    Models.update_ice_albedo!(sim.integrator.u, sim.integrator.p)
    α_0 = Array(parent(sim.integrator.p.α_sfc))[1]
    Interfacer.update_field!(sim, Val(:SW_d), CC.Fields.ones(space) .* FT(300))
    Interfacer.update_field!(sim, Val(:LW_d), CC.Fields.ones(space) .* FT(320))

    Interfacer.step!(sim.integrator)
    @test Array(parent(sim.integrator.p.α_sfc))[1] == α_0

    Interfacer.step!(sim, 1200.0)
    T_sfc = Array(parent(Interfacer.get_field(sim, Val(:surface_temperature))))[1]
    (; T_freeze, ΔT_melt, α_cold, α_melt_nh) = sim.params
    @test Array(parent(sim.integrator.p.α_sfc))[1] ≈
          Models.ice_surface_albedo(T_sfc, T_freeze, ΔT_melt, α_cold, α_melt_nh)
    @test Array(parent(sim.integrator.p.α_sfc))[1] < α_0
end

@testset "PrescribedIceSimulation ice-free reset to T_base and SIC sanitization" begin
    FT = Float32
    coupled_param_dict = CP.create_toml_dict(FT)
    comms_ctx = ClimaComms.context()
    space = CC.CommonSpaces.CubedSphereSpace(
        FT;
        radius = coupled_param_dict["planet_radius"],
        n_quad_points = 4,
        h_elem = 4,
    )
    sim = Models.PrescribedIceSimulation(
        FT;
        tspan = (0.0, 86400.0),
        dt = 600.0,
        saveat = [0.0, 86400.0],
        boundary_space = space,
        coupled_param_dict,
        thermo_params = TDP.ThermodynamicsParameters(coupled_param_dict),
        comms_ctx,
        # Start inside the time range of the lowres SIC artifact (1870-1871) so the SIC
        # varies seasonally and new ice forms between the two reads below.
        start_date = Dates.DateTime(1870, 2, 1),
        land_fraction = CC.Fields.zeros(space),
        binary_area_fraction = false,
        domain_type = "global",
    )

    # Ice-free cells initialize at T_base, reset to T_base on read_sic_data! with
    # immediate land clamping, and stay at T_base when new seasonal ice forms (Aug 1).
    # `≈` rather than `==`: the constructor DSSes the state, and the weighted average at
    # shared element nodes need not reproduce `T_base` bit-for-bit.
    ice_free = Array(parent(sim.integrator.p.area_fraction)) .<= zero(FT)
    @test any(ice_free) &&
          all(T -> T ≈ sim.params.T_base, Array(parent(sim.integrator.u.T_bulk))[ice_free])

    T_stale = sim.params.T_base - FT(15)
    sim.integrator.u.T_bulk .= T_stale
    sim.integrator.p.land_fraction .= FT(0.25)
    Models.read_sic_data!(sim.integrator)
    sic_feb = Array(parent(sim.integrator.p.area_fraction))
    T_feb = Array(parent(sim.integrator.u.T_bulk))
    @test all(zero(FT) .<= sic_feb .<= FT(0.75))
    @test all(T_feb[sic_feb .<= zero(FT)] .== sim.params.T_base)
    @test all(T_feb[sic_feb .> zero(FT)] .== T_stale)

    sim.integrator.p.land_fraction .= zero(FT)
    sim.integrator.t = 181.0 * 86400.0
    Models.read_sic_data!(sim.integrator)
    sic_aug = Array(parent(sim.integrator.p.area_fraction))
    newly_formed = (sic_feb .<= zero(FT)) .& (sic_aug .> zero(FT))
    @test any(newly_formed) &&
          all(Array(parent(sim.integrator.u.T_bulk))[newly_formed] .== sim.params.T_base)
end
