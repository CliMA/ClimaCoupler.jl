import Test: @test, @testset
import ClimaCore as CC
import ClimaParams as CP # required for coupled_param_dict
import ClimaCoupler: Interfacer, Models

for FT in (Float32, Float64)
    coupled_param_dict = CP.create_toml_dict(FT)
    params_ice = Models.EisenmanIceParameters{FT}(coupled_param_dict)
    params_ocean = Models.EisenmanOceanParameters{FT}()
    params = (; p_i = params_ice, p_o = params_ocean)
    (; T_base, T_freeze, L_ice, k_ice, C0_base, σ, ϵ) = params_ice
    hρc_ml = params_ocean.h * params_ocean.ρ * params_ocean.c

    Δt = Float64(1e6)

    # create a boundary space
    boundary_space = CC.CommonSpaces.CubedSphereSpace(
        FT;
        radius = coupled_param_dict["planet_radius"], # in meters
        n_quad_points = 4,
        h_elem = 4,
    )

    @testset "No net fluxes for FT=$FT" begin
        Y, Ya = Models.eisenman_state_init(params_ice, boundary_space)

        # ice covered
        Y.h_ice .= 1

        # no conductive flux, or basal interface flux
        Y.T_sfc .= T_base
        Y.T_ml .= T_base

        # zero net atmospheric flux: downwelling longwave balances surface emission
        @. Ya.LW_d = σ * Y.T_sfc^4

        Models.solve_eisenman_model!(Y, Ya, params, Δt)

        @test all(parent(Y.e_base) .≈ 0)
        @test all(parent(Y.T_ml) .≈ T_base)
        @test all(parent(Y.T_sfc) .≈ T_base)
        @test all(parent(Y.h_ice) .≈ 1)
    end

    @testset "Radiative fluxes: outgoing longwave -> ice growth for FT=$FT" begin
        Y, Ya = Models.eisenman_state_init(params_ice, boundary_space)

        # ice covered
        Y.h_ice .= 1
        h_ice_0 = deepcopy(Y.h_ice)

        # no conductive flux, or basal interface flux
        Y.T_sfc .= T_base
        Y.T_ml .= T_base

        # ice growth due to outgoing longwave (no downwelling radiation, no
        # turbulent fluxes: F_atm = ϵ σ T_sfc^4)
        F_atm = @. ϵ * σ * Y.T_sfc^4
        ∂F_atm∂T_sfc = @. 4 * ϵ * σ * Y.T_sfc^3

        Models.solve_eisenman_model!(Y, Ya, params, Δt)

        @test all(parent(Y.e_base) .≈ 0)
        @test all(parent(Y.T_ml) .≈ T_base) # ocean temperature below ice remains unchanged
        h_ice_new = @. h_ice_0 + F_atm * FT(Δt) / L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new)) # ice growth
        T_sfc_new = @. T_base + (-F_atm) / (k_ice / h_ice_new + ∂F_atm∂T_sfc)
        @test all(parent(Y.T_sfc) .≈ parent(T_sfc_new)) # surface temperature decreases
    end

    @testset "Radiative fluxes balance turbulent fluxes -> no ice growth for FT=$FT" begin
        # no net atmospheric flux: net radiative loss balanced by turbulent gain
        Y, Ya = Models.eisenman_state_init(params_ice, boundary_space)

        # ice covered
        Y.h_ice .= 1
        h_ice_0 = deepcopy(Y.h_ice)

        # no conductive flux, or basal interface flux
        Y.T_sfc .= T_base
        Y.T_ml .= T_base

        # prescribe radiative and turbulent fluxes: net radiative flux 100 W/m2
        # upward, turbulent flux 100 W/m2 downward
        @. Ya.LW_d = σ * Y.T_sfc^4 - 100 / ϵ
        Ya.F_turb .= -100

        Models.solve_eisenman_model!(Y, Ya, params, Δt)

        @test all(parent(Y.e_base) .≈ 0) # no contribution from basal fluxes
        @test all(parent(Y.T_ml) .≈ T_base) # ocean temperature stays at freezing point
        @test all(parent(Y.h_ice) .≈ parent(h_ice_0)) # no ice growth
        @test all(parent(Y.T_sfc) .≈ T_base) # surface temperature unchanged
    end

    @testset "Ice surface temp below freezing for FT=$FT" begin
        # Ensure ice surface temperature never exceeds the freezing point (even
        # when positive radiative and conductive fluxes)
        Y, Ya = Models.eisenman_state_init(params_ice, boundary_space)

        # ice covered
        Y.h_ice .= 10
        h_ice_0 = deepcopy(Y.h_ice)

        # no basal interface fluxes
        Y.T_ml .= T_base

        # non-zero conductive flux
        Y.T_sfc .= T_base .- 1
        T_sfc_0 = deepcopy(Y.T_sfc)

        # net incoming energy: downwelling longwave exceeds surface emission by 300 W/m2
        Ya.F_turb .= 0
        @. Ya.LW_d = σ * Y.T_sfc^4 + 300 / ϵ

        Models.solve_eisenman_model!(Y, Ya, params, Δt)

        F_atm = -300
        ∂F_atm∂T_sfc = @. 4 * ϵ * σ * T_sfc_0^3

        @test all(parent(Y.e_base) .≈ 0) # no contribution from basal fluxes
        @test all(parent(Y.T_ml) .≈ T_base) # ocean temperature stays at freezing point
        h_ice_new = @. h_ice_0 + F_atm * FT(Δt) / L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new)) # ice melt
        δT_sfc = @. (-F_atm + k_ice / h_ice_new * (T_base - T_sfc_0)) /
           (k_ice / h_ice_new + ∂F_atm∂T_sfc)
        @test all(parent(Y.T_sfc) .≈ T_freeze) # ice surface temperature doesn't exceed freezing
        @test all(parent(Y.T_sfc) .< parent(T_sfc_0 .+ δT_sfc))
    end

    @testset "Transition to ice free due to incoming > outgoing for FT=$FT" begin
        Y, Ya = Models.eisenman_state_init(params_ice, boundary_space)

        # ice covered
        Y.h_ice .= 0.1
        h_ice_0 = deepcopy(Y.h_ice)

        # no conductive flux, or basal interface flux
        Y.T_sfc .= T_base
        Y.T_ml .= T_base
        T_ml_0 = deepcopy(Y.T_ml)

        # net incoming energy: downwelling longwave exceeds surface emission by 300 W/m2
        Ya.F_turb .= 0
        @. Ya.LW_d = σ * Y.T_sfc^4 + 300 / ϵ

        Models.solve_eisenman_model!(Y, Ya, params, Δt)

        F_atm = -300

        @test all(parent(Y.h_ice) .≈ 0) # ice melts
        T_ml_new = @. T_ml_0 - F_atm * FT(Δt) / hρc_ml - h_ice_0 * L_ice / hρc_ml
        # ocean temperature increases due to F_atm (reduced by the latent heat of melting)
        @test all(parent(Y.T_ml) .≈ parent(T_ml_new))
        @test all(parent(Y.T_sfc) .≈ parent(T_ml_new)) # ice-free: T_sfc = T_ml
    end

    @testset "Mixed layer warming due to net incoming radiation fluxes for FT=$FT" begin
        Y, Ya = Models.eisenman_state_init(params_ice, boundary_space)

        # ice free
        Y.h_ice .= 0

        # no conductive flux, or basal interface flux
        Y.T_sfc .= T_base
        Y.T_ml .= T_base
        T_ml_0 = deepcopy(Y.T_ml)

        # net incoming energy: downwelling longwave exceeds surface emission by 300 W/m2
        Ya.F_turb .= 0
        @. Ya.LW_d = σ * Y.T_sfc^4 + 300 / ϵ

        Models.solve_eisenman_model!(Y, Ya, params, Δt)

        F_atm = -300

        @test all(parent(Y.e_base) .≈ 0) # no contribution from basal fluxes
        T_ml_new = @. T_ml_0 - F_atm * FT(Δt) / hρc_ml
        @test all(parent(Y.T_ml) .≈ parent(T_ml_new)) # ocean temperature increases
        @test all(parent(Y.h_ice) .≈ 0) # no ice
        @test all(parent(Y.T_sfc) .≈ parent(T_ml_new)) # surface temperature = ocean temperature
    end

    @testset "Mixed layer freezing due to net outgoing radiation fluxes for FT=$FT" begin
        Y, Ya = Models.eisenman_state_init(params_ice, boundary_space)

        # ice free, mixed layer at the freezing point
        Y.h_ice .= 0
        Y.T_sfc .= T_freeze
        Y.T_ml .= T_freeze
        T_sfc_0 = deepcopy(Y.T_sfc)
        T_ml_0 = deepcopy(Y.T_ml)

        # net outgoing longwave (no downwelling radiation, no turbulent fluxes)
        F_atm = @. ϵ * σ * Y.T_sfc^4
        ∂F_atm∂T_sfc = @. 4 * ϵ * σ * Y.T_sfc^3

        Models.solve_eisenman_model!(Y, Ya, params, Δt)

        @test all(parent(Y.e_base) .≈ 0) # no contribution from basal fluxes
        @test all(parent(Y.T_ml) .≈ T_freeze) # ocean temperature remains at the freezing point
        h_ice_new = @. (F_atm - (T_ml_0 - T_freeze) * hρc_ml / FT(Δt)) * FT(Δt) / L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new)) # frazil ice growth
        @test all(parent(Y.h_ice) .> 0)
        T_sfc_new = @. T_sfc_0 - F_atm / (k_ice / h_ice_new + ∂F_atm∂T_sfc)
        @test all(parent(Y.T_sfc) .≈ parent(T_sfc_new)) # surface temperature decreases
    end

    # conductive flux and q flux
    @testset "Non-zero conductive flux for FT=$FT" begin
        Y, Ya = Models.eisenman_state_init(params_ice, boundary_space)
        Δt_conductive = Float64(100)

        # ice covered
        Y.h_ice .= 10
        h_ice_0 = deepcopy(Y.h_ice)

        # non-zero interface basal flux
        ΔT_ml = 10
        Y.T_ml .= T_base .+ ΔT_ml
        T_ml_0 = deepcopy(Y.T_ml)

        # zero conductive fluxes
        Y.T_sfc .= T_base

        # zero net atmospheric flux
        Ya.F_turb .= 0
        @. Ya.LW_d = σ * Y.T_sfc^4

        Models.solve_eisenman_model!(Y, Ya, params, Δt_conductive)

        # non-zero contribution from basal flux
        @test all(parent(Y.e_base) .≈ C0_base * ΔT_ml * FT(Δt_conductive))
        T_ml_new = @. T_ml_0 - C0_base * ΔT_ml * FT(Δt_conductive) / hρc_ml
        @test all(parent(Y.T_ml) .≈ parent(T_ml_new)) # ocean temperature decreases
        @test all(parent(Y.T_ml) .< parent(T_ml_0))
        h_ice_new = @. h_ice_0 - C0_base * ΔT_ml * FT(Δt_conductive) / L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new)) # basal flux melts ice
        # surface temperature unchanged (T_ml can affect it only if ice free)
        @test all(parent(Y.T_sfc) .≈ T_base)
    end

    @testset "Nonzero Q-flux (~horizontal ocean transport) for FT=$FT" begin
        Y, Ya = Models.eisenman_state_init(params_ice, boundary_space)
        Δt_qflux = Float64(100)

        # ice covered
        Y.h_ice .= 10
        h_ice_0 = deepcopy(Y.h_ice)

        # positive q_flux (positive toward current column)
        Ya.ocean_qflux .= 100

        # no conductive flux, or basal interface flux
        Y.T_sfc .= T_base
        Y.T_ml .= T_base
        T_ml_0 = deepcopy(Y.T_ml)

        # zero net atmospheric flux
        Ya.F_turb .= 0
        @. Ya.LW_d = σ * Y.T_sfc^4

        Models.solve_eisenman_model!(Y, Ya, params, Δt_qflux)

        @test all(parent(Y.e_base) .≈ 0) # no contribution from basal flux
        T_ml_new = @. T_ml_0 + Ya.ocean_qflux * FT(Δt_qflux) / hρc_ml
        @test all(parent(Y.T_ml) .≈ parent(T_ml_new)) # qflux increases ocean temperature
        @test all(parent(Y.T_ml) .> parent(T_ml_0))
        # surface temperature unchanged (T_ml can affect it only if ice free)
        @test all(parent(Y.T_sfc) .≈ T_base)
        h_ice_new = @. h_ice_0 - Ya.ocean_qflux * FT(Δt_qflux) / L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new)) # qflux decreases sea-ice thickness
        @test all(parent(Y.h_ice) .< parent(h_ice_0))
    end

    @testset "timestep update + total energy calculation for FT=$FT" begin
        Δt_sim = Float64(1000)

        sim = Interfacer.SeaIceSimulation(
            FT,
            Val(:eisenman);
            tspan = (0.0, 2e6),
            dt = Δt_sim,
            saveat = [2e6],
            boundary_space = boundary_space,
            coupled_param_dict = coupled_param_dict,
        )
        cache = sim.integrator.p
        u = sim.integrator.u

        # init conditions: ice-free, ocean at the freezing point, constant net
        # upward atmospheric flux of 300 W/m2
        u.T_ml .= T_freeze
        cache.F_turb .= 0
        @. cache.LW_d = σ * u.T_sfc^4 - 300 / ϵ

        total_energy_0 = deepcopy(Interfacer.get_field(sim, Val(:energy)))

        # Solve for one timestep and check the results (frazil ice formation)
        Models.solve_eisenman_model!(u, cache, cache.params, Δt_sim)
        @test all(parent(u.h_ice) .≈ 300 * Δt_sim / L_ice)

        # Re-prescribe the downwelling longwave so the net flux stays at
        # 300 W/m2 after the Newton update of T_sfc, then solve a second
        # timestep (bottom ice growth)
        @. cache.LW_d = σ * u.T_sfc^4 - 300 / ϵ
        Models.solve_eisenman_model!(u, cache, cache.params, Δt_sim)
        @test all(abs.(parent(u.h_ice) .- 2 * 300 * Δt_sim / L_ice) .< 10eps(FT))

        # the energy gained by the column equals the accumulated atmospheric flux
        total_energy_calc = Interfacer.get_field(sim, Val(:energy)) .- total_energy_0
        total_energy_expected = 300 .* ones(boundary_space) .* 2 .* FT(Δt_sim)
        @test all(parent(total_energy_calc) .≈ parent(total_energy_expected))
    end
end
