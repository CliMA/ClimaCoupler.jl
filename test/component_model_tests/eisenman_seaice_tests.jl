import Test: @test, @testset
import ClimaCore as CC
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import ClimaCoupler
import ClimaCoupler: Interfacer, Regridder, TestHelper

include("../../experiments/ClimaEarth/components/ocean/eisenman_seaice.jl")


for FT in (Float32, Float64)
    params_ice = EisenmanIceParameters{FT}()
    params_ocean = EisenmanOceanParameters{FT}()
    params = (; p_i = params_ice, p_o = params_ocean)

    Δt = Float64(1e6)

    # create a boundary space
    boundary_space = TestHelper.create_space(FT)

    # thermodynamic parameter set
    thermo_params = TP.ThermodynamicsParameters(FT)

    @testset "No net fluxes for FT=$FT" begin
        Y, Ya = state_init(params_ice, boundary_space)

        # ice covered
        Y.h_ice .= 1

        # no conductive flux, or basal interface flux
        Y.T_sfc .= params_ice.T_base
        Y.T_ml .= params_ice.T_base

        # no atmos fluxes, F_a
        Ya.∂F_turb_energy∂T_sfc .= .-get_∂F_rad_energy∂T_sfc.(Y.T_sfc, Ref(params_ice)) # ∂F_turb_energy∂T_sfc + ∂F_rad_energy∂T_sfc = 0
        Ya.F_rad .= 0
        Ya.F_turb .= 0

        CC.Fields.bycolumn(boundary_space) do colidx
            solve_eisenman_model!(Y[colidx], Ya[colidx], params, thermo_params, Δt)
        end
        @test all(parent(Ya.e_base) .≈ 0)
        @test all(parent(Y.T_ml) .≈ params_ice.T_base)
        @test all(parent(Y.T_sfc) .≈ params_ice.T_base)
        @test all(parent(Y.h_ice) .≈ 1)
    end

    @testset "Radiative fluxes: outgoing longwave -> ice growth for FT=$FT" begin
        Y, Ya = state_init(params_ice, boundary_space)

        # ice covered
        Y.h_ice .= 1
        h_ice_0 = deepcopy(Y.h_ice)

        # no conductive flux, or basal interface flux
        Y.T_sfc .= params_ice.T_base
        Y.T_ml .= params_ice.T_base

        # ice growth due to outgoing longwave
        # only radiative (no aerodynamic) atmos fluxes considered here (F_atm = F_rad)
        Ya.∂F_turb_energy∂T_sfc .= 0
        ∂F_atm∂T_sfc = get_∂F_rad_energy∂T_sfc(Y.T_sfc, params_ice) .+ Ya.∂F_turb_energy∂T_sfc
        @. Ya.F_rad = (1 - params_ice.α) * params_ice.σ * Y.T_sfc .^ 4 # outgoing longwave

        CC.Fields.bycolumn(boundary_space) do colidx
            solve_eisenman_model!(Y[colidx], Ya[colidx], params, thermo_params, Δt)
        end

        F_atm = @. Ya.F_rad + Ya.F_turb

        @test all(parent(Ya.e_base) .≈ 0)
        @test all(parent(Y.T_ml) .≈ params_ice.T_base) # ocean temperature below ice remains unchanged
        h_ice_new = @. h_ice_0 + F_atm * FT(Δt) / params_ice.L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new)) # ice growth
        T_sfc_new = @. params_ice.T_base .+ (-F_atm) / (params_ice.k_ice / (h_ice_new) + ∂F_atm∂T_sfc)
        @test all(parent(Y.T_sfc) .≈ parent(T_sfc_new)) # surface temperature decreases
    end

    @testset "Radiative fluxes balance turbulent fluxes -> no ice growth for FT=$FT" begin
        # no ice growth due to outgoing longwave = incoming longwave + incoming shortwave
        # (NB: incoming independent of T_sfc)
        Y, Ya = state_init(params_ice, boundary_space)

        # ice covered
        Y.h_ice .= 1
        h_ice_0 = deepcopy(Y.h_ice)

        # prescribe radiative and tubulent fluxes
        Ya.F_rad .= 100
        Ya.F_turb .= -100

        # no conductive flux, or basal interface flux
        Y.T_sfc .= params_ice.T_base
        Y.T_ml .= params_ice.T_base

        CC.Fields.bycolumn(boundary_space) do colidx
            solve_eisenman_model!(Y[colidx], Ya[colidx], params, thermo_params, Δt)
        end

        F_atm = @. Ya.F_rad + Ya.F_turb
        ∂F_atm∂T_sfc = 0

        @test all(parent(Ya.e_base) .≈ 0) # no contribution from basal fluxes
        @test all(parent(Y.T_ml) .≈ params_ice.T_base) # ocean temperature stays at freezing point
        h_ice_new = @. h_ice_0 + F_atm * FT(Δt) / params_ice.L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new)) # no ice growth
        T_sfc_new = @. params_ice.T_base .+ (-F_atm) / (params_ice.k_ice / (h_ice_new) + ∂F_atm∂T_sfc)
        @test all(parent(Y.T_sfc) .≈ parent(T_sfc_new)) # surface temperature decreases
    end

    @testset "Ice surface temp below freezing for FT=$FT" begin
        # Ensure ice surface temperature never exceeds the freezing point (even when positive radiative and conductive fluxes)

        Y, Ya = state_init(params_ice, boundary_space)

        # ice covered
        Y.h_ice .= 10
        h_ice_0 = deepcopy(Y.h_ice)

        # no basal interface fluxes
        Y.T_ml .= params_ice.T_base

        # net incoming energy
        Ya.F_turb .= 0
        Ya.F_rad .= @. (1 - params_ice.α) * params_ice.σ * Y.T_sfc^4 - 300 # outgoing longwave < incoming longwave + incoming shortwave

        # non-zero conductive flux
        Y.T_sfc .= params_ice.T_base .- 1
        T_sfc_0 = deepcopy(Y.T_sfc)

        CC.Fields.bycolumn(boundary_space) do colidx
            solve_eisenman_model!(Y[colidx], Ya[colidx], params, thermo_params, Δt)
        end

        F_atm = @. Ya.F_rad + Ya.F_turb
        ∂F_atm∂T_sfc = 0

        F_conductive = @. params_ice.k_ice / (Y.h_ice) * (params_ice.T_base - T_sfc_0)
        ΔT_sfc = @. (-F_atm + F_conductive) / (params_ice.k_ice / (Y.h_ice) + ∂F_atm∂T_sfc)

        @test all(parent(Ya.e_base) .≈ 0) # no contribution from basal fluxes
        @test all(parent(Y.T_ml) .≈ params_ice.T_base) # ocean temperature stays at freezing point
        h_ice_new = @. h_ice_0 + F_atm * FT(Δt) / params_ice.L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new)) # ice growth
        @test all(parent(Y.T_sfc) .≈ params_ice.T_base) # ice surface temperature doesn't exceed freezing
        @test all(parent(Y.T_sfc) .< parent(T_sfc_0 .+ ΔT_sfc))
    end

    @testset "Transition to ice free due to incoming > outgoing for FT=$FT" begin
        Y, Ya = state_init(params_ice, boundary_space)

        # ice covered
        Y.h_ice .= 0.1
        h_ice_0 = deepcopy(Y.h_ice)

        # prescribe radiative and tubulent fluxes
        Ya.F_rad .= @. (1 - params_ice.α) * params_ice.σ * Y.T_sfc .^ 4 - 300  # outgoing longwave < incoming longwave + incoming shortwave
        Ya.F_turb .= 0
        ∂F_atm∂T_sfc = 0

        # no conductive flux, or basal interface flux
        Y.T_sfc .= params_ice.T_base
        Y.T_ml .= params_ice.T_base
        T_ml_0 = deepcopy(Y.T_ml)

        CC.Fields.bycolumn(boundary_space) do colidx
            solve_eisenman_model!(Y[colidx], Ya[colidx], params, thermo_params, Δt)
        end

        F_atm = @. Ya.F_rad + Ya.F_turb

        @test all(parent(Y.h_ice) .≈ 0) # ice melts
        T_ml_new = @. T_ml_0 - (F_atm) * FT(Δt) / (params_ocean.h * params_ocean.ρ * params_ocean.c) -
           h_ice_0 * params_ice.L_ice / (params_ocean.h * params_ocean.ρ * params_ocean.c)
        @test all(parent(Y.T_ml) .≈ parent(T_ml_new)) # ocean temperature increases due to F_atm (reduced by the latent heat of melting)

    end

    @testset "Mixed layer warming due to net negative radiation fluxes for FT=$FT" begin
        Y, Ya = state_init(params_ice, boundary_space)

        # ice free
        Y.h_ice .= 0
        h_ice_0 = deepcopy(Y.h_ice)

        # no conductive flux, or basal interface flux
        Y.T_sfc .= params_ice.T_base
        Y.T_ml .= params_ice.T_base
        T_ml_0 = deepcopy(Y.T_ml)

        # net positive radiative fluxes
        Y.T_sfc .= params_ice.T_base
        Ya.F_turb .= 0
        Ya.F_rad .= @. (1 - params_ice.α) * params_ice.σ * Y.T_sfc^4 - 300

        CC.Fields.bycolumn(boundary_space) do colidx
            solve_eisenman_model!(Y[colidx], Ya[colidx], params, thermo_params, Δt)
        end

        F_atm = @. Ya.F_rad + Ya.F_turb

        @test all(parent(Ya.e_base) .≈ 0) # no contribution from basal fluxes
        T_ml_new = @. T_ml_0 - (F_atm) * FT(Δt) / (params_ocean.h * params_ocean.ρ * params_ocean.c)
        @test all(parent(Y.T_ml) .≈ parent(T_ml_new)) # ocean temperature increases
        @test all(parent(Y.h_ice) .≈ 0) # no ice
        @test all(parent(Y.T_sfc) .≈ parent(T_ml_new)) # surface temperature = ocean temperature

    end

    @testset "Mixed layer freezing due to net positive radiation fluxes for FT=$FT" begin
        Y, Ya = state_init(params_ice, boundary_space)

        # ice free
        Y.h_ice .= 0
        h_ice_0 = deepcopy(Y.h_ice)

        # no conductive flux, or basal interface flux
        Y.T_sfc .= params_ice.T_base
        Y.T_ml .= params_ice.T_base
        T_sfc_0 = deepcopy(Y.T_sfc)
        T_ml_0 = deepcopy(Y.T_ml)

        # net outgoing longwave
        Ya.F_rad .= @. (1 - params_ice.α) * params_ice.σ * Y.T_sfc^4

        CC.Fields.bycolumn(boundary_space) do colidx
            solve_eisenman_model!(Y[colidx], Ya[colidx], params, thermo_params, Δt)
        end

        F_atm = @. Ya.F_rad + Ya.F_turb
        ∂F_atm∂T_sfc = get_∂F_rad_energy∂T_sfc(T_sfc_0, params_ice) .+ Ya.∂F_turb_energy∂T_sfc

        FT = eltype(F_atm)
        @test all(parent(Ya.e_base) .≈ 0) # no contribution from basal fluxes
        @test all(parent(Y.T_ml) .≈ params_ice.T_freeze) # ocean temperature remains at the freezing point
        h_ice_new =
            @. (F_atm - (T_ml_0 - params_ice.T_freeze) * params_ocean.h * params_ocean.ρ * params_ocean.c) * FT(Δt) /
               params_ice.L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new)) # ice growth
        @test all(parent(Y.h_ice) .> 0)
        T_sfc_new = @. T_sfc_0 - F_atm / (params_ice.k_ice / (h_ice_new) + ∂F_atm∂T_sfc)
        @test all(parent(Y.T_sfc) .≈ parent(T_sfc_new)) # surface temperature decreases
    end

    # conductive flux and q flux
    @testset "Non-zero conductive flux for FT=$FT" begin
        Y, Ya = state_init(params_ice, boundary_space)
        Δt = Float64(100)

        # ice covered
        Y.h_ice .= 10
        h_ice_0 = deepcopy(Y.h_ice)

        # non-zero interface basal
        ΔT_ml = 10
        Y.T_ml .= params_ice.T_base .+ ΔT_ml
        T_ml_0 = deepcopy(Y.T_ml)

        # zero conductive fluxes
        Y.T_sfc .= params_ice.T_base

        # zero atmos fluxes and their derivatives
        Ya.∂F_turb_energy∂T_sfc .= .-get_∂F_rad_energy∂T_sfc.(Y.T_sfc, Ref(params_ice)) # ∂F_turb_energy∂T_sfc + ∂F_rad_energy∂T_sfc = 0
        Ya.F_turb .= 0
        Ya.F_rad .= 0

        CC.Fields.bycolumn(boundary_space) do colidx
            solve_eisenman_model!(Y[colidx], Ya[colidx], params, thermo_params, Δt)
        end

        @test all(parent(Ya.e_base) .≈ params_ice.C0_base * ΔT_ml * FT(Δt)) # non-zero contribution from basal flux
        T_ml_new = T_ml_0 .- params_ice.C0_base .* ΔT_ml * FT(Δt) / (params_ocean.h * params_ocean.ρ * params_ocean.c)
        @test all(parent(Y.T_ml) .≈ parent(T_ml_new)) # ocean temperature decreases
        @test all(parent(Y.T_ml) .< parent(T_ml_0))
        h_ice_new = @. h_ice_0 - params_ice.C0_base * 10 * FT(Δt) / params_ice.L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new))
        @test all(parent(Y.T_sfc) .≈ params_ice.T_base) # surface temperature unchanged (T_ml can affect in only if ice free)
    end

    @testset "Nonzero Q-flux (~horizontal ocean transport) for FT=$FT" begin
        Y, Ya = state_init(params_ice, boundary_space)
        Δt = Float64(100)

        # ice covered
        Y.h_ice .= 10
        h_ice_0 = deepcopy(Y.h_ice)

        # positive q_flux (positive toward current column)
        Ya.ocean_qflux .= 100

        # no conductive flux, or basal interface flux
        Y.T_sfc .= params_ice.T_base
        Y.T_ml .= params_ice.T_base
        T_ml_0 = deepcopy(Y.T_ml)

        # no atmos fluxes or their derivatives
        Ya.∂F_turb_energy∂T_sfc .= .-get_∂F_rad_energy∂T_sfc.(Y.T_sfc, Ref(params_ice)) # ∂F_turb_energy∂T_sfc + ∂F_rad_energy∂T_sfc = 0
        Ya.F_turb .= 0
        Ya.F_rad .= 0

        CC.Fields.bycolumn(boundary_space) do colidx
            solve_eisenman_model!(Y[colidx], Ya[colidx], params, thermo_params, Δt)
        end

        @test all(parent(Ya.e_base) .≈ 0) # no contribution from basal flux
        T_ml_new = @. T_ml_0 + Ya.ocean_qflux * FT(Δt) / (params_ocean.h * params_ocean.ρ * params_ocean.c)
        @test all(parent(Y.T_ml) .≈ parent(T_ml_new)) # qflux increases ocean temperature
        @test all(parent(Y.T_ml) .> parent(T_ml_0))
        @test all(parent(Y.T_sfc) .≈ params_ice.T_base) # surface temperature unchanged (T_ml can affect in only if ice free)
        h_ice_new = @. h_ice_0 - Ya.ocean_qflux * FT(Δt) / params_ice.L_ice
        @test all(parent(Y.h_ice) .≈ parent(h_ice_new)) # # qflux decreases sea-ice thickness
        @test all(parent(Y.h_ice) .< parent(h_ice_0))
    end

    include("../../experiments/ClimaEarth/components/slab_utils.jl")
    @testset "step! update + total energy calculation for FT=$FT" begin
        Δt = Float64(1000)

        sim = eisenman_seaice_init(
            FT,
            (0, 2e6),
            space = boundary_space,
            area_fraction = ones(boundary_space),
            thermo_params = thermo_params,
            stepper = CTS.RK4(),
            dt = Δt,
            saveat = 1.0e10,
        )
        sim.integrator.p.Ya.F_turb .= 0
        sim.integrator.p.Ya.F_rad .= 300
        sim.integrator.u.T_ml .= sim.integrator.p.params.p_i.T_freeze # init conditions for ocean temperature

        total_energy_0 = Interfacer.get_field(sim, Val(:energy))
        h_ice_0 = deepcopy(sim.integrator.u.h_ice)

        Interfacer.step!(sim, Δt)
        h_ice = sim.integrator.u.h_ice
        @test all(parent(h_ice) .≈ 0.001)
        Interfacer.step!(sim, 2 * Δt)
        h_ice = sim.integrator.u.h_ice
        @test all(abs.(parent(h_ice) .- 0.002) .< 10eps(FT))

        total_energy_calc = (Interfacer.get_field(sim, Val(:energy)) .- total_energy_0)
        total_energy_expeted = 300 .* ones(boundary_space) .* 2 .* FT(Δt)
        @test all(parent(total_energy_calc) .≈ parent(total_energy_expeted))

    end
end
