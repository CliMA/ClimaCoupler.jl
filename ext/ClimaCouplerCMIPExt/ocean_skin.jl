"""
    Ocean skin temperature (DiffusiveFlux) following ClimaOcean's SkinTemperature{<:DiffusiveFlux}.

Single-equation flux balance: diffusive flux from ocean into the skin layer equals
the total flux out of the surface (turbulent + net radiative). So
  (ρ*c)*κ*(T_ocean - T_sfc)/δ = Q_turb + εσT_sfc^4 - (1-α)*SW_d - LW_d
with κ = thermal diffusivity (m²/s), δ = skin layer thickness (m), ρ = reference_density, c = heat_capacity.
All fluxes positive upward (W/m²).

When used with SurfaceFluxes.jl's `update_T_sfc` callback (via `ocean_update_T_sfc_callback`),
the skin temperature is solved inside the Monin-Obukhov iteration: at each ζ step the
callback returns T_sfc that satisfies the ocean flux balance with Q_turb from the
current MOST state (ζ, u_star, z0m, z0h).
"""

using RootSolvers: find_zero, BrentMethod, CompactSolution
import SurfaceFluxes as SF
import Thermodynamics as TD

"""
    ocean_flux_balance_residual(T_sfc, Q_turb, SW_d, LW_d, α, ε, σ, T_ocean, δ, κ, ρ, c)

Flux balance residual for ocean skin (DiffusiveFlux). Zero when
  (ρ*c)*κ*(T_ocean - T_sfc)/δ = Q_turb + εσT_sfc^4 - (1-α)*SW_d - LW_d.

All arguments are scalars. ρ = reference_density (kg/m³), c = heat_capacity (J/(kg*K)).
Returns residual in W/m².
"""
function ocean_flux_balance_residual(
    T_sfc,
    Q_turb,
    SW_d,
    LW_d,
    α,
    ε,
    σ,
    T_ocean,
    δ,
    κ,
    ρ,
    c,
)
    diffusive_flux = (ρ * c) * κ * (T_ocean - T_sfc) / max(δ, eps(eltype(T_sfc))) # avoid division by zero
    Q_LW_out = ε * σ * T_sfc^4
    flux_out = Q_turb + (Q_LW_out - (1 - α) * SW_d - LW_d)
    return diffusive_flux - flux_out
end

# ------------------------------------------------------------------------------
# Ocean skin T_sfc from flux balance
# ------------------------------------------------------------------------------
# Rearranging the balance (ρ*c)*κ*(T_ocean - T_sfc)/δ = Q_turb(T_sfc) + εσT_sfc^4 - (1-α)*SW_d - LW_d
# gives the implicit equation for T_sfc:
#
#   T_sfc = T_ocean - (δ / ((ρ*c)*κ)) * ( Q_turb(T_sfc) + εσT_sfc^4 - (1-α)*SW_d - LW_d )
#
# Define flux_out(T_sfc) = Q_turb(T_sfc) + εσT_sfc^4 - (1-α)*SW_d - LW_d (outgoing flux, W/m²).
# Then T_sfc = T_ocean - (δ / ((ρ*c)*κ)) * flux_out(T_sfc). T_sfc has no closed form because
# Q_turb(T_sfc) and T_sfc^4 are nonlinear; we solve residual(T_sfc) = 0 with Brent's method.

"""
    ocean_update_T_sfc(param_set, thermo_params, inputs, g_h, SW_d, LW_d, α, ε, σ, T_ocean, δ, κ, ρ, c)

Return ocean skin temperature T_sfc [K] satisfying the flux balance
  (ρ*c)*κ*(T_ocean - T_sfc)/δ = Q_turb(T_sfc) + εσT_sfc^4 - (1-α)*SW_d - LW_d,
i.e. the implicit equation
  T_sfc = T_ocean - (δ/((ρ*c)*κ)) * ( Q_turb(T_sfc) + εσT_sfc^4 - (1-α)*SW_d - LW_d ).

Uses the current MOST state via `g_h` and `inputs` to compute Q_turb(T_sfc) = shf + lhf.
Solves for T_sfc with RootSolvers (Brent). All arguments are scalars for one column.
"""
function ocean_update_T_sfc(
    param_set,
    thermo_params,
    inputs,
    g_h,
    SW_d,
    LW_d,
    α,
    ε,
    σ,
    T_ocean,
    δ,
    κ,
    ρ,
    c,
)
    function Q_turb(T_sfc)
        q_vap_sfc_approx =
            TD.q_vap_saturation_generic(thermo_params, T_sfc, inputs.ρ_int, TD.Liquid())
        ρ_sfc = SF.surface_density(
            param_set,
            inputs.T_int,
            inputs.ρ_int,
            T_sfc,
            inputs.Δz,
            inputs.q_tot_int,
            inputs.q_liq_int,
            inputs.q_ice_int,
            q_vap_sfc_approx,
        )
        q_vap_sfc =
            TD.q_vap_saturation_generic(thermo_params, T_sfc, ρ_sfc, TD.Liquid())
        E = SF.evaporation(
            param_set,
            inputs,
            g_h,
            inputs.q_tot_int,
            q_vap_sfc,
            ρ_sfc,
            SF.MoistModel(),
        )
        lhf = SF.latent_heat_flux(param_set, inputs, E, SF.MoistModel())
        shf = SF.sensible_heat_flux(
            param_set,
            inputs,
            g_h,
            inputs.T_int,
            T_sfc,
            ρ_sfc,
            E,
        )
        return shf + lhf
    end

    residual(T_sfc) = ocean_flux_balance_residual(
        T_sfc,
        Q_turb(T_sfc),
        SW_d,
        LW_d,
        α,
        ε,
        σ,
        T_ocean,
        δ,
        κ,
        ρ,
        c,
    )
    FT = eltype(param_set)
    T_lo = T_ocean - FT(5)
    T_hi = T_ocean + FT(5)
    sol = find_zero(residual, BrentMethod(T_lo, T_hi), CompactSolution())
    return sol.root
end

"""
    ocean_update_T_sfc_callback(sw_down, lw_down, albedo, emissivity, stefan_boltzmann,
        T_ocean_bulk, skin_thickness, thermal_diffusivity, reference_density, heat_capacity)

Return an `update_T_sfc` callback for SurfaceFluxes.surface_fluxes.

The returned function takes the current MOST state (ζ, param_set, thermo_params, inputs,
scheme, u_star, z0m, z0h), computes the heat exchange coefficient from that state, and
solves for the ocean skin temperature satisfying the flux balance. Use when coupling to
SurfaceFluxes (e.g. Oceananigans with iterative ocean skin).
"""
function ocean_update_T_sfc_callback(
    sw_down,
    lw_down,
    albedo,
    emissivity,
    stefan_boltzmann,
    T_ocean_bulk,
    skin_thickness,
    thermal_diffusivity,
    reference_density,
    heat_capacity,
)
    function update_T_sfc(ζ, param_set, thermo_params, inputs, scheme, u_star, z0m, z0h)
        FT = eltype(param_set)
        effective_height = inputs.Δz - inputs.d
        Ch = SF.heat_exchange_coefficient(param_set, ζ, z0m, z0h, effective_height, scheme)
        Cd = SF.drag_coefficient(param_set, ζ, z0m, effective_height, scheme)
        wind_speed_scale = u_star / sqrt(max(Cd, eps(FT)))
        heat_exchange_coef = Ch * wind_speed_scale
        return ocean_update_T_sfc(
            param_set,
            thermo_params,
            inputs,
            heat_exchange_coef,
            sw_down,
            lw_down,
            albedo,
            emissivity,
            stefan_boltzmann,
            T_ocean_bulk,
            skin_thickness,
            thermal_diffusivity,
            reference_density,
            heat_capacity,
        )
    end
    return update_T_sfc
end
