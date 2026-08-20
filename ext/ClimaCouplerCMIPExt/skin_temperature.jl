"""
    Skin Temperature Calculation

Solve the steady-state flux balance for surface temperature Tₛ:

Jᵃ(Tₛ) + (Tₛ - Tᵢ)/R = 0

where Jᵃ is the net upward surface flux and (Tₛ - Tᵢ)/R is the conductive flux
through the column. `R` is the conductive resistance [m² K W⁻¹]: for bare ice
`R = δ/κ`, and with a snow layer the snow and ice resistances add in series,
`R = h_snow/k_snow + h_ice/k_ice`.
"""

import SurfaceFluxes as SF
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

"""
    surface_turbulent_energy_flux(ζ, param_set, thermo_params, inputs, scheme,
                                  u_star, z0m, z0s, T_sfc)

Total upward turbulent energy flux `F_sh + F_lh` at a trial surface temperature
`T_sfc`, with surface density and saturation specific humidity evaluated
consistently at `T_sfc`. Prefer [`turbulent_energy_flux_and_sensitivity`](@ref)
when the skin-temperature Newton step also needs `∂F_turb/∂Tₛ`.
"""
@inline function surface_turbulent_energy_flux(
    ζ,
    param_set,
    thermo_params,
    inputs,
    scheme,
    u_star,
    z0m,
    z0s,
    T_sfc,
)
    F_turb, _ = turbulent_energy_flux_and_sensitivity(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        u_star,
        z0m,
        z0s,
        T_sfc,
    )
    return F_turb
end

"""
    turbulent_energy_flux_and_sensitivity(ζ, param_set, thermo_params, inputs,
                                          scheme, u_star, z0m, z0s, T_sfc)

Return `(F_turb, ∂F_turb/∂Tₛ)` at `T_sfc` from a single conductance evaluation.

`F_turb = F_sh + F_lh` uses the SurfaceFluxes bulk formulas. The sensitivity
holds `g_h` and `ρ_sfc` fixed (standard Newton linearization of the bulk
scheme):

```
∂F_turb/∂Tₛ = ρ_sfc g_h (cₚ,d + (Lᵥ₀ + VSE_sfc) ∂q*/∂T) + cₚ,v E
```

clamped to `≥ 0` so it can only damp the skin-temperature update.
"""
@inline function turbulent_energy_flux_and_sensitivity(
    ζ,
    param_set,
    thermo_params,
    inputs,
    scheme,
    u_star,
    z0m,
    z0s,
    T_sfc,
)
    FT = typeof(T_sfc)
    g_h = SF.heat_conductance(param_set, ζ, u_star, inputs, z0m, z0s, scheme)
    ρ_sfc = SF.surface_density(
        param_set,
        inputs.T_int,
        inputs.ρ_int,
        T_sfc,
        inputs.Δz,
        inputs.q_tot_int,
        inputs.q_liq_int,
        inputs.q_ice_int,
    )
    q_vap_sfc = TD.q_vap_saturation(
        thermo_params,
        T_sfc,
        ρ_sfc,
        inputs.q_liq_int,
        inputs.q_ice_int,
    )
    q_vap_int = inputs.q_tot_int - inputs.q_liq_int - inputs.q_ice_int
    E = SF.evaporation(
        param_set,
        inputs,
        g_h,
        q_vap_int,
        q_vap_sfc,
        ρ_sfc,
        inputs.moisture_model,
    )
    F_sh = SF.sensible_heat_flux(
        param_set,
        inputs,
        g_h,
        inputs.T_int,
        T_sfc,
        ρ_sfc,
        E,
    )
    F_lh = SF.latent_heat_flux(param_set, inputs, E, inputs.moisture_model)
    F_turb = F_sh + F_lh

    # Analytic ∂F_turb/∂Tₛ from the bulk formulas (g_h, ρ_sfc held fixed).
    cp_d = TDP.cp_d(thermo_params)
    cp_v = TDP.cp_v(thermo_params)
    LH_v0 = TDP.LH_v0(thermo_params)
    dq_dT = TD.∂q_vap_sat_∂T(
        thermo_params,
        T_sfc,
        ρ_sfc,
        inputs.q_liq_int,
        inputs.q_ice_int,
    )
    VSE_sfc = TD.vapor_static_energy(thermo_params, T_sfc, inputs.Φ_sfc)
    dF_dT = ρ_sfc * g_h * (cp_d + (LH_v0 + VSE_sfc) * dq_dT) + cp_v * E
    return F_turb, max(dF_dT, zero(FT))
end

"""
    UpdateTSfc{FT}

Skin-temperature Newton step for SurfaceFluxes (functor, not a capturing
closure): `Jᵃ(Tₛ) + (Tₛ - Tᵢ)/R = 0` with
`Λ = 4σϵTₛ³ + ∂F_turb/∂Tₛ` from
[`turbulent_energy_flux_and_sensitivity`](@ref). Clamp ±5 K, cap at `T_melt`.
"""
struct UpdateTSfc{FT}
    R::FT
    T_i::FT
    σ::FT
    ϵ::FT
    SW_d::FT
    LW_d::FT
    α_albedo::FT
    T_melt::FT
end

# Broadcast-friendly constructor (`update_T_sfc.(…)` on the boundary path).
update_T_sfc(R, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt) =
    UpdateTSfc(promote(R, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt)...)

@inline function (cb::UpdateTSfc)(
    ζ,
    param_set,
    thermo_params,
    inputs,
    scheme,
    u_star,
    z0m,
    z0s,
)
    T_sfc_n = inputs.T_sfc_guess
    FT = typeof(T_sfc_n)

    F_turb, dF_turb_dT = turbulent_energy_flux_and_sensitivity(
        ζ,
        param_set,
        thermo_params,
        inputs,
        scheme,
        u_star,
        z0m,
        z0s,
        T_sfc_n,
    )

    # Net upward flux: Jᵃ = σϵT⁴ - (1-α)SW↓ - ϵLW↓ + F_sh + F_lh
    J_a =
        cb.σ * cb.ϵ * T_sfc_n^4 - (1 - cb.α_albedo) * cb.SW_d - cb.ϵ * cb.LW_d +
        F_turb

    # Newton step on Jᵃ(T) + (T - Tᵢ)/R = 0 with Jᵃ(T) ≈ Jᵃ(Tₙ) + Λ (T - Tₙ)
    Λ = 4 * cb.σ * cb.ϵ * T_sfc_n^3 + dF_turb_dT
    T_sfc_new = (cb.T_i - cb.R * (J_a - Λ * T_sfc_n)) / (1 + cb.R * Λ)

    T_sfc_new = ifelse(isnan(T_sfc_new), T_sfc_n, T_sfc_new)
    ΔT = T_sfc_new - T_sfc_n
    ΔT_iter_max = FT(5)
    T_sfc_new = T_sfc_n + min(abs(ΔT), ΔT_iter_max) * sign(ΔT)
    T_sfc_new = min(T_sfc_new, cb.T_melt)

    return T_sfc_new
end

# Fixed MOST iteration count (GPU-uniform; no early exit).
@inline ice_surface_flux_solver_opts(::Type{FT}) where {FT} =
    SF.SolverOptions{FT}(forced_fixed_iters = true)
