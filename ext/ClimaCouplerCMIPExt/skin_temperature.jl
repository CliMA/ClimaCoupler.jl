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

"""
    surface_turbulent_energy_flux(ζ, param_set, thermo_params, inputs, scheme,
                                  u_star, z0m, z0s, T_sfc)

Total upward turbulent energy flux `F_sh + F_lh` at a trial surface temperature
`T_sfc`, with surface density and saturation specific humidity evaluated
consistently at `T_sfc`.
"""
function surface_turbulent_energy_flux(
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
    q_vap_sfc =
        TD.q_vap_saturation(thermo_params, T_sfc, ρ_sfc, inputs.q_liq_int, inputs.q_ice_int)

    F_sh = SF.sensible_heat_flux(
        param_set,
        ζ,
        u_star,
        inputs,
        z0m,
        z0s,
        T_sfc,
        q_vap_sfc,
        inputs.ρ_int,
        scheme,
    )
    F_lh = SF.latent_heat_flux(
        param_set,
        ζ,
        u_star,
        inputs,
        z0m,
        z0s,
        q_vap_sfc,
        inputs.ρ_int,
        scheme,
    )
    return F_sh + F_lh
end

"""
    update_T_sfc(R, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt)

Create a callback for `SurfaceFluxes.jl` that updates surface temperature using
a fully-linearized (Newton) step of the skin balance Jᵃ(Tₛ) + (Tₛ - Tᵢ)/R = 0:

    Tₛⁿ⁺¹ = (Tᵢ - R · (Jᵃ - Λ Tₛⁿ)) / (1 + R Λ),    Λ = 4σϵTₛⁿ³ + ∂F_turb/∂Tₛ

where Jᵃ = σϵTₛⁿ⁴ - (1-α)SW↓ - ϵLW↓ + F_sh + F_lh  (positive upward) and
F_turb = F_sh + F_lh. 

Single iteration update is safeguarded by
±ΔT_iter_max, and the result is capped at the melting
temperature T_melt to prevent the surface temperature from exceeding the
melting point under heating fluxes.

# Arguments
- `R`: Conductive resistance of the column [m² K W⁻¹] (snow and ice in series)
- `T_i`: Internal (ice-ocean interface) temperature [K]
- `σ`: Stefan-Boltzmann constant [W m⁻² K⁻⁴]
- `ϵ`: Surface emissivity [-]
- `SW_d`: Downward shortwave radiation [W m⁻²]
- `LW_d`: Downward longwave radiation [W m⁻²]
- `α_albedo`: Surface albedo [-]
- `T_melt`: Melting temperature [K] (typically 273.15 K for freshwater ice)
"""
function update_T_sfc(R, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt)
    return function (ζ, param_set, thermo_params_callback, inputs, scheme, u_star, z0m, z0s)
        T_sfc_n = inputs.T_sfc_guess
        FT = typeof(T_sfc_n)

        F_turb = surface_turbulent_energy_flux(
            ζ,
            param_set,
            thermo_params_callback,
            inputs,
            scheme,
            u_star,
            z0m,
            z0s,
            T_sfc_n,
        )

        # Turbulent sensitivity ∂F_turb/∂Tₛ by one-sided finite difference,
        # clamped to ≥ 0 so it can only damp the update (the physical
        # sensitivity is positive: warmer surface ⇒ larger upward fluxes).
        δT = FT(1)
        F_turb_δ = surface_turbulent_energy_flux(
            ζ,
            param_set,
            thermo_params_callback,
            inputs,
            scheme,
            u_star,
            z0m,
            z0s,
            T_sfc_n + δT,
        )
        dF_turb_dT = max((F_turb_δ - F_turb) / δT, zero(FT))

        # Net upward flux: Jᵃ = σϵT⁴ - (1-α)SW↓ - ϵLW↓ + F_sh + F_lh
        J_a = σ * ϵ * T_sfc_n^4 - (1 - α_albedo) * SW_d - ϵ * LW_d + F_turb

        # Newton step on Jᵃ(T) + (T - Tᵢ)/R = 0 with
        # Jᵃ(T) ≈ Jᵃ(Tₙ) + Λ (T - Tₙ)
        Λ = 4 * σ * ϵ * T_sfc_n^3 + dF_turb_dT
        T_sfc_new = (T_i - R * (J_a - Λ * T_sfc_n)) / (1 + R * Λ)

        T_sfc_new = ifelse(isnan(T_sfc_new), T_sfc_n, T_sfc_new)
        ΔT = T_sfc_new - T_sfc_n
        ΔT_iter_max = FT(5)
        T_sfc_new = T_sfc_n + min(abs(ΔT), ΔT_iter_max) * sign(ΔT)
        T_sfc_new = min(T_sfc_new, T_melt)

        return T_sfc_new
    end
end
