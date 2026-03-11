"""
    Skin Temperature Calculation

Solve the steady-state flux balance for surface temperature Tₛ:

            κ
Jᵃ(Tₛ) + --- (Tₛ - Tᵢ) = 0
            δ

where Jᵃ is the net upward surface flux and κ/δ·(Tₛ - Tᵢ) is the conductive flux.

For sea ice, the penetrating part of shortwave is removed from the surface balance
following OIFES (Komori et al., J. Earth Simulator Vol. 4, 2005): only the
non-penetrating fraction (1 - I0*exp(-1.5*δ)) of absorbed SW heats the surface.
"""

import SurfaceFluxes as SF
import Thermodynamics as TD

"""
    update_T_sfc(κ, δ, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt, hc;
                 I0 = 0.17, max_ΔT = 5)

Create a callback for `SurfaceFluxes.jl` that updates surface temperature using
the same semi-implicit linearization of the LW emission term used in
`ClimaOcean.OceanSeaIceModels.InterfaceComputations.SkinTemperature`:

    Tₛⁿ⁺¹ = (Tᵢ - δ/κ · (Jᵃ - 4σϵTₛⁿ⁴)) / (1 + 4δσϵTₛⁿ³/κ)

where Jᵃ uses only the non-penetrating shortwave at the surface (OIFES Eq. 58, 60):
    (1-α)SW↓_surface = (1-α)SW↓ * (1 - I0*exp(-1.5*δ))
so Jᵃ = σϵTₛⁿ⁴ - (1-α)SW↓_surface - ϵLW↓ + F_sh + F_lh  (positive upward).

The result is:
- capped at the melting temperature `T_melt` (OIFES Table 2: 273.05 K for sea ice),
- limited so that the increment from the initial guess `Tₛⁿ` does not exceed `max_ΔT`,
- and, following the consolidated/thin-ice logic in
  `flux_balance_temperature(::SkinTemperature{<:ClimaSeaIce.ConductiveFlux}, ...)`,
  forced to the basal temperature `T_i` when the ice is thinner than the
  consolidation thickness `hc`.

# Arguments
- `κ`: Thermal conductivity [W m⁻¹ K⁻¹]
- `δ`: Ice thickness [m]
- `T_i`: Internal (ice-ocean interface) temperature [K]
- `σ`: Stefan-Boltzmann constant [W m⁻² K⁻⁴]
- `ϵ`: Surface emissivity [-]
- `SW_d`: Downward shortwave radiation [W m⁻²]
- `LW_d`: Downward longwave radiation [W m⁻²]
- `α_albedo`: Surface albedo [-]
- `T_melt`: Melting temperature [K] (273.05 K for sea ice, OIFES Table 2)
- `hc`: Ice consolidation thickness [m]; when `δ < hc` we use `T_i` instead of a
        diagnosed skin temperature.
- `I0`: Fraction of absorbed SW that penetrates ice (default 0.17, OIFES Eq. 58)
- `max_ΔT`: Maximum allowed change in skin temperature per call [K]
"""
function update_T_sfc(κ, δ, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt, hc;
                      I0 = 0.17, max_ΔT = 5.0)
    return function (ζ, param_set, thermo_params_callback, inputs, scheme, u_star, z0m, z0s)
        T_sfc_n = inputs.T_sfc_guess
        I0 = convert(eltype(T_sfc_n), I0)

        # Surface density and saturation specific humidity at T_sfc_n
        ρ_sfc = SF.surface_density(
            param_set,
            inputs.T_int,
            inputs.ρ_int,
            T_sfc_n,
            inputs.Δz,
            inputs.q_tot_int,
            inputs.q_liq_int,
            inputs.q_ice_int,
        )
        q_vap_sfc = TD.q_vap_saturation(
            thermo_params_callback,
            T_sfc_n,
            ρ_sfc,
            inputs.q_liq_int,
            inputs.q_ice_int,
        )

        F_sh = SF.sensible_heat_flux(
            param_set,
            ζ,
            u_star,
            inputs,
            z0m,
            z0s,
            T_sfc_n,
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

        # Net upward flux: Jᵃ = σϵT⁴ - (1-α)SW↓_surface - ϵLW↓ + F_sh + F_lh
        # OIFES Eq. 58, 60: only non-penetrating SW heats surface; (1-α)SW↓_surface = (1-α)SW↓ * (1 - I0*exp(-1.5*δ))
        sw_surface = (1 - α_albedo) * SW_d * (1 - I0 * exp(FT(-1.5) * δ))
        J_a = σ * ϵ * T_sfc_n^4 - sw_surface - ϵ * LW_d + F_sh + F_lh

        # Semi-implicit solve: linearize σϵT⁴ ≈ -3σϵTₙ⁴ + 4σϵTₙ³T
        numerator = T_i - (δ / κ) * (J_a - 4 * σ * ϵ * T_sfc_n^4)
        denominator = 1 + 4 * δ * σ * ϵ * T_sfc_n^3 / κ
        T_sfc_new = numerator / denominator

        # Limit the change in surface temperature to avoid instabilities in the
        # fixed-point iteration, mimicking SkinTemperature.max_ΔT in
        # ClimaOcean's interface_states.jl.
        ΔT = T_sfc_new - T_sfc_n
        max_ΔT_T = convert(typeof(T_sfc_new), max_ΔT)
        abs_ΔT = min(max_ΔT_T, abs(ΔT))
        T_sfc_limited = T_sfc_n + abs_ΔT * sign(ΔT)

        # Cap surface temperature at melting temperature
        T_sfc_limited = min(T_sfc_limited, T_melt)

        # Thin/unconsolidated ice: if thickness δ is less than consolidation
        # thickness hc, use the basal temperature T_i instead of a diagnosed
        # skin temperature (h/hc logic from ClimaOcean's conductive case).
        h = δ
        T_sfc_final = ifelse(h ≥ hc, T_sfc_limited, T_i)

        return T_sfc_final
    end
end
