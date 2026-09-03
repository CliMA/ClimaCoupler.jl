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
    update_T_sfc(R, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt, maximum_temperature_step)

Create a callback for `SurfaceFluxes.jl` that updates surface temperature from a
semi-implicit linearization of the surface flux balance

    Qᵃ(Tₛ) + Ωᵀ (Tᵃ - Tₛ) + (Tₛ - Tᵢ)/R = 0

where `Qᵃ = σϵTₛ⁴ - (1-α)SW↓ - ϵLW↓ + F_lh` collects every upward flux except the sensible heat, and
`Ωᵀ = F_sh/(Tᵃ - Tₛ)` is the linearized sensible heat coefficient. Linearizing the emitted longwave as
`σϵTₛ⁴ ≈ σϵTₙ⁴ + β(Tₛ - Tₙ)` with `β = 4σϵTₙ³` gives the Newton-like update

    Tₛⁿ⁺¹ = (Tᵢ + βR Tₙ - Ωᵀ R Tᵃ - Qᵃ R) / (1 + βR - Ωᵀ R)

Both the radiative and the sensible response to Tₛ appear in the denominator. Carrying the sensible term
implicitly is what makes the iteration a contraction: `Ωᵀ ≈ -ρ cₚ Cₕ |U|` is an order of magnitude larger
than `β` under Arctic conditions, so a scheme that holds it fixed amplifies the error by
`|Ωᵀ| R / (1 + βR)` per iteration and oscillates between the melting point and arbitrarily cold values.

The step is capped at `maximum_temperature_step` per iteration, and the result at the melting temperature
`T_melt`.

# Arguments
- `R`: Conductive resistance of the column [m² K W⁻¹] (snow and ice in series)
- `T_i`: Internal (ice-ocean interface) temperature [K]
- `σ`: Stefan-Boltzmann constant [W m⁻² K⁻⁴]
- `ϵ`: Surface emissivity [-]
- `SW_d`: Downward shortwave radiation [W m⁻²]
- `LW_d`: Downward longwave radiation [W m⁻²]
- `α_albedo`: Surface albedo [-]
- `T_melt`: Melting temperature [K] (typically 273.15 K for freshwater ice)
- `maximum_temperature_step`: Largest surface temperature change allowed in a single iteration [K]
"""
function update_T_sfc(R, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt, maximum_temperature_step)
    return function (ζ, param_set, thermo_params_callback, inputs, scheme, u_star, z0m, z0s)
        Tₛ⁻ = inputs.T_sfc_guess

        # Surface density and saturation specific humidity at Tₛ⁻
        ρₛ = SF.surface_density(
            param_set,
            inputs.T_int,
            inputs.ρ_int,
            Tₛ⁻,
            inputs.Δz,
            inputs.q_tot_int,
            inputs.q_liq_int,
            inputs.q_ice_int,
        )
        qₛ = TD.q_vap_saturation(
            thermo_params_callback,
            Tₛ⁻,
            ρₛ,
            inputs.q_liq_int,
            inputs.q_ice_int,
        )

        𝒬ᵀ = SF.sensible_heat_flux(
            param_set,
            ζ,
            u_star,
            inputs,
            z0m,
            z0s,
            Tₛ⁻,
            qₛ,
            ρₛ,
            scheme,
        )
        𝒬ᵛ = SF.latent_heat_flux(param_set, ζ, u_star, inputs, z0m, z0s, qₛ, ρₛ, scheme)

        # Everything upward except the sensible heat, which is carried implicitly below
        Qᵃ = σ * ϵ * Tₛ⁻^4 - (1 - α_albedo) * SW_d - ϵ * LW_d + 𝒬ᵛ

        # Secant slope of the sensible heat flux. Physically Ωᵀ ≤ 0; the geopotential
        # offset between T_int and the dry static energy makes the quotient change sign
        # for |Δ| below ~0.1 K, where 𝒬ᵀ itself is negligible.
        Tᵃ = inputs.T_int
        Δ = Tᵃ - Tₛ⁻
        Ωᵀ = min(𝒬ᵀ / Δ, zero(Δ))
        Ωᵀ = ifelse(isfinite(Ωᵀ), Ωᵀ, zero(Δ))

        # Newton linearization of the emitted longwave: σϵTₛ⁴ ≈ σϵTₛ⁻⁴ + β (Tₛ - Tₛ⁻)
        β = 4 * σ * ϵ * Tₛ⁻^3

        D = 1 + β * R - Ωᵀ * R
        Tₛ = (T_i + β * R * Tₛ⁻ - Ωᵀ * R * Tᵃ - Qᵃ * R) / D
        Tₛ = ifelse((D == zero(D)) | isnan(Tₛ), Tₛ⁻, Tₛ)

        # Cap the step for iteration stability, then cap at the melting temperature
        Tₛ = Tₛ⁻ + clamp(Tₛ - Tₛ⁻, -maximum_temperature_step, maximum_temperature_step)
        return min(Tₛ, T_melt)
    end
end
