"""
    Skin Temperature Calculation

Solve the steady-state flux balance for surface temperature Tₛ:

            κ
Jᵃ(Tₛ) + --- (Tₛ - Tᵢ) = 0
            δ

where Jᵃ is the net upward surface flux and κ/δ·(Tₛ - Tᵢ) is the conductive flux.
"""

import SurfaceFluxes as SF
import Thermodynamics as TD

"""
    ice_skin_J_a(T_sfc, SW_d, LW_d, α_albedo, F_sh, F_lh, σ, ϵ)

Net upward atmospheric flux at the ice surface, ``J^a`` [W m⁻²] (**positive upward**), using the
same definition as `update_T_sfc`:

``J^a = \\sigma\\epsilon T_\\mathrm{sfc}^4 - (1-\\alpha)\\mathrm{SW}_\\downarrow - \\epsilon\\,\\mathrm{LW}_\\downarrow + F_{sh} + F_{lh}``.

Use with `T_sfc`, `F_sh`, and `F_lh` from the **same** `get_surface_fluxes` / SurfaceFluxes evaluation
you want to compare (e.g. `T_sfc_new` after the flux call). For diagnostics on the coupler boundary
space, see `csf.sea_ice_skin_J_a` filled in `compute_surface_fluxes!` for `ClimaSeaIceSimulation`.
"""
@inline function ice_skin_J_a(T_sfc, SW_d, LW_d, α_albedo, F_sh, F_lh, σ, ϵ)
    return σ * ϵ * T_sfc^4 - (1 - α_albedo) * SW_d - ϵ * LW_d + F_sh + F_lh
end

"""
    update_T_sfc(κ, δ, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt)

Create a callback for `SurfaceFluxes.jl` that updates surface temperature using
a semi-implicit linearization of the LW emission term:

    Tₛⁿ⁺¹ = (Tᵢ - δ/κ · (Jᵃ - 4σϵTₛⁿ⁴)) / (1 + 4δσϵTₛⁿ³/κ)

where Jᵃ = σϵTₛⁿ⁴ - (1-α)SW↓ - ϵLW↓ + F_sh + F_lh  (positive upward).

The result is capped at the melting temperature T_melt to prevent the surface
temperature from exceeding the melting point under heating fluxes.

# Arguments
- `κ`: Thermal conductivity [W m⁻¹ K⁻¹]
- `δ`: Ice thickness [m]
- `T_i`: Internal (ice-ocean interface) temperature [K]
- `σ`: Stefan-Boltzmann constant [W m⁻² K⁻⁴]
- `ϵ`: Surface emissivity [-]
- `SW_d`: Downward shortwave radiation [W m⁻²]
- `LW_d`: Downward longwave radiation [W m⁻²]
- `α_albedo`: Surface albedo [-]
- `T_melt`: Melting temperature [K] (typically 273.15 K for freshwater ice)
"""
function update_T_sfc(κ, δ, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt)
    return function (ζ, param_set, thermo_params_callback, inputs, scheme, u_star, z0m, z0s)
        T_sfc_n = inputs.T_sfc_guess

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

        J_a = ice_skin_J_a(T_sfc_n, SW_d, LW_d, α_albedo, F_sh, F_lh, σ, ϵ)

        # Semi-implicit solve: linearize σϵT⁴ ≈ -3σϵTₙ⁴ + 4σϵTₙ³T
        numerator = T_i - (δ / κ) * (J_a - 4 * σ * ϵ * T_sfc_n^4)
        denominator = 1 + 4 * δ * σ * ϵ * T_sfc_n^3 / κ
        T_sfc_new = numerator / denominator

        # Cap surface temperature at melting temperature 
        T_sfc_new = min(T_sfc_new, T_melt)

        return T_sfc_new
    end
end
