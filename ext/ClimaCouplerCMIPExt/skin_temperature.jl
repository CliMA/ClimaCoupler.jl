"""
    Skin Temperature Calculation Utilities

This module provides functions to compute skin temperature using flux balance equations.
The flux balance is solved by computing

            κ
Jᵃ(Tₛⁿ) + --- (Tₛⁿ⁺¹ - Tᵢ) = 0
            δ

where Jᵃ is the external flux impinging on the surface from above and
Jᵢ = - κ (Tₛ - Tᵢ) / δ is the "internal flux" coming up from below.
"""

import SurfaceFluxes as SF
import Thermodynamics as TD

"""
    update_T_sfc(κ, δ, T_i, σ, ϵ, SW_d, LW_d, α_albedo)

Create a callback function for updating surface temperature that can be passed to
SurfaceFluxes.jl's `surface_fluxes` function.

The returned callback matches the SurfaceFluxes.jl callback signature:
`(ζ, param_set, thermo_params, inputs, scheme, u_star, z0m, z0s) -> T_sfc_new`

The callback computes the updated surface temperature based on the semi-implicit flux balance equation:

    Tₛⁿ⁺¹ = (Tᵢ - δ / κ * (Jᵃ - 4 σ ϵ Tₛⁿ⁴)) / (1 + 4 δ σ ϵ Tₛⁿ³ / κ)

This formulation linearizes the outgoing longwave radiation term σϵT⁴ around the current
surface temperature Tₛⁿ, providing stability for the iterative solution of Jᵃ + κ/δ*(T-Tᵢ) = 0.

Inside the callback, sensible heat flux (F_sh) and latent heat flux (F_lh) are computed directly
using SurfaceFluxes.jl helper functions, and the net upward flux J_a is computed as:

    J_a = σ * ϵ * T_sfc_n^4 - (1 - α_albedo) * SW_d - ϵ * LW_d + F_sh + F_lh

where F_sh and F_lh follow the SurfaceFluxes.jl sign convention (positive = upward).

# Arguments
- `κ`: Thermal conductivity [W m⁻¹ K⁻¹]
- `δ`: Thickness/depth of the layer [m]
- `T_i`: Internal temperature (temperature at the base of the layer) [K]
- `σ`: Stefan-Boltzmann constant [W m⁻² K⁻⁴]
- `ϵ`: Surface emissivity [-]
- `SW_d`: Downward shortwave radiation [W m⁻²]
- `LW_d`: Downward longwave radiation [W m⁻²]
- `α_albedo`: Surface albedo [-]

# Returns
- `callback`: A function matching SurfaceFluxes.jl callback signature that can be passed to
              `SF.surface_fluxes(..., update_T_sfc_callback, ...)`.
"""
function update_T_sfc(
    κ,
    δ,
    T_i,
    σ,
    ϵ,
    SW_d,
    LW_d,
    α_albedo,
)
    return function(ζ, param_set, thermo_params_callback, inputs, scheme, u_star, z0m, z0s)
        T_sfc_n = inputs.T_sfc_guess
        
        # Compute surface humidity from T_sfc_n
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
        q_vap_sfc = TD.q_vap_saturation(thermo_params_callback, T_sfc_n, ρ_sfc, inputs.q_liq_int, inputs.q_ice_int)
        
        # Compute sensible and latent heat fluxes using SurfaceFluxes.jl helpers
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
        
        # Compute net upward flux Jᵃ at the surface (positive = heat leaving surface)
        # Jᵃ = (LW emission) - (absorbed SW) - (absorbed LW) + (sensible up) + (latent up)
        # where F_sh, F_lh are positive upward (SurfaceFluxes.jl convention)
        # Old (incorrect sign convention for radiation terms):
        # J_a = (1 - α_albedo) * SW_d + LW_d + F_sh + F_lh - σ * ϵ * T_sfc_n^4
        J_a = σ * ϵ * T_sfc_n^4 - (1 - α_albedo) * SW_d - ϵ * LW_d + F_sh + F_lh
        
        # Apply semi-implicit flux balance equation.
        # Linearize the LW emission σϵT⁴ around T_n:
        #   σϵT⁴ ≈ σϵTₙ⁴ + 4σϵTₙ³(T - Tₙ) = -3σϵTₙ⁴ + 4σϵTₙ³T
        # Substituting into Jᵃ + κ/δ*(T - Tᵢ) = 0 and solving for T:
        numerator = T_i - (δ / κ) * (J_a - 4 * σ * ϵ * T_sfc_n^4)
        denominator = 1 + 4 * δ * σ * ϵ * T_sfc_n^3 / κ
        
        return numerator / denominator
    end
end
