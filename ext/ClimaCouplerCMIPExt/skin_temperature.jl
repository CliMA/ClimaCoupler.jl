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

"""
    update_T_sfc(κ, δ, T_i, σ, ϵ, ρ, c, J_a_fn)

Create a callback function for updating surface temperature that can be passed to
SurfaceFluxes.jl or used in iterative flux calculations.

The returned function takes the current surface temperature `T_sfc_n` and returns
the updated surface temperature `T_sfc_new` based on the semi-implicit flux balance equation:

    Tₛⁿ⁺¹ = (Tᵢ - δ / κ * (Jᵃ - 4 α Tₛⁿ⁴)) / (1 + 4 δ σ ϵ Tₛⁿ³ / ρ c κ)

where α = σ ϵ / (ρ c). This formulation linearizes the outgoing longwave radiation term
around the current surface temperature, providing better stability for iterative solutions.

# Arguments
- `κ`: Thermal conductivity [W m⁻¹ K⁻¹]
- `δ`: Thickness/depth of the layer [m]
- `T_i`: Internal temperature (temperature at the base of the layer) [K]
- `σ`: Stefan-Boltzmann constant [W m⁻² K⁻⁴]
- `ϵ`: Surface emissivity [-]
- `ρ`: Density [kg m⁻³]
- `c`: Specific heat capacity [J kg⁻¹ K⁻¹]
- `J_a_fn`: Function that computes the external flux `J_a` given `T_sfc_n`.
            Signature: `J_a_fn(T_sfc_n) -> J_a [W m⁻²]`
            This function should compute the net flux impinging on the surface,
            typically including: `(1 - α) * SW_d + LW_d + F_sh + F_lh - σ * ϵ * T_sfc_n^4`

# Returns
- `callback`: A function `callback(T_sfc_n) -> T_sfc_new` that can be used as a callback
              in SurfaceFluxes.jl or iterative flux calculations.
"""
function update_T_sfc(
    κ,
    δ,
    T_i,
    σ,
    ϵ,
    ρ,
    c,
    J_a_fn,
)
    α = σ * ϵ / (ρ * c)
    
    return function(T_sfc_n)
        T_sfc_n³ = T_sfc_n^3
        T_sfc_n⁴ = T_sfc_n^4

        numerator = T_i - (δ / κ) * (J_a_fn(T_sfc_n) - 4 * α * T_sfc_n⁴)
        denominator = 1 + 4 * δ * σ * ϵ * T_sfc_n³ / (ρ * c * κ)

        return numerator / denominator
    end
end
