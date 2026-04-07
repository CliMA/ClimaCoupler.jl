"""
    Skin Temperature Calculation

Solve the steady-state flux balance for surface temperature ``T_s``.

**Conductive (Fourier) flux through ice** — [`update_T_sfc`](@ref):

```math
J^a(T_s) + \\frac{\\kappa}{\\delta} (T_s - T_i) = 0
```

with ``\\kappa`` thermal conductivity [``\\mathrm{W\\,m^{-1}\\,K^{-1}}``] and ``\\delta`` ice thickness [m].

**Diffusive boundary-layer flux** (same closure as
`ClimaOcean.OceanSeaIceModels.InterfaceComputations.SkinTemperature` with
`DiffusiveFlux`) — [`update_T_sfc_diffusive`](@ref):

```math
\\frac{\\kappa}{\\delta} (T_s - T_i) + \\frac{1}{\\lambda} \\left( Q_a + \\Omega_c (T_a - T_s) \\right) = 0 ,
\\quad \\lambda = \\frac{1}{\\rho c}
```

with ``\\kappa`` kinematic diffusivity [``\\mathrm{m^2\\,s^{-1}}``], ``\\delta`` a boundary-layer thickness [m],
``Q_a`` the net upward non-sensible surface heat flux [``\\mathrm{W\\,m^{-2}}``],
``\\Omega_c = Q_c / (T_a - T_s)`` the sensible exchange coefficient in ``(\\mathrm{W\\,m^{-2}\\,K^{-1}})\\cdot\\lambda`` form,
and reference ``\\rho``, ``c`` for converting heat flux to the diffusive balance.

See `ClimaOcean` `interface_states.jl` (`flux_balance_temperature` for `SkinTemperature{<:DiffusiveFlux}`).
"""

import SurfaceFluxes as SF
import Thermodynamics as TD

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

        # Net upward flux: Jᵃ = σϵT⁴ - (1-α)SW↓ - ϵLW↓ + F_sh + F_lh
        J_a = σ * ϵ * T_sfc_n^4 - (1 - α_albedo) * SW_d - ϵ * LW_d + F_sh + F_lh

        # Semi-implicit solve: linearize σϵT⁴ ≈ -3σϵTₙ⁴ + 4σϵTₙ³T
        numerator = T_i - (δ / κ) * (J_a - 4 * σ * ϵ * T_sfc_n^4)
        denominator = 1 + 4 * δ * σ * ϵ * T_sfc_n^3 / κ
        T_sfc_new = numerator / denominator

        # Cap surface temperature at melting temperature 
        T_sfc_new = min(T_sfc_new, T_melt)

        return T_sfc_new
    end
end

"""
    update_T_sfc_diffusive(κ, δ_bl, ρ_ref, c_ref, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt)

Like [`update_T_sfc`](@ref), but use the **diffusive closure** from ClimaOcean
`SkinTemperature(DiffusiveFlux(δ, κ))`: fixed-point surface temperature

```math
T_s = \\frac{T_i \\kappa - (J^T + \\Omega_c T_a) \\, \\delta}{\\kappa - \\Omega_c \\, \\delta}
```

with ``J^T = Q_a \\, /(\\rho_{\\mathrm{ref}} c_{\\mathrm{ref}})``,
``Q_a = \\sigma\\epsilon T_s^n{}^4 - (1-\\alpha)SW^- - \\epsilon LW^- + F_{lh}`` (non-sensible, upward positive),
``\\Omega_c = F_{sh} / (T_a - T_s^n) \\, /(\\rho_{\\mathrm{ref}} c_{\\mathrm{ref}})`` when ``T_a \\neq T_s^n``, else zero,
and ``T_a = T_{\\mathrm{int}} + g\\, \\Delta z / c_{p,m}`` (same lapse correction as ClimaOcean `surface_atmosphere_temperature`).

# Arguments
- `κ`: Kinematic diffusivity [``\\mathrm{m^2\\,s^{-1}}``]
- `δ_bl`: Boundary-layer thickness [m]
- `ρ_ref`: Reference density [``\\mathrm{kg\\,m^{-3}}``] (typical seawater ``~1025``)
- `c_ref`: Reference heat capacity [``\\mathrm{J\\,kg^{-1}\\,K^{-1}}``] (typical seawater ``~4000``)
- `T_i`: Interior (ice–ocean interface) temperature [K]
- Remaining arguments: as in [`update_T_sfc`](@ref)
"""
function update_T_sfc_diffusive(κ, δ_bl, ρ_ref, c_ref, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt)
    return function (ζ, param_set, thermo_params_callback, inputs, scheme, u_star, z0m, z0s)
        T_sfc_n = inputs.T_sfc_guess

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

        g = TD.Parameters.grav(thermo_params_callback)
        c_a = TD.cp_m(thermo_params_callback, inputs.q_tot_int)
        T_a = inputs.T_int + g * inputs.Δz / c_a

        Qa =
            σ * ϵ * T_sfc_n^4 -
            (1 - α_albedo) * SW_d - ϵ * LW_d +
            F_lh
        λ = inv(ρ_ref * c_ref)
        Jᴛ = Qa * λ
        ΔT_atm = T_a - T_sfc_n
        Ωc = ifelse(ΔT_atm == 0, zero(ΔT_atm), F_sh / ΔT_atm * λ)

        denom = κ - Ωc * δ_bl
        T_sfc_new = (T_i * κ - (Jᴛ + Ωc * T_a) * δ_bl) / denom
        T_sfc_new = ifelse(isnan(T_sfc_new), T_sfc_n, T_sfc_new)
        T_sfc_new = min(T_sfc_new, T_melt)

        return T_sfc_new
    end
end
