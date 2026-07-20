# Ported from NumericalEarth v0.5.8, src/EarthSystemModels/InterfaceComputations/{sea_ice_ocean_fluxes,
# sea_ice_ocean_heat_flux_formulations,friction_velocity,InterfaceComputations}.jl, src/EarthSystemModels/
# earth_system_model.jl, and src/SeaIces/sea_ice_simulation.jl

import Oceananigans.Architectures: architecture
import Oceananigans.Grids: Flat, topology
import Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, Δzᶜᶜᶜ
import Oceananigans.Utils: KernelParameters, launch!, worksize
import ClimaSeaIce.SeaIceThermodynamics: LinearLiquidus, melting_temperature
import ClimaSeaIce.SeaIceDynamics: x_momentum_stress, y_momentum_stress

#####
##### Friction velocity
#####

"""
    MomentumBasedFrictionVelocity

A friction velocity formulation that computes the friction velocity from momentum stresses.

The friction velocity is computed as:
```math
u_* = \\sqrt{\\frac{|\\boldsymbol{\\tau}|}{\\rho_o}}
```
where τ is the magnitude of the momentum stress vector and ρᵒᶜ is the ocean reference density.
"""
struct MomentumBasedFrictionVelocity end

# ϕ² is shared with ocean_simulation.jl, which upstream keeps in a separate module.
@inline τᶜᶜᶜ(i, j, k, grid, τˣ, τʸ) =
    @inbounds sqrt(ℑxᶜᵃᵃ(i, j, k, grid, ϕ², τˣ) + ℑyᵃᶜᵃ(i, j, k, grid, ϕ², τʸ))

Base.summary(::MomentumBasedFrictionVelocity) = "MomentumBasedFrictionVelocity"

function Base.show(io::IO, ::MomentumBasedFrictionVelocity)
    print(io, "MomentumBasedFrictionVelocity (computed from momentum stresses)")
end

"""
    get_friction_velocity(u★, i, j, grid, τˣ, τʸ, ρᵒᶜ)

Return the friction velocity at grid point `(i, j)`.

For a constant friction velocity (`u★::Number`), returns the value directly.
For `MomentumBasedFrictionVelocity`, computes ``u_* = \\sqrt{|\\tau| / \\rho_o}`` from momentum stresses.
"""
@inline get_friction_velocity(u★::Number, i, j, grid, τˣ, τʸ, ρᵒᶜ) = u★
@inline get_friction_velocity(::MomentumBasedFrictionVelocity, i, j, grid, τˣ, τʸ, ρᵒᶜ) =
    sqrt(τᶜᶜᶜ(i, j, 1, grid, τˣ, τʸ) / ρᵒᶜ)

#####
##### Three-equation heat flux
#####

"""
    ThreeEquationHeatFlux{F, T, FT, U}

Three-equation formulation for sea ice-ocean heat flux.

This formulation solves a coupled system for the interface temperature and salinity:
1. Heat balance: ``\\rho c_p \\gamma_T (T - T_b) = ℰ q``
2. Salt balance: ``\\gamma_S (S - S_b) = q (S_b - S_i)``
3. Freezing point: ``T_b = T_m(S_b)``

where ``T_b`` and ``S_b`` are the interface temperature and salinity,
``\\gamma_T = \\alpha_h u_*`` and ``\\gamma_S = \\alpha_s u_*`` are turbulent exchange velocities,
``L`` is the latent heat of fusion, and ``q`` is the melt rate (computed, not input).

Fields
======

- `conductive_flux::F`: diffusive flux inside the sea ice (`ConductiveFlux`)
- `internal_temperature::T`: sea ice internal temperature field
- `heat_transfer_coefficient::FT`: turbulent heat exchange coefficient ``\\alpha_h`` (dimensionless)
- `salt_transfer_coefficient::FT`: turbulent salt exchange coefficient ``\\alpha_s`` (dimensionless)
- `friction_velocity::U`: friction velocity value or formulation (constant `Number` or `MomentumBasedFrictionVelocity`)

References
==========

- Holland, D. M., & Jenkins, A. (1999). Modeling thermodynamic ice–ocean interactions at the base of an ice
  shelf. *Journal of Physical Oceanography*, 29(8), 1787-1800.
- Shi, X., Notz, D., Liu, J., Yang, H., & Lohmann, G. (2021). Sensitivity of Northern Hemisphere climate to
  ice-ocean interface heat flux parameterizations. *Geosci. Model Dev.*, 14, 4891-4908.
"""
struct ThreeEquationHeatFlux{F, T, FT, U}
    conductive_flux::F
    internal_temperature::T
    heat_transfer_coefficient::FT
    salt_transfer_coefficient::FT
    friction_velocity::U
end

Adapt.adapt_structure(to, f::ThreeEquationHeatFlux) = ThreeEquationHeatFlux(
    Adapt.adapt(to, f.conductive_flux),
    Adapt.adapt(to, f.internal_temperature),
    f.heat_transfer_coefficient,
    f.salt_transfer_coefficient,
    Adapt.adapt(to, f.friction_velocity),
)

"""
    ThreeEquationHeatFlux(sea_ice::Simulation{<:SeaIceModel}, FT::DataType = Oceananigans.defaults.FloatType;
                          heat_transfer_coefficient = 0.0095,
                          salt_transfer_coefficient = heat_transfer_coefficient / 35,
                          friction_velocity = 0.002)

Construct a `ThreeEquationHeatFlux` whose conductive flux and internal temperature are taken from the
thermodynamics of `sea_ice`.

Default values follow Shi et al. (2021) with ``R = \\alpha_h / \\alpha_s = 35``.

Keyword Arguments
=================

- `heat_transfer_coefficient`: turbulent heat exchange coefficient ``\\alpha_h``. Default: 0.0095.
- `salt_transfer_coefficient`: turbulent salt exchange coefficient ``\\alpha_s``. Default: ``\\alpha_h / 35 \\approx 0.000271``.
- `friction_velocity`: friction velocity value or formulation. Default: 0.002.
"""
function ThreeEquationHeatFlux(
    sea_ice::OC.Simulation{<:CSI.SeaIceModel},
    FT::DataType = OC.defaults.FloatType;
    heat_transfer_coefficient = 0.0095,
    salt_transfer_coefficient = heat_transfer_coefficient / 35,
    friction_velocity = convert(FT, 0.002),
)

    conductive_flux = sea_ice.model.ice_thermodynamics.internal_heat_flux
    ice_temperature = sea_ice.model.ice_thermodynamics.top_surface_temperature

    return ThreeEquationHeatFlux(
        conductive_flux,
        ice_temperature,
        convert(FT, heat_transfer_coefficient),
        convert(FT, salt_transfer_coefficient),
        friction_velocity,
    )
end

@inline extract_internal_temperature(flux::ThreeEquationHeatFlux, i, j) =
    @inbounds flux.internal_temperature[i, j, 1]

@inline function store_interface_state!(::ThreeEquationHeatFlux, T★, S★, i, j, Tᵦ, Sᵦ)
    @inbounds T★[i, j, 1] = Tᵦ
    @inbounds S★[i, j, 1] = Sᵦ
end

"""
    compute_interface_heat_flux(flux::ThreeEquationHeatFlux, ocean_state, ice_state, liquidus, ocean_properties, ℰ, u★)

Compute the heat flux and melt rate at the sea ice-ocean interface using the three-equation formulation.

Returns `(Q, q, Tᵦ, Sᵦ)` where:
- `Q > 0` means heat flux from ocean to ice (ocean cooling)
- `q > 0` means melting (ice volume loss)
- `Tᵦ, Sᵦ` are the interface temperature and salinity
"""
@inline function compute_interface_heat_flux(
    flux::ThreeEquationHeatFlux,
    ocean_state,
    ice_state,
    liquidus,
    ocean_properties,
    ℰ,
    u★,
)
    Tᵒᶜ = ocean_state.T
    Sᵒᶜ = ocean_state.S
    ℵ = ice_state.ℵ

    ρᵒᶜ = ocean_properties.reference_density
    cᵒᶜ = ocean_properties.heat_capacity

    αₕ = flux.heat_transfer_coefficient
    αₛ = flux.salt_transfer_coefficient

    T★, S★, q = solve_interface_conditions(
        flux,
        Tᵒᶜ,
        Sᵒᶜ,
        ice_state,
        αₕ,
        αₛ,
        u★,
        ℰ,
        ρᵒᶜ,
        cᵒᶜ,
        liquidus,
    )

    # Scale by ice concentration
    q = q * ℵ
    Qᵢₒ = ℰ * q

    return Qᵢₒ, q, T★, S★
end

@inline function conductive_flux_parameters(flux::ThreeEquationHeatFlux, ice_state, ℰ)
    h = ice_state.h
    hc = ice_state.hc
    Tˢⁱ = ice_state.T
    k = flux.conductive_flux.conductivity
    # Set κ to zero when h < hc (ice not consolidated)
    consolidated = h ≥ hc
    κ = ifelse(consolidated, k / (h * ℰ), zero(h))
    return κ, Tˢⁱ
end

"""
    solve_interface_conditions(flux::ThreeEquationHeatFlux, Tᵒᶜ, Sᵒᶜ, ice_state, αₕ, αₛ, u★, ℰ, ρᵒᶜ, cᵒᶜ, liquidus)

Solve the three-equation system for interface temperature, salinity, and melt rate.

The three equations are:
1. Heat balance: ``ρᵒᶜ cᵒᶜ αₕ u★ (Tᵒᶜ - T★) + κ (Tˢⁱ - T★) = ℰ q``
2. Salt balance: ``ρᵒᶜ αₛ u★ (Sᵒᶜ - S★) = q (S★ - Sˢⁱ)``
3. Freezing point: ``T★ = Tₘ(S★)``

where `κ = k/(h ℰ)` is the conductive heat transfer coefficient, zero where the ice is not consolidated.

Arguments
=========
- `ice_state`: NamedTuple with fields `S`, `h`, `hc`, `ℵ`, `T` (internal temperature)

Returns `(T★, S★, q)` where q is the melt rate (positive for melting).
"""
@inline function solve_interface_conditions(
    flux::ThreeEquationHeatFlux,
    Tᵒᶜ,
    Sᵒᶜ,
    ice_state,
    αₕ,
    αₛ,
    u★,
    ℰ,
    ρᵒᶜ,
    cᵒᶜ,
    liquidus::LinearLiquidus,
)
    Sˢⁱ = ice_state.S

    κ, Tˢⁱ = conductive_flux_parameters(flux, ice_state, ℰ)

    λ₁ = -liquidus.slope
    λ₂ = liquidus.freshwater_melting_temperature

    # Transfer coefficients
    η = ρᵒᶜ * cᵒᶜ * αₕ * u★ / ℰ  # turbulent heat
    γ = ρᵒᶜ * αₛ * u★           # turbulent salt
    θ = η + κ                  # total heat

    # Quadratic coefficients: a S★² + b S★ + c = 0
    a = θ * λ₁
    b = -γ - η * Tᵒᶜ - κ * Tˢⁱ + θ * (λ₂ - λ₁ * Sˢⁱ)
    c = γ * Sᵒᶜ + (η * Tᵒᶜ + κ * Tˢⁱ - θ * λ₂) * Sˢⁱ

    # Solve quadratic with zero-safe reciprocal (MITgcm approach)
    ξ = ifelse(a == zero(a), zero(a), one(a) / (2a))
    Δ = max(b^2 - 4a * c, zero(a))
    S★ = (-b - sqrt(Δ)) * ξ
    S★ = ifelse(S★ < zero(S★), (-b + sqrt(Δ)) * ξ, S★)

    T★ = melting_temperature(liquidus, S★)

    q = η * (Tᵒᶜ - T★) + κ * (Tˢⁱ - T★)

    return T★, S★, q
end

Base.summary(
    ::ThreeEquationHeatFlux{<:Any, <:Any, FT},
) where {FT} = "ThreeEquationHeatFlux{$FT}"

function Base.show(io::IO, flux::ThreeEquationHeatFlux)
    print(io, summary(flux), '\n')
    print(io, "├── heat_transfer_coefficient: ", flux.heat_transfer_coefficient, '\n')
    print(io, "├── salt_transfer_coefficient: ", flux.salt_transfer_coefficient, '\n')
    print(io, "└── friction_velocity: ", flux.friction_velocity)
end

#####
##### Utilities
#####

function interface_kernel_parameters(grid)
    Sx, Sy, _ = worksize(grid)
    TX, TY, _ = topology(grid)
    single_column_grid = Sx == 1 && Sy == 1

    if single_column_grid
        kernel_parameters = KernelParameters(1:1, 1:1)
    else
        # Compute fluxes into halo regions (0:N+1) for non-Flat dimensions.
        # Flat dimensions have no halo cells, so only iterate over the interior.
        x_range = TX === Flat ? (1:Sx) : (0:(Sx + 1))
        y_range = TY === Flat ? (1:Sy) : (0:(Sy + 1))
        kernel_parameters = KernelParameters(x_range, y_range)
    end

    return kernel_parameters
end

#####
##### Sea ice-ocean fluxes
#####

"""
    compute_sea_ice_ocean_fluxes!(interface, ocean, sea_ice, ocean_properties; Δt)

Compute heat, salt, and momentum fluxes at the sea ice-ocean interface over the time interval `Δt`.

This function computes:
- Frazil heat flux: heat released when ocean temperature drops below freezing
- Interface heat flux: heat flux from ocean to ice, computed with `interface.flux_formulation`
- Salt flux: salt exchange due to ice growth/melt
- Momentum stresses: ice-ocean momentum transfer
"""
function compute_sea_ice_ocean_fluxes!(interface, ocean, sea_ice, ocean_properties; Δt)
    Tᵒᶜ = ocean.model.tracers.T
    Sᵒᶜ = ocean.model.tracers.S
    Sⁱ = sea_ice.model.tracers.S
    ℵ = sea_ice.model.ice_concentration
    hˢⁱ = sea_ice.model.ice_thickness
    hc = sea_ice.model.ice_consolidation_thickness

    phase_transitions = sea_ice.model.phase_transitions
    liquidus = phase_transitions.liquidus
    L = phase_transitions.reference_latent_heat

    grid = sea_ice.model.grid
    clock = sea_ice.model.clock
    arch = architecture(grid)

    uˢⁱ, vˢⁱ = sea_ice.model.velocities
    dynamics = sea_ice.model.dynamics

    fluxes = interface.fluxes
    flux_formulation = interface.flux_formulation
    Tˢⁱ = interface.temperature
    Sˢⁱ = interface.salinity

    if !isnothing(dynamics)
        kernel_parameters = interface_kernel_parameters(grid)
        τₛ = dynamics.external_momentum_stresses.bottom
        launch!(
            arch,
            grid,
            kernel_parameters,
            _compute_sea_ice_ocean_stress!,
            fluxes,
            grid,
            clock,
            hˢⁱ,
            ℵ,
            uˢⁱ,
            vˢⁱ,
            τₛ,
        )
    else
        τₛ = nothing
    end

    launch!(
        arch,
        grid,
        :xy,
        _compute_sea_ice_ocean_fluxes!,
        flux_formulation,
        fluxes,
        Tˢⁱ,
        Sˢⁱ,
        grid,
        clock,
        hˢⁱ,
        hc,
        ℵ,
        Sⁱ,
        Tᵒᶜ,
        Sᵒᶜ,
        uˢⁱ,
        vˢⁱ,
        τₛ,
        liquidus,
        ocean_properties,
        L,
        Δt,
    )

    return nothing
end

@kernel function _compute_sea_ice_ocean_stress!(
    fluxes,
    grid,
    clock,
    ice_thickness,
    ice_concentration,
    sea_ice_u_velocity,
    sea_ice_v_velocity,
    sea_ice_ocean_stress,
)
    i, j = @index(Global, NTuple)

    τˣ = fluxes.x_momentum
    τʸ = fluxes.y_momentum
    Nz = size(grid, 3)

    uˢⁱ = sea_ice_u_velocity
    vˢⁱ = sea_ice_v_velocity
    hˢⁱ = ice_thickness
    ℵ = ice_concentration
    sea_ice_fields = (; u = uˢⁱ, v = vˢⁱ, h = hˢⁱ, ℵ = ℵ)

    # Momentum stresses
    @inbounds begin
        τˣ[i, j, 1] =
            x_momentum_stress(i, j, Nz, grid, sea_ice_ocean_stress, clock, sea_ice_fields)
        τʸ[i, j, 1] =
            y_momentum_stress(i, j, Nz, grid, sea_ice_ocean_stress, clock, sea_ice_fields)
    end
end

@kernel function _compute_sea_ice_ocean_fluxes!(
    flux_formulation,
    fluxes,
    interface_temperature,
    interface_salinity,
    grid,
    clock,
    ice_thickness,
    ice_consolidation_thickness,
    ice_concentration,
    ice_salinity,
    ocean_temperature,
    ocean_salinity,
    sea_ice_u_velocity,
    sea_ice_v_velocity,
    sea_ice_ocean_stresses,
    liquidus,
    ocean_properties,
    latent_heat,
    Δt,
)

    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)
    𝒬ᶠʳᶻ = fluxes.frazil_heat
    𝒬ⁱⁿ = fluxes.interface_heat
    Jˢ = fluxes.salt
    τˣ = fluxes.x_momentum
    τʸ = fluxes.y_momentum
    T★ = interface_temperature
    S★ = interface_salinity
    Tᵒᶜ = ocean_temperature
    Sᵒᶜ = ocean_salinity
    hc = ice_consolidation_thickness
    ℰ = latent_heat

    ρᵒᶜ = ocean_properties.reference_density
    cᵒᶜ = ocean_properties.heat_capacity

    # =============================================
    # Part 1: Frazil ice formation
    # =============================================
    # When ocean temperature drops below freezing, frazil ice forms
    # and heat is released to the ice component.

    δ𝒬ᶠʳᶻ = zero(grid)

    for k in Nz:-1:1
        @inbounds begin
            Δz = Δzᶜᶜᶜ(i, j, k, grid)
            Tᵏ = Tᵒᶜ[i, j, k]
            Sᵏ = Sᵒᶜ[i, j, k]
        end

        # Melting/freezing temperature at this depth
        Tₘ = melting_temperature(liquidus, Sᵏ)
        freezing = Tᵏ < Tₘ

        # Compute change in ocean heat energy due to freezing.
        # When Tᵏ < Tₘ, we heat the ocean back to melting temperature
        # by extracting heat from the ice.
        δE = freezing * ρᵒᶜ * cᵒᶜ * (Tₘ - Tᵏ)

        # Perform temperature adjustment
        @inbounds Tᵒᶜ[i, j, k] = ifelse(freezing, Tₘ, Tᵏ)

        # Compute the heat flux from ocean into ice during frazil formation.
        # A negative value δ𝒬ᶠʳᶻ < 0 implies heat is fluxed from the ice into
        # the ocean (frazil ice formation).
        δ𝒬ᶠʳᶻ -= δE * Δz / Δt
    end

    # Store frazil heat flux
    @inbounds 𝒬ᶠʳᶻ[i, j, 1] = δ𝒬ᶠʳᶻ

    # Freezing rate
    qᶠ = δ𝒬ᶠʳᶻ / ℰ

    @inbounds begin
        Tᴺ = Tᵒᶜ[i, j, Nz]
        Sᴺ = Sᵒᶜ[i, j, Nz]
        Sˢⁱ = ice_salinity[i, j, 1]
        hˢⁱ = ice_thickness[i, j, 1]
        ℵᵢ = ice_concentration[i, j, 1]
        hc = ice_consolidation_thickness[i, j, 1]
    end

    Tˢⁱ = extract_internal_temperature(flux_formulation, i, j)

    # Package states
    ocean_surface_state = (; T = Tᴺ, S = Sᴺ)
    ice_state = (; S = Sˢⁱ, h = hˢⁱ, hc = hc, ℵ = ℵᵢ, T = Tˢⁱ)

    # Compute friction velocity
    u★ = get_friction_velocity(flux_formulation.friction_velocity, i, j, grid, τˣ, τʸ, ρᵒᶜ)

    # =============================================
    # Part 3: Interface heat flux
    # =============================================
    # Returns interfacial heat flux, melt rate qᵐ, and interface T, S
    𝒬ⁱᵒ, qᵐ, Tᵦ, Sᵦ = compute_interface_heat_flux(
        flux_formulation,
        ocean_surface_state,
        ice_state,
        liquidus,
        ocean_properties,
        ℰ,
        u★,
    )

    # Store interface values and heat flux
    @inbounds 𝒬ⁱⁿ[i, j, 1] = 𝒬ⁱᵒ
    store_interface_state!(flux_formulation, T★, S★, i, j, Tᵦ, Sᵦ)

    # =============================================
    # Part 4: Salt flux
    # =============================================
    # Salt flux from melting/freezing:
    # - during ice melt   (qᵐ > 0), fresh meltwater dilutes the ocean
    # - during ice growth (qᶠ < 0), brine rejection adds salt to ocean
    @inbounds Jˢ[i, j, 1] = (qᵐ + qᶠ) / ρᵒᶜ * (Sᴺ - Sˢⁱ)
end

#####
##### Above-freezing ocean temperature
#####

@kernel function _above_freezing_ocean_temperature!(T, grid, S, liquidus)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        for k in 1:Nz
            Tm = melting_temperature(liquidus, S[i, j, k])
            T[i, j, k] = max(T[i, j, k], Tm)
        end
    end
end

"""
    above_freezing_ocean_temperature!(ocean, grid, sea_ice)

Raise the ocean temperature to the melting temperature of the local salinity wherever it lies below it.
"""
function above_freezing_ocean_temperature!(ocean, grid, sea_ice)
    T = ocean.model.tracers.T
    S = ocean.model.tracers.S
    liquidus = sea_ice.model.phase_transitions.liquidus

    arch = architecture(grid)
    launch!(arch, grid, :xy, _above_freezing_ocean_temperature!, T, grid, S, liquidus)

    return nothing
end

above_freezing_ocean_temperature!(ocean, grid, ::Nothing) = nothing
