
Base.@kwdef struct FluxAccumulator{B} <: AbstractPhysicsComponent
    bctype::B
end

function calc_component!(source, accumulator::FluxAccumulator, state, aux, physics)

    # save flux in atmos source for integration in time using atmos timestepping
    if boundary_mask(physics.parameters, aux.x, aux.y, aux.z)
        E, H = calc_ocean_sfc_fluxes(physics, accumulator.bctype, state, aux)
        source.F_ρe_accum = (E + H) # latent + sensible heat fluxes [W/m^2]
    end
end


"""
    calculate_land_sfc_fluxes(model::DryAtmosModel, state, aux)
- calculate furface fluxes using the bulk gradient diffusion theory
"""

function calc_ocean_sfc_fluxes(physics, bctype::BulkFormulaTemperature, state⁻, aux⁻)
    # Apply bulks laws using the tangential velocity as energy flux
    ρ = state⁻.ρ
    ρu = state⁻.ρu
    ρq = state⁻.ρq
    eos = physics.eos
    parameters = physics.parameters

    # vertical unit vector
    n̂ = aux⁻.∇Φ / parameters.g

    # obtain surface fields from bcs
    ϕ = lat(aux⁻.x, aux⁻.y, aux⁻.z)
    Cₕ = bctype.drag_coef_temperature(parameters, ϕ)
    Cₑ = bctype.drag_coef_moisture(parameters, ϕ)
    T_sfc = bctype.surface_temperature(parameters, ϕ)
    LH_v0 = parameters.LH_v0

    # magnitude of tangential velocity (usually called speed)
    u = ρu / ρ
    speed_tangential = norm((I - n̂ ⊗ n̂) * u)

    # sensible heat flux
    cp = calc_heat_capacity_at_constant_pressure(eos, state⁻, parameters)
    T = calc_air_temperature(eos, state⁻, aux⁻, parameters)
    H = ρ * Cₕ * speed_tangential * cp * (T_sfc - T)

    # latent heat flux
    q = ρq / ρ
    q_tot_sfc = calc_saturation_specific_humidity(ρ, T_sfc, parameters)
    E = ρ * Cₑ * speed_tangential * LH_v0 * (q_tot_sfc - q)

    return E, H
end

function calc_ocean_sfc_fluxes(physics, bctype::CoupledPrimarySlabOceanBC, state⁻, aux⁻; MO_params = nothing) # should pass in the coupler state (also move to coupler), so can access states of both models derectly -e.g. callback?

    # Apply bulks laws using the tangential velocity as energy flux
    ρ = state⁻.ρ
    ρu = state⁻.ρu
    ρq = state⁻.ρq
    eos = physics.eos
    parameters = physics.parameters

    # vertical unit vector
    n̂ = aux⁻.∇Φ / parameters.g

    # obtain surface fields from bcs
    Cₕ = parameters.Cₗ
    Cₑ = parameters.Cₑ
    LH_v0 = parameters.LH_v0
    T_sfc = aux⁻.T_sfc

    # magnitude of tangential velocity (usually called speed)
    u = ρu / ρ
    speed_tangential = norm((I - n̂ ⊗ n̂) * u)

    # sensible heat flux
    cp = calc_heat_capacity_at_constant_pressure(eos, state⁻, parameters)
    T = calc_air_temperature(eos, state⁻, aux⁻, parameters)
    H = ρ * Cₕ * speed_tangential * cp * (T_sfc - T)

    # latent heat flux
    q = ρq / ρ
    q_tot_sfc = calc_saturation_specific_humidity(ρ, T_sfc, parameters)
    E = ρ * Cₑ * speed_tangential * LH_v0 * (q_tot_sfc - q)

    E = isnan(E) ? zero(E) : E # set the fluxes of the non-sfc layers to 0 (otherwise NaNs if moisture) TODO: if calc_component imports BL, then could apply the FluxAccumulator on the boundary only and this line won't be needed 

    return E, H
end
