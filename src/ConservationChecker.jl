"""
    ConservationChecker

This module contains machinery for checking that conserved quantities are in fact conserved.
"""
module ConservationChecker

import ..Interfacer, ..Utilities
import ..Utilities: integral
import ClimaCore as CC

export ConservedQuantity, TotalWater, TotalEnergy, ConservationCheck, check_conservation!

abstract type ConservedQuantity end

struct TotalWater <: ConservedQuantity end
struct TotalEnergy <: ConservedQuantity end

name(::TotalWater) = "water"
units(::TotalWater) = "kg"
extra_terms(::TotalWater) = (:runoff_net_source, :nonflux_water_source)

name(::TotalEnergy) = "energy"
units(::TotalEnergy) = "J"
extra_terms(::TotalEnergy) = (:toa_net_source,)

"""
    ConservationCheck(cq::ConservedQuantity, model_sims)

A log of how much of `cq` the coupled system holds, with one timeseries per
component simulation plus one per entry of `extra_terms(cq)`.
"""
struct ConservationCheck{CQ <: ConservedQuantity, C}
    conserved_quantity::CQ
    components::C

    function ConservationCheck(cq::ConservedQuantity, model_sims)
        sim_names = map(sim -> Symbol(Base.nameof(sim)), Tuple(model_sims))
        ks = (sim_names..., extra_terms(cq)...)
        components = NamedTuple{ks}(ntuple(_ -> Float64[], length(ks)))
        return new{typeof(cq), typeof(components)}(cq, components)
    end
end

"""
    area_weighted_integral(field, area_fraction)

Integrate `field` over the boundary space, weighted by `area_fraction`, ignoring
the contribution from points that the surface does not cover.
"""
function area_weighted_integral(field, area_fraction)
    FT = eltype(area_fraction)
    return integral(area_fraction .* ifelse.(area_fraction .≈ 0, zero(FT), field))
end

function preprocess!(::TotalWater, cs::Interfacer.CoupledSimulation)
    # update coupler net precipitation (for surfaces that don't collect water)
    E = cs.fields.F_turb_moisture # kg / m^2 / s / layer depth
    FT = eltype(E)
    P_liq = cs.fields.P_liq
    P_snow = cs.fields.P_snow
    P_net = @. -(E + P_liq + P_snow) * FT(float(cs.Δt_cpl)) # kg / m^2 / layer depth
    Interfacer.remap!(cs.fields.scalar_temp1, P_net)
    cs.fields.P_net .+= cs.fields.scalar_temp1
    return nothing
end
preprocess!(::ConservedQuantity, ::Interfacer.CoupledSimulation) = nothing

"""
    component_contribution(cq::ConservedQuantity,
                           sim::Interfacer.AbstractComponentSimulation,
                           cs::Interfacer.CoupledSimulation)

Compute global integral of `cq` currently held by `sim`.
"""
component_contribution(
    ::ConservedQuantity,
    ::Interfacer.AbstractComponentSimulation,
    ::Interfacer.CoupledSimulation,
) = 0

function component_contribution(
    ::TotalEnergy,
    sim::Interfacer.AbstractAtmosSimulation,
    ::Interfacer.CoupledSimulation,
)
    return integral(Interfacer.get_field(sim, Val(:energy))) # J (∫ J / m^3 dV)
end

function component_contribution(
    ::TotalEnergy,
    sim::Interfacer.AbstractSurfaceSimulation,
    cs::Interfacer.CoupledSimulation,
)
    isnothing(Interfacer.get_field(sim, Val(:energy))) && return 0
    # regrid onto the boundary space before integrating
    return area_weighted_integral(
        Interfacer.get_field(Interfacer.boundary_space(cs), sim, Val(:energy)),
        Interfacer.get_field(sim, Val(:area_fraction)),
    ) # J (∫ J / m^2 dA)
end

function component_contribution(
    ::TotalWater,
    sim::Interfacer.AbstractAtmosSimulation,
    ::Interfacer.CoupledSimulation,
)
    water = Interfacer.get_field(sim, Val(:water))
    water isa Number && return 0 # Case with no moisture. TODO: Handle this better
    return integral(water) # kg (∫ kg of water / m^3 dV)
end

function component_contribution(
    ::TotalWater,
    sim::Interfacer.AbstractSurfaceSimulation,
    cs::Interfacer.CoupledSimulation,
)
    area_fraction = Interfacer.get_field(sim, Val(:area_fraction))
    if isnothing(Interfacer.get_field(sim, Val(:water)))
        # surfaces that don't collect water hold the net precipitation instead
        return area_weighted_integral(cs.fields.P_net, area_fraction) # kg (∫ kg / m^2 dA)
    end
    return area_weighted_integral(
        Interfacer.get_field(Interfacer.boundary_space(cs), sim, Val(:water)),
        area_fraction,
    ) # kg (∫ kg / m^2 dA)
end

"""
    update_extra_terms!(cq::ConservedQuantity,
                        components,
                        cs::Interfacer.CoupledSimulation)

Append one value to each of the `extra_terms(cq)` timeseries in `components`.

These terms are collected for bookkeeping purposes to close the budget; ideally
we should not need them!
"""
function update_extra_terms!(::TotalEnergy, components, cs::Interfacer.CoupledSimulation)
    FT = CC.Spaces.undertype(Interfacer.boundary_space(cs))

    # TOA radiation imbalance
    accum = isempty(components.toa_net_source) ? FT(0) : components.toa_net_source[end]
    for sim in cs.model_sims
        sim isa Interfacer.AbstractAtmosSimulation || continue
        radiative_energy_flux_toa =
            Interfacer.get_field(sim, Val(:radiative_energy_flux_toa))
        radiative_energy_flux_toa isa Number && continue # Case with no radiation. TODO: Handle this better
        accum += integral(radiative_energy_flux_toa) * FT(float(cs.Δt_cpl)) # ∫ J / m^2 dA
    end
    push!(components.toa_net_source, accum)
    return nothing
end

function update_extra_terms!(::TotalWater, components, cs::Interfacer.CoupledSimulation)
    FT = CC.Spaces.undertype(Interfacer.boundary_space(cs))

    # Runoff is not currently sent to the ocean :( we accumulate it here
    runoff_accum =
        isempty(components.runoff_net_source) ? FT(0) : components.runoff_net_source[end]

    # Water that components report holding but that never arrived as a flux
    # (this happens in the land model becuase LAI is prescribed)
    # (no need to accumulate this, just store the current value to close the budget)
    nonflux_water = FT(0)

    for sim in cs.model_sims
        sim isa Interfacer.AbstractSurfaceSimulation || continue
        area_fraction = Interfacer.get_field(sim, Val(:area_fraction))

        # accumulate the water this surface has shed as runoff
        if !isnothing(Interfacer.get_field(sim, Val(:runoff)))
            runoff_accum +=
                area_weighted_integral(
                    Interfacer.get_field(Interfacer.boundary_space(cs), sim, Val(:runoff)),
                    area_fraction,
                ) * FT(float(cs.Δt_cpl)) # kg (∫ kg / m^2 / s dA dt)
        end

        # collect any water this surface reports that isn't flux-driven
        if !isnothing(Interfacer.get_field(sim, Val(:nonflux_water)))
            nonflux_water += area_weighted_integral(
                Interfacer.get_field(
                    Interfacer.boundary_space(cs),
                    sim,
                    Val(:nonflux_water),
                ),
                area_fraction,
            ) # kg (∫ kg / m^2 dA)
        end
    end

    push!(components.runoff_net_source, runoff_accum)
    push!(components.nonflux_water_source, -nonflux_water) # note the minus!
    return nothing
end

"""
    check_conservation!(coupler_sim::Interfacer.CoupledSimulation; runtime_check = false)

Iterates over all specified conservation checks and returns values.
"""
function check_conservation!(
    coupler_sim::Interfacer.CoupledSimulation;
    runtime_check = false,
)
    isnothing(coupler_sim.conservation_checks) && return nothing
    return map(
        x -> check_conservation!(x, coupler_sim, runtime_check),
        coupler_sim.conservation_checks,
    )
end

"""
    check_conservation!(cc::ConservationCheck,
                        cs::Interfacer.CoupledSimulation,
                        runtime_check = false)

Compute and append global integral of `cc.conserved_quantity` to every timeseries in `cc.components`.
"""
function check_conservation!(
    cc::ConservationCheck,
    cs::Interfacer.CoupledSimulation,
    runtime_check = false,
)
    cq = cc.conserved_quantity
    components = cc.components

    preprocess!(cq, cs)

    for sim in cs.model_sims
        push!(
            getproperty(components, Symbol(Base.nameof(sim))),
            component_contribution(cq, sim, cs),
        )
    end
    update_extra_terms!(cq, components, cs)

    total = sum(last, values(components))
    if runtime_check
        total_initial = sum(first, values(components))
        @assert abs((total - total_initial) / total) < 1e-4
    end
    return total
end

end # module
