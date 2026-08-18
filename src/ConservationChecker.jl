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
extra_terms(::TotalWater) = (:runoff, :lai_leakage)

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
        keys = (sim_names..., extra_terms(cq)...)
        components = NamedTuple{keys}(ntuple(_ -> Float64[], length(keys))) # just using Float64 here because it's a diagnostic
        return new{typeof(cq), typeof(components)}(cq, components)
    end
end

# Compute global integrals of `ConservedQuantity` currently held by `sim`
component_contribution(::TotalEnergy, sim::Interfacer.AbstractComponentSimulation) =
    Interfacer.get_field(sim, Val(:total_energy)) # J
component_contribution(::TotalWater, sim::Interfacer.AbstractComponentSimulation) =
    Interfacer.get_field(sim, Val(:total_water)) # kg

"""
    sum_over_sims(cs::Interfacer.CoupledSimulation, val::Val)

Sum `val` over every component simulation that reports it, treating components that
return `nothing` as contributing zero.
"""
function sum_over_sims(cs::Interfacer.CoupledSimulation, val::Val)
    total = 0.0
    for sim in cs.model_sims
        contribution = Interfacer.get_field(sim, val)
        isnothing(contribution) || (total += contribution)
    end
    return total
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
    # TOA radiation imbalance accumulation over a timestep
    toa_rad_accum =
        isempty(components.toa_net_source) ? 0.0 : components.toa_net_source[end]
    toa_rad_accum += float(cs.Δt_cpl) * sum_over_sims(cs, Val(:radiative_energy_flux_toa)) # J
    push!(components.toa_net_source, toa_rad_accum)
    return nothing
end

function update_extra_terms!(::TotalWater, components, cs::Interfacer.CoupledSimulation)
    # Runoff is not currently sent to the ocean :( we accumulate it here
    runoff_accum = isempty(components.runoff) ? 0.0 : components.runoff[end]
    runoff_accum += float(cs.Δt_cpl) * sum_over_sims(cs, Val(:runoff)) # kg
    push!(components.runoff, runoff_accum)

    # Water that the land model aquires due to a prescribed LAI (this is already accumulated by the getter)
    lai_leakage = sum_over_sims(cs, Val(:lai_leakage)) # kg
    push!(components.lai_leakage, -lai_leakage) # minus sign because this is a source for the land model but a sink for the coupled system
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

    for sim in cs.model_sims
        contribution = component_contribution(cq, sim)
        push!(
            getproperty(components, Symbol(Base.nameof(sim))),
            isnothing(contribution) ? 0.0 : contribution,
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
