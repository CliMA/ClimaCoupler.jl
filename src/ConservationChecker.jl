"""
    ConservationChecker

This module contains machinery for checking that conserved quantities are in fact conserved.
"""
module ConservationChecker

import ..Interfacer, ..Utilities
import ..Utilities: integral
import ClimaCore as CC

export ConservedQuantity,
    TotalWater,
    TotalEnergy,
    ConservationCheck,
    contributions,
    Accumulated,
    check_conservation!

abstract type ConservedQuantity end

struct TotalWater <: ConservedQuantity end
struct TotalEnergy <: ConservedQuantity end

name(::TotalWater) = "water"
units(::TotalWater) = "kg"

name(::TotalEnergy) = "energy"
units(::TotalEnergy) = "J"

"""
    Accumulated(rate)

Wrapper marking a contribution that is a rate rather than an amount, so that the
checker integrates it in time instead of logging it directly.
"""
struct Accumulated{T}
    rate::T
end

"""
    contributions(cq::ConservedQuantity, sim::Interfacer.AbstractComponentSimulation)

Return everything `sim` contributes to the budget of `cq`, as a `NamedTuple` of named
terms (default: `(;)`, no contribution).

Every term is expressed so that **the sum of all terms over all components is the
total**, in the units of `cq`. A term describing something that has left the coupled
system is therefore reported with the sign that adds it back, and a term describing
something that appeared from outside is reported negated.

Terms are either amounts, used as they are, or an [`Accumulated`](@ref) rate, which
the checker integrates in time. By convention the amount a component holds is named
`reservoir` and is logged under the component's own name; any other term `x` reported
by component `C` is logged as `C_x`.

For example, the integrated land model reports the water it holds, the runoff rate
leaving the coupled system, and the water its prescribed LAI has conjured up:

```julia
ConservationChecker.contributions(::TotalWater, sim::ClimaLandSimulation) = (;
    reservoir = <kg held>,
    runoff = Accumulated(<kg / s leaving>),
    lai_leakage = <-kg gained from the prescribed LAI>,
)
```
"""
contributions(::ConservedQuantity, ::Interfacer.AbstractComponentSimulation) = (;)

"""
    term_key(sim, term::Symbol)

Return the key that a component's term is logged under: the component's own name for
its `reservoir`, and `<component>_<term>` for anything else.
"""
term_key(sim, term::Symbol) =
    term === :reservoir ? Symbol(Base.nameof(sim)) : Symbol(Base.nameof(sim), :_, term)

"""
    ConservationCheck(cq::ConservedQuantity, model_sims)

A log of how much of `cq` the coupled system holds, with one timeseries per term
reported by the component simulations. See [`contributions`](@ref).
"""
struct ConservationCheck{CQ <: ConservedQuantity, C}
    conserved_quantity::CQ
    components::C

    function ConservationCheck(cq::ConservedQuantity, model_sims)
        keys = Symbol[]
        for sim in model_sims, term in propertynames(contributions(cq, sim))
            push!(keys, term_key(sim, term))
        end
        # just using Float64 here because it's a diagnostic
        components = NamedTuple{Tuple(keys)}(ntuple(_ -> Float64[], length(keys)))
        return new{typeof(cq), typeof(components)}(cq, components)
    end
end

"""
    log_value(contribution, timeseries, Δt_cpl)

Return the value to append to `timeseries` for `contribution`.

Amounts are logged as they are; an [`Accumulated`](@ref) rate is integrated onto the
running total, so that its timeseries holds the time integral of the rate.
"""
log_value(contribution, _timeseries, _Δt_cpl) = contribution
log_value(contribution::Accumulated, timeseries, Δt_cpl) =
    (isempty(timeseries) ? 0.0 : timeseries[end]) + float(Δt_cpl) * contribution.rate

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

Append one value to every timeseries in `cc.components`, and return the total amount
of `cc.conserved_quantity` held by the coupled system.
"""
function check_conservation!(
    cc::ConservationCheck,
    cs::Interfacer.CoupledSimulation,
    runtime_check = false,
)
    cq = cc.conserved_quantity
    components = cc.components

    for sim in cs.model_sims
        for (term, contribution) in pairs(contributions(cq, sim))
            timeseries = getproperty(components, term_key(sim, term))
            push!(timeseries, log_value(contribution, timeseries, cs.Δt_cpl))
        end
    end

    # `init` so that a quantity no component reports totals to zero rather than erroring
    total = sum(last, values(components); init = 0.0)
    if runtime_check
        total_initial = sum(first, values(components); init = 0.0)
        @assert abs((total - total_initial) / total) < 1e-4
    end
    return total
end

end # module
