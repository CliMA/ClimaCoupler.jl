"""
    ConservationChecker

This module contains functions that check global conservation of energy and water (and momentum - TODO).
"""
module ConservationChecker

using ClimaCore: ClimaCore, Geometry, Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaComms
using NCDatasets
using ClimaCoreTempestRemap
using Dates
using JLD2
using Plots
using ClimaAtmos: RRTMGPI
using ClimaLSM
using ClimaCoupler.Utilities: swap_space!
import ClimaCoupler: Interfacer

export AbstractConservationCheck,
    EnergyConservationCheck, WaterConservationCheck, check_conservation!, plot_global_conservation

abstract type AbstractConservationCheck end

"""
    EnergyConservationCheck{A} <: AbstractConservationCheck

Struct of type `AbstractConservationCheck` containing global energy conservation logs.
"""
mutable struct EnergyConservationCheck <: AbstractConservationCheck
    sums::NamedTuple
    function EnergyConservationCheck(sums)
        all_sims = (;)
        for sim in sums
            all_sims = merge(all_sims, [Symbol(Interfacer.name(sim)) => []])
        end
        all_sims = (all_sims..., toa_net_source = [], total = [])
        return new(all_sims)
    end
end
Interfacer.name(::EnergyConservationCheck) = "energy [J]"

"""
    WaterConservationCheck{A} <: AbstractConservationCheck

Struct of type `AbstractConservationCheck` containing global water mass conservation logs.
"""
mutable struct WaterConservationCheck <: AbstractConservationCheck
    sums::NamedTuple
    function WaterConservationCheck(sums)
        all_sims = (;)
        for sim in sums
            all_sims = merge(all_sims, [Symbol(Interfacer.name(sim)) => []])
        end
        all_sims = (all_sims..., total = [])
        return new(all_sims)
    end
end
Interfacer.name(::WaterConservationCheck) = "water [kg]"

"""
    check_conservation!(coupler_sim::Interfacer.CoupledSimulation; runtime_check = false)

itertes over all specified conservation checks.
"""
check_conservation!(coupler_sim::Interfacer.CoupledSimulation; runtime_check = false) =
    map(x -> check_conservation!(x, coupler_sim, runtime_check), coupler_sim.conservation_checks)

"""
        check_conservation!(
        cc::EnergyConservationCheck,
        coupler_sim::Interfacer.CoupledSimulation,
        runtime_check = false,
        )

computes the total energy, ∫ ρe dV, of the model components
of the coupled simulations and the TOA radiation, and updates
`cc` with these values.
"""
function check_conservation!(
    cc::EnergyConservationCheck,
    coupler_sim::Interfacer.CoupledSimulation,
    runtime_check = false,
)

    ccs = cc.sums
    (; model_sims) = coupler_sim

    boundary_space = coupler_sim.boundary_space # thin shell approx (boundary_space[z=0] = boundary_space[z_top])

    FT = eltype(coupler_sim.fields[1])

    total = 0

    # save surfaces
    for sim in model_sims
        sim_name = Symbol(Interfacer.name(sim))
        if sim isa Interfacer.AtmosModelSimulation
            # save radiation source
            parent(coupler_sim.fields.F_radiative_TOA) .= parent(Interfacer.get_field(sim, Val(:F_radiative_TOA)))

            if isempty(ccs.toa_net_source)
                radiation_sources_accum = sum(coupler_sim.fields.F_radiative_TOA .* coupler_sim.Δt_cpl) # ∫ J / m^2 dA
            else
                radiation_sources_accum =
                    sum(coupler_sim.fields.F_radiative_TOA .* coupler_sim.Δt_cpl) .+ ccs.toa_net_source[end] # ∫ J / m^2 dA
            end
            push!(ccs.toa_net_source, radiation_sources_accum)

            # save atmos
            previous = getproperty(ccs, sim_name)
            current = sum(Interfacer.get_field(sim, Val(:energy))) # # ∫ J / m^3 dV

            push!(previous, current)
            total += current + radiation_sources_accum

        elseif sim isa Interfacer.SurfaceModelSimulation || sim isa Interfacer.SurfaceStub
            # save surfaces
            area_fraction = Interfacer.get_field(sim, Val(:area_fraction))
            if isnothing(Interfacer.get_field(sim, Val(:energy)))
                previous = getproperty(ccs, sim_name)
                current = FT(0)
            else
                previous = getproperty(ccs, sim_name)
                current = sum(Interfacer.get_field(sim, Val(:energy)) .* area_fraction) # # ∫ J / m^3 dV
            end
            push!(previous, current)
            total += current
        end

    end
    push!(ccs.total, total)

    if runtime_check
        @assert abs((total[end] - total[1]) / total[end]) < 1e-4
    end
    return total

end

"""
    check_conservation!(
    cc::WaterConservationCheck,
    coupler_sim::Interfacer.CoupledSimulation,
    runtime_check = false,
    )

computes the total water, ∫ ρq_tot dV, of the various components
of the coupled simulations, and updates `cc` with the values.

Note: in the future this should not use `push!`.
"""
function check_conservation!(
    cc::WaterConservationCheck,
    coupler_sim::Interfacer.CoupledSimulation,
    runtime_check = false,
)

    ccs = cc.sums
    (; model_sims) = coupler_sim

    boundary_space = coupler_sim.boundary_space # thin shell approx (boundary_space[z=0] = boundary_space[z_top])

    total = 0

    # net precipitation (for surfaces that don't collect water)
    PE_net = coupler_sim.fields.P_net .+= swap_space!(zeros(boundary_space), surface_water_gain_from_rates(coupler_sim))

    # save surfaces
    for sim in model_sims
        sim_name = Symbol(Interfacer.name(sim))
        if sim isa Interfacer.AtmosModelSimulation

            # save atmos
            previous = getproperty(ccs, sim_name)
            current = sum(Interfacer.get_field(sim, Val(:water))) # kg (∫kg of water / m^3 dV)
            push!(previous, current)

        elseif sim isa Interfacer.SurfaceModelSimulation
            # save surfaces
            area_fraction = Interfacer.get_field(sim, Val(:area_fraction))
            if isnothing(Interfacer.get_field(sim, Val(:water)))
                previous = getproperty(ccs, sim_name)
                current = sum(PE_net .* area_fraction) # kg (∫kg of water / m^3 dV)
                push!(previous, current)
            else
                previous = getproperty(ccs, sim_name)
                current = sum(Interfacer.get_field(sim, Val(:water)) .* area_fraction) # kg (∫kg of water / m^3 dV)
                push!(previous, current)
            end
        end
        total += current


    end
    push!(ccs.total, total)

    if runtime_check
        @assert abs((total[end] - total[1]) / total[end]) < 1e-4
    end

    return total
end

"""
    surface_water_gain_from_rates(cs::Interfacer.CoupledSimulation)

Determines the total water content gain/loss of a surface from the begining of the simulation based on evaporation and precipitation rates.
"""
function surface_water_gain_from_rates(cs::Interfacer.CoupledSimulation)
    evaporation = cs.fields.F_turb_moisture # kg / m^2 / s / layer depth
    precipitation_l = cs.fields.P_liq
    precipitation_s = cs.fields.P_snow
    @. -(evaporation + precipitation_l + precipitation_s) * cs.Δt_cpl # kg / m^2 / layer depth
end

# setup the GKS socket application environment as nul for better performance and to avoid GKS connection errors while plotting
ENV["GKSwstype"] = "nul"

"""
    plot_global_conservation(
        cc::EnergyConservationCheck,
        coupler_sim::Interfacer.CoupledSimulation;
        figname1 = "total_energy.png",
        figname2 = "total_energy_log.png",
    )

Creates two plots of the globally integrated quantity (energy, ``\\rho e``):
1. global quantity of each model component as a function of time,
relative to the initial value;
2. fractional change in the sum of all components over time on a log scale.
"""
function plot_global_conservation(
    cc::AbstractConservationCheck,
    coupler_sim::Interfacer.CoupledSimulation;
    figname1 = "total.png",
    figname2 = "total_log.png",
)

    model_sims = coupler_sim.model_sims
    ccs = cc.sums

    days = collect(1:length(ccs[1])) * coupler_sim.Δt_cpl / 86400

    # evolution of energy of each component relative to initial value
    total = ccs.total  # total

    var_name = Interfacer.name(cc)
    Plots.plot(
        days,
        total .- total[1],
        label = "total",
        xlabel = "time [days]",
        ylabel = "$var_name: (t) - (t=0)",
        linewidth = 3,
    )
    for sim in model_sims
        sim_name = Interfacer.name(sim)
        global_field = getproperty(ccs, Symbol(sim_name))
        diff_global_field = (global_field .- global_field[1])
        Plots.plot!(days, diff_global_field[1:length(days)], label = sim_name)
    end
    if cc isa EnergyConservationCheck
        global_field = ccs.toa_net_source
        diff_global_field = (global_field .- global_field[1])
        Plots.plot!(days, diff_global_field[1:length(days)], label = "toa_net")
    end
    Plots.savefig(figname1)

    # evolution of log error of total
    Plots.plot(
        days,
        log.(abs.(total .- total[1]) / abs(total[1])),
        label = "tot",
        xlabel = "time [days]",
        ylabel = "log( | e(t) - e(t=0)| / e(t=0))",
    )
    Plots.savefig(figname2)
end

end # module
