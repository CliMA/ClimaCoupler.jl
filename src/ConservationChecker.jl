"""
    ConservationChecker

This module contains functions that check global conservation of energy and water (and momentum - TODO).
"""
module ConservationChecker

import ..Interfacer, ..Utilities
import ClimaCore as CC

export AbstractConservationCheck, EnergyConservationCheck, WaterConservationCheck, check_conservation!

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

    FT = CC.Spaces.undertype(coupler_sim.boundary_space)
    total = FT(0)

    # save surfaces
    for sim in model_sims
        sim_name = Symbol(Interfacer.name(sim))
        if sim isa Interfacer.AtmosModelSimulation
            radiative_energy_flux_toa = Interfacer.get_field(sim, Val(:radiative_energy_flux_toa))

            # ClimaCore #1578, if radiative_energy_flux_toa comes from Fields.level
            # TODO: Fix this in ClimaCore
            rad_grid = CC.Spaces.grid(axes(radiative_energy_flux_toa))
            if rad_grid isa CC.Grids.LevelGrid
                if rad_grid.level isa CC.Utilities.PlusHalf
                    # FaceSpace
                    surface_integral =
                        sum(radiative_energy_flux_toa ./ CC.Fields.Δz_field(radiative_energy_flux_toa) .* 2)
                else
                    # Center
                    surface_integral = sum(radiative_energy_flux_toa ./ CC.Fields.Δz_field(radiative_energy_flux_toa))
                end
            else
                surface_integral = sum(radiative_energy_flux_toa)
            end

            if isempty(ccs.toa_net_source)
                radiation_sources_accum = surface_integral * FT(float(coupler_sim.Δt_cpl)) # ∫ J / m^2 dA
            else
                radiation_sources_accum = surface_integral * FT(float(coupler_sim.Δt_cpl)) + ccs.toa_net_source[end] # ∫ J / m^2 dA
            end
            push!(ccs.toa_net_source, radiation_sources_accum)

            # save atmos
            previous = getproperty(ccs, sim_name)
            current = sum(Interfacer.get_field(sim, Val(:energy))) # # ∫ J / m^3 dV

            push!(previous, current)
            total += current + radiation_sources_accum
        elseif sim isa Interfacer.SurfaceModelSimulation
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
        @assert abs((ccs.total[end] - ccs.total[1]) / ccs.total[end]) < 1e-4
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
    PE_net =
        coupler_sim.fields.P_net .+= Utilities.swap_space!(boundary_space, surface_water_gain_from_rates(coupler_sim))

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
        @assert abs((ccs.total[end] - ccs.total[1]) / ccs.total[end]) < 1e-4
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
    FT = eltype(evaporation)
    @. -(evaporation + precipitation_l + precipitation_s) * FT(float(cs.Δt_cpl)) # kg / m^2 / layer depth
end

end # module
