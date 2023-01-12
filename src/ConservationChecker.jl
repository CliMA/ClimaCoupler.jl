"""
    ConservationChecker

This module contains functions that check global conservation of energy and water (and momentum - TODO).
"""
module ConservationChecker

using ClimaCore: ClimaCore, Geometry, Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaCore.Utilities: half
using ClimaComms
using NCDatasets
using ClimaCoreTempestRemap
using Dates
using JLD2
using UnPack
using Plots
using ClimaAtmos: RRTMGPI
using ClimaLSM
using ..Utilities: CoupledSimulation, swap_space!

export AbstractCheck, EnergyConservationCheck, WaterConservationCheck, check_conservation!, plot_global_conservation

abstract type AbstractCheck end

"""
    EnergyConservationCheck{A} <: AbstractCheck

Struct of type `AbstractCheck` containing global energy conservation logs.
"""
struct EnergyConservationCheck{A} <: AbstractCheck
    ρe_tot_atmos::A
    ρe_tot_land::A
    ρe_tot_ocean::A
    ρe_tot_seaice::A
    toa_net_source::A
    ice_base_source::A
end

"""
    WaterConservationCheck{A} <: AbstractCheck

Struct of type `AbstractCheck` containing global water mass conservation logs.
"""
struct WaterConservationCheck{A} <: AbstractCheck
    ρq_tot_atmos::A
    ρq_tot_ocean::A
    ρq_tot_land::A
    ρq_tot_seaice::A
end

"""
    check_conservation!(coupler_sim::CoupledSimulation, get_slab_energy, get_land_energy)

itertes over all specified conservation checks.
"""
check_conservation!(coupler_sim::CoupledSimulation, get_slab_energy, get_land_energy) =
    map(x -> check_conservation!(x, coupler_sim, get_slab_energy, get_land_energy), coupler_sim.conservation_checks)

"""
        check_conservation!(
        cc::EnergyConservationCheck,
        coupler_sim,
        get_slab_energy,
        get_land_energy,
        )

computes the total energy, ∫ ρe dV, of the various components
of the coupled simulations, and updates `cc` with the values.

TODO: move `get_slab_energy` and `get_land_energy` to their respective sims upon optimization refactor.
"""
function check_conservation!(
    cc::EnergyConservationCheck,
    coupler_sim::CoupledSimulation,
    get_slab_energy,
    get_land_energy,
)
    @unpack model_sims, surface_masks = coupler_sim
    @unpack atmos_sim, land_sim, ocean_sim, ice_sim = model_sims
    radiation = atmos_sim.integrator.p.radiation_model # TODO: take out of global scope in ClimaAtmos
    boundary_space = coupler_sim.boundary_space # thin shell approx (boundary_space[z=0] = boundary_space[z_top])

    FT = eltype(coupler_sim.surface_masks.land)

    # save radiation source
    if radiation != nothing
        face_space = axes(atmos_sim.integrator.u.f)
        z = parent(Fields.coordinate_field(face_space).z)
        Δz_top = round(FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1]))
        n_faces = length(z[:, 1, 1, 1, 1])

        LWd_TOA = Fields.level(
            RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_dn), face_space),
            n_faces - half,
        )
        LWu_TOA = Fields.level(
            RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_up), face_space),
            n_faces - half,
        )
        SWd_TOA = Fields.level(
            RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_dn), face_space),
            n_faces - half,
        )
        SWu_TOA = Fields.level(
            RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_up), face_space),
            n_faces - half,
        )

        coupler_sim.fields.F_R_TOA .-=
            swap_space!(LWd_TOA .+ SWd_TOA .- LWu_TOA .- SWu_TOA, boundary_space) .* coupler_sim.Δt_cpl
        radiation_sources_accum = sum(coupler_sim.fields.F_R_TOA) # accumulated radiation sources + sinks [J]
        push!(cc.toa_net_source, radiation_sources_accum)
    else
        push!(cc.toa_net_source, FT(0))
    end

    # save atmos
    push!(cc.ρe_tot_atmos, sum(atmos_sim.integrator.u.c.ρe_tot))

    # save land
    if land_sim !== nothing
        e_per_area_land = zeros(axes(land_sim.integrator.u.bucket.W))
        get_land_energy(land_sim, e_per_area_land)
        push!(cc.ρe_tot_land, sum(e_per_area_land .* surface_masks.land))
    else
        push!(cc.ρe_tot_land, FT(0))
    end

    # save sea ice
    if ice_sim !== nothing
        push!(cc.ρe_tot_seaice, sum(get_slab_energy(ice_sim, ice_sim.integrator.u.T_sfc) .* surface_masks.ice))
    else
        push!(cc.ρe_tot_seaice, FT(0))
    end

    # save ocean
    if ocean_sim !== nothing
        push!(cc.ρe_tot_ocean, sum(get_slab_energy(ocean_sim, ocean_sim.integrator.u.T_sfc) .* surface_masks.ocean))
    else
        push!(cc.ρe_tot_ocean, FT(0))
    end

    push!(cc.ice_base_source, FT(0))

    tot = @. cc.ρe_tot_atmos + cc.ρe_tot_ocean + cc.ρe_tot_land + cc.ρe_tot_seaice + cc.toa_net_source
    # @assert abs((tot[end] - tot[1]) / tot[end]) < 1e-4
end

"""
    check_conservation!(
    cc::WaterConservationCheck,
    coupler_sim,
    get_slab_energy,
    get_land_energy,
    )

computes the total water, ∫ ρq_tot dV, of the various components
of the coupled simulations, and updates `cc` with the values.

Note: in the future this should not use `push!`.
"""
function check_conservation!(
    cc::WaterConservationCheck,
    coupler_sim::CoupledSimulation,
    get_slab_energy,
    get_land_energy,
)
    @unpack model_sims, surface_masks = coupler_sim
    @unpack atmos_sim, land_sim, ocean_sim, ice_sim = model_sims

    boundary_space = coupler_sim.boundary_space
    FT = eltype(coupler_sim.surface_masks.land)

    # save atmos
    push!(cc.ρq_tot_atmos, sum(atmos_sim.integrator.u.c.ρq_tot))

    # save land
    if land_sim !== nothing
        ρ_cloud_liq = ClimaLSM.LSMP.ρ_cloud_liq(land_sim.params.earth_param_set)
        water_content =
            @. (land_sim.integrator.u.bucket.σS + land_sim.integrator.u.bucket.W + land_sim.integrator.u.bucket.Ws) # m^3 water / land area / layer height
        parent(water_content) .= parent(water_content .* surface_masks.land) * ρ_cloud_liq # kg / land area / layer height
        push!(cc.ρq_tot_land, sum(water_content)) # kg (∫ water_content dV)
    else
        push!(cc.ρq_tot_land, FT(0))
    end

    # save sea ice
    coupler_sim.fields.P_net .-= swap_space!(surface_water_gain_from_rates(coupler_sim), boundary_space) # accumulated surface water gain
    if ice_sim !== nothing
        push!(cc.ρq_tot_seaice, sum(coupler_sim.fields.P_net .* surface_masks.ice)) # kg (∫ P_net dV)
    else
        push!(cc.ρq_tot_seaice, FT(0))
    end

    # save ocean
    if ocean_sim !== nothing
        push!(cc.ρq_tot_ocean, sum(coupler_sim.fields.P_net .* surface_masks.ocean))  # kg (∫ P_net dV)
    else
        push!(cc.ρq_tot_ocean, FT(0))
    end

    tot = @. cc.ρq_tot_atmos + cc.ρq_tot_land + cc.ρq_tot_seaice + cc.ρq_tot_ocean
    # @assert abs((tot[end] - tot[1]) / tot[end]) < 1e-2
end

"""
    surface_water_gain_from_rates(cs)

Determines the total water content gain/loss of a surface from the begining of the simulation based on evaporation and precipitation rates.
"""
function surface_water_gain_from_rates(cs::CoupledSimulation)
    evaporation = cs.fields.F_E # kg / m^2 / s / layer depth
    precipitation_l = cs.fields.P_liq
    precipitation_s = cs.fields.P_snow
    @. (evaporation + precipitation_l + precipitation_s) * cs.Δt_cpl # kg / m^2 / layer depth
end

# setup the GKS socket application environment as nul for better performance and to avoid GKS connection errors while plotting
ENV["GKSwstype"] = "nul"

"""
    plot_global_conservation(
        cc::EnergyConservationCheck,
        coupler_sim::CoupledSimulation;
        figname1 = "total_energy.png",
        figname2 = "total_energy_log.png",
    )

Creates two plots of the globally integrated quantity (energy, ``\\rho e``):
1. global quantity of each model component as a function of time,
relative to the initial value;
2. fractional change in the sum of all components over time on a log scale.
"""
function plot_global_conservation(
    cc::EnergyConservationCheck,
    coupler_sim::CoupledSimulation;
    figname1 = "total_energy.png",
    figname2 = "total_energy_log.png",
)

    times = collect(1:length(cc.ρe_tot_atmos)) * coupler_sim.Δt_cpl
    diff_ρe_tot_atmos = (cc.ρe_tot_atmos .- cc.ρe_tot_atmos[1])
    diff_ρe_tot_slab = (cc.ρe_tot_land .- cc.ρe_tot_land[1])
    diff_ρe_tot_slab_seaice = (cc.ρe_tot_seaice .- cc.ρe_tot_seaice[1])
    diff_toa_net_source = (cc.toa_net_source .- cc.toa_net_source[1])
    diff_ρe_tot_slab_ocean = (cc.ρe_tot_ocean .- cc.ρe_tot_ocean[1])

    times_days = times ./ (24 * 60 * 60)
    Plots.plot(times_days, diff_ρe_tot_atmos[1:length(times_days)], label = "atmos")
    Plots.plot!(times_days, diff_ρe_tot_slab[1:length(times_days)], label = "land")
    Plots.plot!(times_days, diff_ρe_tot_slab_seaice[1:length(times_days)], label = "seaice")
    Plots.plot!(times_days, diff_toa_net_source[1:length(times_days)], label = "toa")
    Plots.plot!(times_days, diff_ρe_tot_slab_ocean[1:length(times_days)], label = "ocean")

    tot = @. cc.ρe_tot_atmos + cc.ρe_tot_ocean + cc.ρe_tot_land + cc.ρe_tot_seaice + cc.toa_net_source

    Plots.plot!(
        times_days,
        tot .- tot[1],
        label = "tot",
        xlabel = "time [days]",
        ylabel = "energy(t) - energy(t=0) [J]",
    )
    Plots.savefig(figname1)
    Plots.plot(
        times_days,
        log.(abs.(tot .- tot[1]) / tot[1]),
        label = "tot",
        xlabel = "time [days]",
        ylabel = "log( | e(t) - e(t=0)| / e(t=0))",
    )
    Plots.savefig(figname2)
end

"""
    plot_global_conservation(
        cc::WaterConservationCheck,
        coupler_sim::CoupledSimulation;
        figname1 = "total_energy.png",
        figname2 = "total_energy_log.png",
    )

Creates two plots of the globally integrated quantity (water, ``\\rho q_{tot}``):
1. global quantity of each model component as a function of time,
relative to the initial value;
2. fractional change in the sum of all components over time on a log scale.
"""
function plot_global_conservation(
    cc::WaterConservationCheck,
    coupler_sim::CoupledSimulation{FT};
    figname1 = "total_water.png",
    figname2 = "total_water_log.png",
) where {FT}
    times = collect(1:length(cc.ρq_tot_atmos)) * coupler_sim.Δt_cpl
    diff_ρe_tot_atmos = (cc.ρq_tot_atmos .- cc.ρq_tot_atmos[1])
    diff_ρe_tot_slab = (cc.ρq_tot_land .- cc.ρq_tot_land[1]) * FT(1e3)
    diff_ρe_tot_slab_seaice = (cc.ρq_tot_seaice .- cc.ρq_tot_seaice[1])
    diff_ρe_tot_slab_ocean = (cc.ρq_tot_ocean .- cc.ρq_tot_ocean[1])

    times_days = times ./ (24 * 60 * 60)
    Plots.plot(times_days, diff_ρe_tot_atmos[1:length(times_days)], label = "atmos")
    Plots.plot!(times_days, diff_ρe_tot_slab[1:length(times_days)], label = "land")
    Plots.plot!(times_days, diff_ρe_tot_slab_seaice[1:length(times_days)], label = "seaice")
    Plots.plot!(times_days, diff_ρe_tot_slab_ocean[1:length(times_days)], label = "ocean")

    tot = @. cc.ρq_tot_atmos + cc.ρq_tot_land * FT(1e3) + cc.ρq_tot_seaice + cc.ρq_tot_ocean

    Plots.plot!(times_days, tot .- tot[1], label = "tot", xlabel = "time [days]", ylabel = "water(t) - water(t=0) [kg]")
    Plots.savefig(figname1)
    Plots.plot(
        times_days,
        log.(abs.(tot .- tot[1]) / tot[1]),
        label = "tot",
        xlabel = "time [days]",
        ylabel = "log( | w(t) - w(t=0)| / w(t=0))",
    )
    Plots.savefig(figname2)
end

end # module
