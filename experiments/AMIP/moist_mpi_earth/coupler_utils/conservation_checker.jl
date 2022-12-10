# generalise this into coupler-specific function

abstract type AbstractCheck end

struct OnlineConservationCheck{A} <: AbstractCheck
    ρe_tot_atmos::A
    ρe_tot_land::A
    ρe_tot_ocean::A
    ρe_tot_seaice::A
    toa_net_source::A
    ice_base_source::A
    friction_sink::A
end

"""
     check_conservation(
        cc::OnlineConservationCheck,
        coupler_sim,
        radiation = true,
        )

computes the total energy ∫ ρe dV of the various components
of the coupled simulations, and updates cc with the values.

Note: in the future this should not use ``push!``.
"""
function check_conservation(cc::OnlineConservationCheck, coupler_sim)
    @unpack model_sims, surface_masks = cs
    @unpack atmos_sim, land_sim, ocean_sim, ice_sim = model_sims
    radiation = integrator.p.radiation_model

    @assert ice_sim != nothing
    @assert atmos_sim != nothing

    FT = eltype(coupler_sim.surface_masks.land)

    u_atm = atmos_sim.integrator.u.c.ρe_tot

    if land_sim !== nothing
        e_per_area_land = zeros(axes(land_sim.integrator.u.bucket.W))
        get_land_energy(land_sim, e_per_area_land)
    end

    u_ocn = ocean_sim !== nothing ? swap_space!(ocean_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    u_ice = ice_sim !== nothing ? swap_space!(ice_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing

    # global sums
    atmos_e = sum(u_atm)

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

        radiation_sources = -sum(LWd_TOA .+ SWd_TOA .- LWu_TOA .- SWu_TOA) ./ Δz_top
        radiation_sources_accum =
            size(cc.toa_net_source)[1] > 0 ? cc.toa_net_source[end] + radiation_sources .* coupler_sim.Δt_cpl :
            radiation_sources .* coupler_sim.Δt_cpl# accumulated radiation sources + sinks
        push!(cc.toa_net_source, radiation_sources_accum)
    end

    # save atmos
    push!(cc.ρe_tot_atmos, atmos_e)

    # save land
    parent(e_per_area_land) .= parent(e_per_area_land .* surface_masks.land)
    land_e = land_sim !== nothing ? sum(e_per_area_land) : FT(0)
    push!(cc.ρe_tot_land, land_e)

    # save sea ice
    parent(u_ice) .= parent(u_ice .* surface_masks.ice)
    seaice_e = ice_sim !== nothing ? sum(get_slab_energy(ice_sim, u_ice)) : FT(0)
    push!(cc.ρe_tot_seaice, seaice_e)

    # save ocean
    if ocean_sim != nothing
        parent(u_ocn) .= parent(u_ocn .* surface_masks.ocean)
        ocean_e = sum(get_slab_energy(ocean_sim, u_ocn))
    else
        ocean_e = FT(0)
    end
    push!(cc.ρe_tot_ocean, ocean_e)

    # save surface friction sink
    push!(cc.friction_sink, sum((ke_dissipation(atmos_sim)) .* coupler_sim.Δt_cpl)) # ρ d ke_friction / dt 

end
function ke_dissipation(sim)
    drag_uv =
        .-Geometry.UVVector.(
            Geometry.Covariant12Vector.(
                sim.integrator.p.dif_flux_uₕ_bc.components.data.:1,
                sim.integrator.p.dif_flux_uₕ_bc.components.data.:2,
            )
        )
    dot.(Geometry.UVVector.(level(sim.integrator.u.c.uₕ, 1)), drag_uv) .* level(sim.integrator.u.c.ρ, 1)
end

import Plots
# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

"""
     plot_global_energy(cc,
                        coupler_sim,
                        figname1 = "total_energy.png",
                        figname2 = "total_energy_log.png")

Creates two plots: one showing fractional total energy change
 over time on a log scale,
and the other showing the energy of each component as a function of time,
relative to the initial value.
"""
function plot_global_energy(cc, coupler_sim, figname1 = "total_energy.png", figname2 = "total_energy_log.png")

    times = 0:(coupler_sim.Δt_cpl):(atmos_sim.integrator.t)
    diff_ρe_tot_atmos = (cc.ρe_tot_atmos .- cc.ρe_tot_atmos[1])
    diff_ρe_tot_slab = (cc.ρe_tot_land .- cc.ρe_tot_land[1])
    diff_ρe_tot_slab_seaice = (cc.ρe_tot_seaice .- cc.ρe_tot_seaice[1])
    diff_toa_net_source = (cc.toa_net_source .- cc.toa_net_source[1])
    # diff_friction_sink = (cc.friction_sink .- cc.friction_sink[1]) # currently, kinetic energy dissipation is neglected in total energy tendency
    diff_ρe_tot_slab_ocean = (cc.ρe_tot_ocean .- cc.ρe_tot_ocean[1])

    times_days = times ./ (24 * 60 * 60)
    Plots.plot(times_days, diff_ρe_tot_atmos[1:length(times_days)], label = "atmos")
    Plots.plot!(times_days, diff_ρe_tot_slab[1:length(times_days)], label = "land")
    Plots.plot!(times_days, diff_ρe_tot_slab_seaice[1:length(times_days)], label = "seaice")
    Plots.plot!(times_days, diff_toa_net_source[1:length(times_days)], label = "toa")
    # Plots.plot!(times_days, diff_friction_sink[1:length(times_days)], label = "friction")
    Plots.plot!(times_days, diff_ρe_tot_slab_ocean[1:length(times_days)], label = "ocean")

    tot = cc.ρe_tot_atmos .+ cc.ρe_tot_ocean .+ cc.ρe_tot_land .+ cc.ρe_tot_seaice .+ cc.toa_net_source #.+ cc.friction_sink

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

    @assert abs(tot[end] - tot[1]) < tot[end] * 1e-3 # TODO make this more stringent once small errors resolved
end
