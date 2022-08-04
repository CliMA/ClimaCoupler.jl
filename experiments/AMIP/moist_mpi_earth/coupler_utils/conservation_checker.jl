# generalise this into coupler-specific function

abstract type AbstractCheck end

struct OnlineConservationCheck{A} <: AbstractCheck
    ρe_tot_atmos::A
    ρe_tot_land::A
    ρe_tot_ocean::A
    ρe_tot_seaice::A
    toa_net_source::A
    ice_base_source::A
end

function array2field(array, space) # reversing the RRTMGP function field2array (TODO: this now exists in ClimaAtmos)
    FT = eltype(array)
    Nq = Spaces.Quadratures.polynomial_degree(space.horizontal_space.quadrature_style) + 1
    ne = space.horizontal_space.topology.mesh.ne
    return Fields.Field(VIJFH{FT, Nq}(reshape(array, size(array, 1), Nq, Nq, 1, ne * ne * 6)), space)
end

"""
     check_conservation(
        cc::OnlineConservationCheck,
        coupler_sim,
        atmos_sim,
        land_sim = nothing,
        ocean_sim = nothing,
        seaice_sim = nothing,
        radiation = true,
        )
computes the total energy ∫ ρe dV of the various components
of the coupled simulations, and updates cc with the values.

Note: in the future this should not use ``push!``.
"""
function check_conservation(
    cc::OnlineConservationCheck,
    coupler_sim,
    atmos_sim,
    land_sim = nothing,
    ocean_sim = nothing,
    seaice_sim = nothing,
    radiation = true,
)
    @assert seaice_sim != nothing
    @assert atmos_sim != nothing
    face_space = axes(atmos_sim.integrator.u.f)
    FT = eltype(coupler_sim.land_mask)
    z = parent(Fields.coordinate_field(face_space).z)
    Δz_bot = FT(0.5) * (z[2, 1, 1, 1, 1] - z[1, 1, 1, 1, 1])

    u_atm = atmos_sim.integrator.u.c.ρe_tot
    T_land = land_sim !== nothing ? get_land_temp(land_sim) : nothing
    u_lnd = land_sim !== nothing ? swap_space!(T_land, coupler_sim.boundary_space) : nothing
    u_ocn = ocean_sim !== nothing ? swap_space!(ocean_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    u_ice = seaice_sim !== nothing ? swap_space!(seaice_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing

    # global sums
    atmos_e = sum(u_atm)

    if radiation

        Δz_top = FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1])
        n_faces = length(z[:, 1, 1, 1, 1])

        LWu_TOA = Fields.level(
            array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_lw_flux_up), face_space),
            n_faces - half,
        )
        SWd_TOA = Fields.level(
            array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_dn), face_space),
            n_faces - half,
        )
        SWu_TOA = Fields.level(
            array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_up), face_space),
            n_faces - half,
        )

        radiation_sources = -sum(SWd_TOA .- LWu_TOA .- SWu_TOA) ./ Δz_top
        radiation_sources_accum =
            size(cc.toa_net_source)[1] > 0 ? cc.toa_net_source[end] + radiation_sources .* coupler_sim.Δt_cpl :
            radiation_sources .* coupler_sim.Δt_cpl# accumulated radiation sources + sinks
        push!(cc.toa_net_source, radiation_sources_accum)
    end


    # Save atmos
    push!(cc.ρe_tot_atmos, atmos_e)

    # Surface masks
    univ_mask = parent(coupler_sim.land_mask) .- parent(seaice_sim.integrator.p.Ya.ice_mask .* FT(2))
    land_mask(u_lnd_1, univ_mask) = (univ_mask ≈ FT(1) ? u_lnd_1 : FT(0))
    ice_mask(u_ice_1, univ_mask) = (univ_mask ≈ FT(-2) ? u_ice_1 : FT(0))
    ocean_mask(u_ocn_1, univ_mask) = (univ_mask ≈ FT(0) ? u_ocn_1 : FT(0))


    

    # Save land
    e_lnd = land_sim !== nothing ?  get_land_energy(land_sim, coupler_sim.boundary_space) : nothing
    parent(e_lnd) .= land_mask.(parent(e_lnd), univ_mask)
    land_e = land_sim !== nothing ? sum(e_lnd) ./ Δz_bot : FT(0)
    push!(cc.ρe_tot_land, land_e)

    parent(u_ice) .= ice_mask.(parent(u_ice), univ_mask)
    seaice_e = seaice_sim !== nothing ? sum(get_slab_energy(seaice_sim, u_ice)) ./ Δz_bot : FT(0)
    push!(cc.ρe_tot_seaice, seaice_e)

    if ocean_sim != nothing
        parent(u_ocn) .= ocean_mask.(parent(u_ocn), univ_mask)
        ocean_e = sum(get_slab_energy(ocean_sim, u_ocn)) ./ Δz_bot
    else
        ocean_e = FT(0)
    end

    push!(cc.ρe_tot_ocean, ocean_e)


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
    diff_ρe_tot_slab_ocean = (cc.ρe_tot_ocean .- cc.ρe_tot_ocean[1])

    times_days = times ./ (24 * 60 * 60)
    Plots.plot(times_days, diff_ρe_tot_atmos[1:length(times_days)], label = "atmos")
    Plots.plot!(times_days, diff_ρe_tot_slab[1:length(times_days)], label = "land")
    Plots.plot!(times_days, diff_ρe_tot_slab_seaice[1:length(times_days)], label = "seaice")
    Plots.plot!(times_days, diff_toa_net_source[1:length(times_days)], label = "toa")
    Plots.plot!(times_days, diff_ρe_tot_slab_ocean[1:length(times_days)], label = "ocean")

    tot = cc.ρe_tot_atmos .+ cc.ρe_tot_ocean .+ cc.ρe_tot_land .+ cc.ρe_tot_seaice .+ cc.toa_net_source

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
