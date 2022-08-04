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
function check_conservation(
    cs::OnlineConservationCheck,
    coupler_sim,
    atmos_sim = nothing,
    land_sim = nothing,
    ocean_sim = nothing,
    seaice_sim = nothing,
    radiation = true,
)
    mask = coupler_sim.land_mask

    z = parent(Fields.coordinate_field(face_space).z)
    Δz_bot = FT(0.5) * (z[2, 1, 1, 1, 1] - z[1, 1, 1, 1, 1])

    u_atm = atmos_sim !== nothing ? atmos_sim.integrator.u.c.ρe : nothing
    u_lnd = land_sim !== nothing ? swap_space!(land_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    u_ocn = ocean_sim !== nothing ? swap_space!(ocean_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    u_ice = seaice_sim !== nothing ? swap_space!(seaice_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing

    # global sums
    land_e = land_sim !== nothing ? sum(get_slab_energy(land_sim, u_lnd)) ./ Δz_bot : FT(0)

    if atmos_sim !== nothing
        atmos_e = sum(u_atm)

        if radiation

            Δz_top = FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1])
            nz = length(z[:, 1, 1, 1, 1])

            LWu_TOA = Fields.level(
                array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_lw_flux_up), face_space),
                nz - half,
            )
            SWd_TOA = Fields.level(
                array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_dn), face_space),
                nz - half,
            )
            SWu_TOA = Fields.level(
                array2field(FT.(atmos_sim.integrator.p.rrtmgp_model.face_sw_flux_up), face_space),
                nz - half,
            )

            radiation_sources = -sum(SWd_TOA .- LWu_TOA .- SWu_TOA) ./ Δz_top
            radiation_sources_accum =
                size(cs.toa_net_source)[1] > 0 ? cs.toa_net_source[end] + radiation_sources .* coupler_sim.Δt :
                radiation_sources .* coupler_sim.Δt# accumulated radiation sources + sinks


            push!(cs.toa_net_source, radiation_sources_accum)
        end
    end

    ocean_mask(univ_mask, value1 = FT(-0.5), value2 = FT(0.5)) =
        ((univ_mask >= FT(value1) && (univ_mask <= FT(value2))) ? FT(1) : FT(0))
    ice_mask_f(univ_mask, value1 = FT(-0.5)) = ((univ_mask <= FT(value1)) ? FT(1) : FT(0))

    if (prescribed_sst !== true) && (prescribed_sic == true)
        # zero out the ocean where there's sea ice
        univ_mask = parent(mask) .- parent(seaice_sim.integrator.p.Ya.ice_mask .* FT(2))
        value1, value2 = (FT(-0.5), FT(0.5))

        u_ocn_ = similar(u_ocn)
        parent(u_ocn_) .= ocean_mask.(parent(univ_mask)) .* parent(u_ocn)
        u_ocn = u_ocn_

        ocean_e = ocean_sim !== nothing ? sum(get_slab_energy(ocean_sim, u_ocn)) ./ Δz_bot : FT(0)
        seaice_e = seaice_sim !== nothing ? sum(get_slab_energy(seaice_sim, u_ice)) ./ Δz_bot : FT(0)

        # save to coupler cache
        atmos_sim !== nothing ? push!(cs.ρe_tot_atmos, atmos_e) : nothing
        land_sim !== nothing ? push!(cs.ρe_tot_land, land_e) : nothing
        ocean_sim !== nothing ? push!(cs.ρe_tot_ocean, ocean_e) : nothing
        seaice_sim !== nothing ? push!(cs.ρe_tot_seaice, seaice_e) : nothing
    end

    if (prescribed_sic == false)

        # ocean (from T_ml)
        u_ocn = seaice_sim !== nothing ? swap_space!(seaice_sim.integrator.u.T_ml, coupler_sim.boundary_space) : nothing
        parent(u_ocn) .= abs.(parent(mask) .- FT(1)) .* parent(u_ocn)
        ocean_e = ocean_sim !== nothing ? sum(get_ml_energy(seaice_sim, u_ocn)) ./ Δz_bot : FT(0)

        # ice (from h_ice)
        u_ice =
            seaice_sim !== nothing ? swap_space!(seaice_sim.integrator.u.h_ice, coupler_sim.boundary_space) : nothing
        parent(u_ice) .= abs.(parent(mask) .- FT(1)) .* parent(u_ice)
        seaice_e = ocean_sim !== nothing ? sum(get_dyn_ice_energy(seaice_sim, u_ice)) ./ Δz_bot : FT(0)

        # save to coupler cache
        atmos_sim !== nothing ? push!(cs.ρe_tot_atmos, atmos_e) : nothing
        land_sim !== nothing ? push!(cs.ρe_tot_land, land_e) : nothing
        ocean_sim !== nothing ? push!(cs.ρe_tot_ocean, ocean_e) : nothing
        seaice_sim !== nothing ? push!(cs.ρe_tot_seaice, seaice_e) : nothing
    end

end

function plot_global_energy(CS, coupler_sim, figname = "total_energy.png")
    times = (coupler_sim.Δt):(coupler_sim.Δt):(atmos_sim.integrator.t)
    diff_ρe_tot_atmos = (CS.ρe_tot_atmos .- CS.ρe_tot_atmos[1])
    diff_ρe_tot_slab = (CS.ρe_tot_land .- CS.ρe_tot_land[1])
    diff_ρe_tot_slab_seaice = (CS.ρe_tot_seaice .- CS.ρe_tot_seaice[1])
    diff_ρe_tot_slab_ocean = (CS.ρe_tot_ocean .- CS.ρe_tot_ocean[1])
    diff_toa_net_source = (CS.toa_net_source .- CS.toa_net_source[1])

    tot = CS.ρe_tot_atmos .+ CS.ρe_tot_ocean .+ CS.ρe_tot_land .+ CS.ρe_tot_seaice .+ CS.toa_net_source

    times_days = floor.(times ./ (24 * 60 * 60))
    Plots.plot(diff_ρe_tot_atmos, label = "atmos")
    Plots.plot!(diff_ρe_tot_slab, label = "land")
    Plots.plot!(diff_ρe_tot_slab_ocean, label = "ocean")
    Plots.plot!(diff_ρe_tot_slab_seaice, label = "seaice")
    Plots.plot!(diff_toa_net_source, label = "toa")

    Plots.plot!(
        tot .- tot[3],
        label = "tot",
        xlabel = "time [days]",
        ylabel = "energy(t) - energy(t=0) [J]",
        xticks = (collect(1:length(times))[1:1000:end], times_days[1:1000:end]),
        color = "black",
    )
    Plots.savefig(figname)
    Plots.plot(
        log.(abs.(tot .- tot[1]) / tot[1]),
        label = "tot",
        xlabel = "time [days]",
        ylabel = "log( | e(t) - e(t=0)| / e(t=0))",
        xticks = (collect(1:length(times))[1:50:end], times_days[1:50:end]),
    )
    Plots.savefig(figname[1:(end - 4)] * "_logerror.png")
end
