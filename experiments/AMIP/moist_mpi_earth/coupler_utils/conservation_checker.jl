# generalise this into coupler-specific function

abstract type AbstractCheck end

struct OnlineConservationCheck{A} <: AbstractCheck
    ρe_tot_atmos::A
    ρe_tot_land::A
    ρe_tot_ocean::A
    ρe_tot_seaice::A
    F_energy_land::A
    F_energy_ocean::A
    F_energy_ice::A
    F_rad_TOA::A
end
function check_conservation(
    cs::OnlineConservationCheck,
    coupler_sim,
    atmos_sim,
    land_sim,
    ocean_sim,
    seaice_sim,
    F_energy_surf,
    univ_mask, 
    radiation = true,
)
    #R_earth = 6.371229e6
    #R_top = 30e3+R_earth
    #A_top = 4π*R_top^2.0
    #A_bottom = 4π*R_earth^2.0
    
    z = parent(Fields.coordinate_field(face_space).z)
    Δz_bot = FT(0.5) * (z[2, 1, 1, 1, 1] - z[1, 1, 1, 1, 1])
    Δz_top = FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1])
    nz = length(z[:, 1, 1, 1, 1])

    u_atm = atmos_sim !== nothing ? atmos_sim.integrator.u.c.ρe : nothing
    u_lnd = land_sim !== nothing ? swap_space!(land_sim.integrator.u.bucket.T_sfc, coupler_sim.boundary_space) : nothing
    u_ocn = ocean_sim !== nothing ? swap_space!(ocean_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    u_ice = seaice_sim !== nothing ? swap_space!(seaice_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    # global sums
    if atmos_sim !== nothing
        atmos_e = sum(u_atm)

        if radiation

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
            radiation_sources_tsrs = radiation_sources .* atmos_sim.integrator.t # * coupler_sim.Δt # accumulate flux over coupling timestep (TODO: accumulate in atmos)
            push!(cs.F_rad_TOA, radiation_sources_tsrs)
            atmos_e = atmos_e .+ radiation_sources_tsrs # J 
        end
    end
    
#    mask = coupler_sim.mask
 #   univ_mask = parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2))
    ice_mask(u_ice_1, univ_mask) = (univ_mask ≈ FT(-2) ? u_ice_1 : FT(0))
    parent(u_ice) .= ice_mask.(parent(u_ice), univ_mask)

    if (prescribed_sst != true) 
        ocean_mask(u_ocn_1, univ_mask) = (univ_mask ≈ FT(0) ? u_ocn_1 : FT(0))
        parent(u_ocn) .= ocean_mask.(parent(u_ocn), univ_mask)
    end
    land_mask(u_lnd_1, univ_mask) = (univ_mask ≈ FT(1) ? u_lnd_1 : FT(0))
    parent(u_lnd) .= land_mask.(parent(u_lnd), univ_mask)

    f_land = swap_space!(F_energy_surf, coupler_sim.boundary_space)
    f_ocean = swap_space!(F_energy_surf, coupler_sim.boundary_space)
    f_ice = swap_space!(F_energy_surf, coupler_sim.boundary_space)
    
    parent(f_land) .= land_mask.(parent(F_energy_surf), univ_mask)
    parent(f_ocean) .= ocean_mask.(parent(F_energy_surf), univ_mask)
    parent(f_ice) .= ice_mask.(parent(F_energy_surf), univ_mask)
    
    land_e = land_sim !== nothing ? sum(get_bucket_energy(land_sim, u_lnd))./ Δz_bot : FT(0)
    ocean_e = ocean_sim !== nothing ? sum(get_slab_energy(ocean_sim, u_ocn))./ Δz_bot : FT(0)
    seaice_e = seaice_sim !== nothing ? sum(get_slab_energy(seaice_sim, u_ice))./ Δz_bot : FT(0)
    # save to coupler cache
    atmos_sim !== nothing ? push!(cs.ρe_tot_atmos, atmos_e) : nothing
    land_sim !== nothing ? push!(cs.ρe_tot_land, land_e) : nothing
    ocean_sim !== nothing ? push!(cs.ρe_tot_ocean, ocean_e) : nothing
    seaice_sim !== nothing ? push!(cs.ρe_tot_seaice, seaice_e) : nothing
    push!(cs.F_energy_land, sum(f_land)/Δz_bot)
    push!(cs.F_energy_ocean, sum(f_ocean)/Δz_bot)
    push!(cs.F_energy_ice, sum(f_ice)/Δz_bot)


end

