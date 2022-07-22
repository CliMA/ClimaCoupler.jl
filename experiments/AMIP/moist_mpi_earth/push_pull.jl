"""
   atmos_push!(cs)
updates F_A, F_R, P_liq, and F_E in place based on values used in the atmos_sim for the current step.
"""
function atmos_push!(cs)
    atmos_sim = cs.model_sims.atm
    csf = cs.fields
    dummmy_remap!(csf.F_A, atmos_sim.integrator.p.dif_flux_energy)
    dummmy_remap!(csf.F_E, atmos_sim.integrator.p.dif_flux_ρq_tot)
    dummmy_remap!(csf.P_liq, atmos_sim.integrator.p.col_integrated_rain)
    dummmy_remap!(csf.P_snow, atmos_sim.integrator.p.col_integrated_snow)
    cs.parsed_args["rad"] == "gray" ? dummmy_remap!(csf.F_R, level(atmos_sim.integrator.p.ᶠradiation_flux, half)) :
    nothing
end

"""
   land_pull!(cs)
Updates the slab_sim cache state in place with the current values of F_A, F_R, F_E, P_liq, P_snow, and ρ_sfc.
The surface air density is computed using the atmospheric state at the first level and making ideal gas
and hydrostatic balance assumptions. The land model does not compute the surface air density so this is
a reasonable stand-in.
"""
function land_pull!(cs)
    slab_sim = cs.model_sims.lnd
    csf = cs.fields
    FT = cs.FT
    @. slab_sim.integrator.p.bucket.ρ_sfc = csf.ρ_sfc
    @. slab_sim.integrator.p.bucket.SHF = csf.F_A
    @. slab_sim.integrator.p.bucket.LHF = FT(0.0)
    ρ_liq = (LSMP.ρ_cloud_liq(slab_sim.params.earth_param_set))
    @. slab_sim.integrator.p.bucket.E = csf.F_E / ρ_liq
    @. slab_sim.integrator.p.bucket.R_n = csf.F_R
    @. slab_sim.integrator.p.bucket.P_liq = FT(-1.0) .* csf.P_liq # land expects this to be positive
    @. slab_sim.integrator.p.bucket.P_snow = FT(-1.0) .* csf.P_snow # land expects this to be positive
end

"""
   ocean_pull!(cs)
Updates the ocean_sim cache state in place with the current values of F_A and F_R.
The ocean model does not require moisture fluxes at the surface, so F_E is not returned.
"""
function ocean_pull!(cs)
    ocean_sim = cs.model_sims.ocn
    csf = cs.fields
    @. ocean_sim.integrator.p.F_aero = csf.F_A
    @. ocean_sim.integrator.p.F_rad = csf.F_R
end

"""
   ice_pull!(cs)
Updates the ice_sim cache state in place with the current values of F_A and F_R.
In the current version, the sea ice has a prescribed thickness, and we assume that it is not
sublimating. That contribution has been zeroed out in the atmos fluxes.
"""
function ice_pull!(cs)
    ice_sim = cs.model_sims.ice
    csf = cs.fields
    @. ice_sim.integrator.p.Ya.F_aero = csf.F_A
    @. ice_sim.integrator.p.Ya.F_rad = csf.F_R
end

"""
   atmos_pull!(cs)
              
Creates the surface fields for temperature, roughness length, albedo, and specific humidity; computes
turbulent surface fluxes; updates the atmosphere boundary flux cache variables in place; updates the
RRTMGP cache variables in place.
"""
function atmos_pull!(cs)

    atmos_sim = cs.model_sims.atm
    slab_sim = cs.model_sims.lnd
    ocean_sim = cs.model_sims.ocn
    ice_sim = cs.model_sims.ice

    csf = cs.fields
    T_sfc_cpl = csf.T_S
    z0m_cpl = csf.z0m_S
    z0b_cpl = csf.z0b_S
    ρ_sfc_cpl = csf.ρ_sfc
    q_sfc_cpl = csf.q_sfc
    albedo_sfc_cpl = csf.albedo

    thermo_params = CAP.thermodynamics_params(atmos_sim.integrator.p.params)

    T_land = get_land_temp(slab_sim)
    z0m_land, z0b_land = get_land_roughness(slab_sim)
    T_ocean = ocean_sim.integrator.u.T_sfc
    z0m_ocean = ocean_sim.integrator.p.params.z0m
    z0b_ocean = ocean_sim.integrator.p.params.z0b
    α_ocean = ocean_sim.integrator.p.params.α
    T_ice = ice_sim.integrator.u.T_sfc
    ice_mask = ice_sim.integrator.p.Ya.ice_mask
    z0m_ice = ice_sim.integrator.p.params.z0m
    z0b_ice = ice_sim.integrator.p.params.z0b

    land_mask = cs.mask

    # mask of all models (1 = land, 0 = ocean, -2 = sea ice)
    univ_mask = parent(land_mask) .- parent(ice_mask .* FT(2))

    # combine models' surfaces onlo one coupler field 
    # surface temperature
    combined_field = zeros(boundary_space)
    parent(combined_field) .= combine_surface.(FT, univ_mask, parent(T_land), parent(T_ocean), parent(T_ice))
    dummmy_remap!(T_sfc_cpl, combined_field)

    # roughness length for momentum
    parent(combined_field) .= combine_surface.(FT, univ_mask, z0m_land, z0m_ocean, z0m_ice)
    dummmy_remap!(z0m_cpl, combined_field)

    # roughness length for tracers
    parent(combined_field) .= combine_surface.(FT, univ_mask, z0b_land, z0b_ocean, z0b_ice)
    dummmy_remap!(z0b_cpl, combined_field)

    # calculate atmospheric surface density 
    set_ρ_sfc!(ρ_sfc_cpl, T_sfc_cpl, atmos_sim.integrator)

    # surface specific humidity
    ocean_q_sfc = TD.q_vap_saturation_generic.(thermo_params, T_ocean, ρ_sfc_cpl, TD.Liquid())
    sea_ice_q_sfc = TD.q_vap_saturation_generic.(thermo_params, T_ice, ρ_sfc_cpl, TD.Ice())
    land_q_sfc = get_land_q(slab_sim, atmos_sim, T_land, ρ_sfc_cpl)
    parent(combined_field) .=
        combine_surface.(FT, univ_mask, parent(land_q_sfc), parent(ocean_q_sfc), parent(sea_ice_q_sfc))
    dummmy_remap!(q_sfc_cpl, combined_field)

    # albedo
    α_land = land_albedo(slab_sim)
    α_ice = ice_sim.integrator.p.params.α
    parent(combined_field) .= combine_surface.(FT, univ_mask, α_land, α_ocean, α_ice)
    dummmy_remap!(albedo_sfc_cpl, combined_field)
    atmos_sim.integrator.p.rrtmgp_model.diffuse_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(albedo_sfc_cpl), 1, length(parent(albedo_sfc_cpl)))
    atmos_sim.integrator.p.rrtmgp_model.direct_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(albedo_sfc_cpl), 1, length(parent(albedo_sfc_cpl)))

    # calculate turbulent fluxes on atmos grid and save in atmos cache
    info_sfc =
        (; T_sfc = T_sfc_cpl, ρ_sfc = ρ_sfc_cpl, q_sfc = q_sfc_cpl, z0m = z0m_cpl, z0b = z0b_cpl, ice_mask = ice_mask)
    calculate_surface_fluxes_atmos_grid!(atmos_sim.integrator, info_sfc)

    atmos_sim.integrator.p.rrtmgp_model.surface_temperature .= RRTMGPI.field2array(T_sfc_cpl)
end
