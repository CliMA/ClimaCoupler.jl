T_S = ClimaCore.Fields.zeros(boundary_space) # temperature
z0m_S = ClimaCore.Fields.zeros(boundary_space)
z0b_S = ClimaCore.Fields.zeros(boundary_space)
ρ_sfc = ClimaCore.Fields.zeros(boundary_space)
q_sfc = ClimaCore.Fields.zeros(boundary_space)
albedo = ClimaCore.Fields.zeros(boundary_space)

F_A = ClimaCore.Fields.zeros(boundary_space) # aerodynamic turbulent fluxes
F_E = ClimaCore.Fields.zeros(boundary_space) # aerodynamic turbulent fluxes
F_R = ClimaCore.Fields.zeros(boundary_space) # radiative fluxes
P_liq = ClimaCore.Fields.zeros(boundary_space) # m/s precip
P_snow = ClimaCore.Fields.zeros(boundary_space) # m/s precip
"""
   atmos_push!(atmos_sim, boundary_space, F_A, F_R, F_E, P_liq, parsed_args)

updates F_A, F_R, P_liq, and F_E in place based on values used in the atmos_sim for the current step.
"""
function atmos_push!(atmos_sim, boundary_space, F_A, F_R, F_E, P_liq, P_snow,parsed_args)
    F_A .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_A, atmos_sim.integrator.p.dif_flux_energy)
    F_E .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_E, atmos_sim.integrator.p.dif_flux_ρq_tot)
    P_liq .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(P_liq, atmos_sim.integrator.p.col_integrated_rain)
    P_snow .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(P_liq, atmos_sim.integrator.p.col_integrated_snow)
    F_R .= ClimaCore.Fields.zeros(boundary_space)
    parsed_args["rad"] == "gray" ? dummmy_remap!(F_R, level(atmos_sim.integrator.p.ᶠradiation_flux, half)) : nothing
end

"""
   land_pull!(slab_sim::SlabSimulation, F_A, F_R, _...)

Updates the slab_sim cache state in place with the current values of F_A and F_R.
"""
function land_pull!(slab_sim::SlabSimulation, F_A, F_R, _...)
    @. slab_sim.integrator.p.F_aero = F_A
    @. slab_sim.integrator.p.F_rad = F_R
end

"""
   land_pull!(slab_sim::BucketSimulation, F_A, F_R, F_E, P_liq, ρ_sfc)

Updates the slab_sim cache state in place with the current values of F_A, F_R, F_E, P_liq, and ρ_sfc.

The surface air density is computed using the atmospheric state at the first level and making ideal gas
and hydrostatic balance assumptions. The land model does not compute the surface air density so this is
a reasonable stand-in.
"""
function land_pull!(slab_sim::BucketSimulation, F_A, F_R, F_E, P_liq, P_snow,ρ_sfc)
    FT = eltype(F_A)
    @. slab_sim.integrator.p.bucket.ρ_sfc = ρ_sfc
    @. slab_sim.integrator.p.bucket.SHF = F_A
    @. slab_sim.integrator.p.bucket.LHF = FT(0.0)
    ρ_liq = (LSMP.ρ_cloud_liq(slab_sim.params.earth_param_set))
    @. slab_sim.integrator.p.bucket.E = F_E / ρ_liq
    @. slab_sim.integrator.p.bucket.R_n = F_R
    @. slab_sim.integrator.p.bucket.P_liq = FT(-1.0) .* P_liq # land expects this to be positive
    @. slab_sim.integrator.p.bucket.P_snow = FT(-1.0) .* P_snow # land expects this to be positive
end

"""
   ocean_pull!(slab_ocean_sim, F_A, F_R)

Updates the slab_ocean_sim cache state in place with the current values of F_A and F_R.

The ocean model does not require moisture fluxes at the surface, so F_E is not returned.
"""
function ocean_pull!(slab_ocean_sim, F_A, F_R)
    @. slab_ocean_sim.integrator.p.F_aero = F_A
    @. slab_ocean_sim.integrator.p.F_rad = F_R
end

"""
   ice_pull!(slab_ice_sim, F_A, F_R)

Updates the slab_ice_sim cache state in place with the current values of F_A and F_R.

In the current version, the sea ice has a prescribed thickness, and we assume that it is not
sublimating. That contribution has been zeroed out in the atmos fluxes.
"""
function ice_pull!(slab_ice_sim, F_A, F_R)
    @. slab_ice_sim.integrator.p.Ya.F_aero = F_A
    @. slab_ice_sim.integrator.p.Ya.F_rad = F_R
end

"""
   atmos_pull!(atmos_sim,
               slab_ice_sim,
               slab_sim,
               slab_ocean_sim,
               boundary_space,
               prescribed_sst,
               z0m_S,
               z0b_S,
               T_S,
               ρ_sfc,
               q_sfc,
               ocean_params,
               SST,
               land_sea_mask
              )
Creates the surface fields for temperature, roughness length, albedo, and specific humidity; computes
turbulent surface fluxes; updates the atmosphere boundary flux cache variables in place; updates the
RRTMGP cache variables in place.
"""
function atmos_pull!(
    atmos_sim,
    slab_ice_sim,
    slab_sim,
    slab_ocean_sim,
    boundary_space,
    prescribed_sst,
    z0m_S,
    z0b_S,
    T_S,
    ρ_sfc,
    q_sfc,
    ocean_params,
    SST,
    land_sea_mask,
)
    thermo_params = CAP.thermodynamics_params(atmos_sim.integrator.p.params)

    univ_mask = parent(land_sea_mask) .- parent(slab_ice_sim.integrator.p.Ya.ice_mask .* FT(2))
    T_land = get_land_temp(slab_sim)
    z0m_land, z0b_land = get_land_roughness(slab_sim)
    combined_field = zeros(boundary_space)
    if prescribed_sst
        T_ocean = SST
        z0m_ocean = ocean_params.z0m
        z0b_ocean = ocean_params.z0b
        α_ocean = ocean_params.α
    else
        T_ocean = slab_ocean_sim.integrator.u.T_sfc
        z0m_ocean = slab_ocean_sim.integrator.p.params.z0m
        z0b_ocean = slab_ocean_sim.integrator.p.params.z0b
        α_ocean = slab_ocean_sim.integrator.p.params.α
    end

    parent(combined_field) .=
        combine_surface.(FT, univ_mask, parent(T_land), parent(T_ocean), parent(slab_ice_sim.integrator.u.T_sfc))
    dummmy_remap!(T_S, combined_field)

    parent(combined_field) .=
        combine_surface.(FT, univ_mask, z0m_land, z0m_ocean, slab_ice_sim.integrator.p.Ya.params.z0m)
    dummmy_remap!(z0m_S, combined_field)

    parent(combined_field) .=
        combine_surface.(FT, univ_mask, z0b_land, z0b_ocean, slab_ice_sim.integrator.p.Ya.params.z0b)
    dummmy_remap!(z0b_S, combined_field)

    set_ρ_sfc!(ρ_sfc, T_S, atmos_sim.integrator)
    # Now compute ocean and sea ice q_sat
    ocean_q_sfc = TD.q_vap_saturation_generic.(thermo_params, T_ocean, ρ_sfc, TD.Liquid())
    sea_ice_q_sfc = TD.q_vap_saturation_generic.(thermo_params, slab_ice_sim.integrator.u.T_sfc, ρ_sfc, TD.Ice())
    land_q_sfc = get_land_q(slab_sim, atmos_sim, T_land, ρ_sfc)
    # Pull q_sfc from land, and set q_sfc on surface with it and the above computed values.
    parent(combined_field) .=
        combine_surface.(FT, univ_mask, parent(land_q_sfc), parent(ocean_q_sfc), parent(sea_ice_q_sfc))
    dummmy_remap!(q_sfc, combined_field)
    α_land = land_albedo(slab_sim)

    α_ice = slab_ice_sim.params.α
    parent(combined_field) .= combine_surface.(FT, univ_mask, α_land, α_ocean, α_ice)
    dummmy_remap!(albedo, combined_field)
    atmos_sim.integrator.p.rrtmgp_model.diffuse_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(albedo), 1, length(parent(albedo)))
    atmos_sim.integrator.p.rrtmgp_model.direct_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(albedo), 1, length(parent(albedo)))

    # calculate turbulent fluxes on atmos grid and save in atmos cache
    info_sfc = (;
        T_sfc = T_S,
        ρ_sfc = ρ_sfc,
        q_sfc = q_sfc,
        z0m = z0m_S,
        z0b = z0b_S,
        ice_mask = slab_ice_sim.integrator.p.Ya.ice_mask,
    )
    calculate_surface_fluxes_atmos_grid!(atmos_sim.integrator, info_sfc)

    atmos_sim.integrator.p.rrtmgp_model.surface_temperature .= RRTMGPI.field2array(T_S)
end
