"""
   atmos_push!(cs)

updates F_A, F_R, P_liq, and F_E in place based on values used in the atmos_sim for the current step.
"""
function atmos_push!(cs)
    atmos_sim = cs.model_sims.atmos_sim
    csf = cs.fields
    dummmy_remap!(csf.F_A, .-atmos_sim.integrator.p.dif_flux_energy_bc)
    dummmy_remap!(csf.F_E, .-atmos_sim.integrator.p.dif_flux_ρq_tot_bc)
    dummmy_remap!(csf.F_R, level(atmos_sim.integrator.p.ᶠradiation_flux, half))
    dummmy_remap!(csf.P_liq, atmos_sim.integrator.p.col_integrated_rain .+ atmos_sim.integrator.p.col_integrated_snow)
end

"""
   land_pull!(cs)

Updates the land_sim cache state in place with the current values of F_A, F_R, F_E, P_liq, and ρ_sfc.
The surface air density is computed using the atmospheric state at the first level and making ideal gas
and hydrostatic balance assumptions. The land model does not compute the surface air density so this is
a reasonable stand-in.
"""
function land_pull!(cs)
    FT = float_type(cs)
    land_sim = cs.model_sims.land_sim
    csf = cs.fields
    land_mask = cs.surface_masks.land
    parent(land_sim.integrator.p.bucket.ρ_sfc) .= parent(csf.ρ_sfc)
    parent(land_sim.integrator.p.bucket.turbulent_energy_flux) .=
        apply_mask.(parent(land_mask), >, parent(csf.F_A), parent(csf.F_A) .* FT(0), FT(0))
    ρ_liq = (LSMP.ρ_cloud_liq(land_sim.model.parameters.earth_param_set))
    parent(land_sim.integrator.p.bucket.evaporation) .=
        apply_mask.(parent(land_mask), >, parent(csf.F_E) ./ ρ_liq, parent(csf.F_E) .* FT(0), FT(0))
    parent(land_sim.integrator.p.bucket.R_n) .=
        apply_mask.(parent(land_mask), >, parent(csf.F_R), parent(csf.F_R) .* FT(0), FT(0))
    parent(land_sim.integrator.p.bucket.P_liq) .= FT(-1.0) .* parent(csf.P_liq) # land expects this to be positive
    parent(land_sim.integrator.p.bucket.P_snow) .= FT(0.0) .* parent(csf.P_snow)

end

"""
   ocean_pull!(cs)

Updates the ocean_sim cache state in place with the current values of F_A and F_R.
The ocean model does not require moisture fluxes at the surface, so F_E is not returned.
"""
function ocean_pull!(cs)
    ocean_sim = cs.model_sims.ocean_sim
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
    FT = float_type(cs)
    ice_sim = cs.model_sims.ice_sim
    csf = cs.fields
    ice_mask = cs.surface_masks.ice
    parent(ice_sim.integrator.p.F_rad) .=
        apply_mask.(parent(ice_mask), >, parent(csf.F_R), parent(csf.F_R) .* FT(0), FT(0))
    parent(ice_sim.integrator.p.F_aero) .=
        apply_mask.(parent(ice_mask), >, parent(csf.F_A), parent(csf.F_A) .* FT(0), FT(0))
end

"""
   atmos_pull!(cs)

Creates the surface fields for temperature, roughness length, albedo, and specific humidity; computes
turbulent surface fluxes; updates the atmosphere boundary flux cache variables in place; updates the
RRTMGP cache variables in place.
"""
function atmos_pull!(cs)

    @unpack model_sims = cs
    @unpack atmos_sim, land_sim, ocean_sim, ice_sim = model_sims
    radiation = atmos_sim.integrator.p.radiation_model

    csf = cs.fields
    T_sfc_cpl = csf.T_S
    z0m_cpl = csf.z0m_S
    z0b_cpl = csf.z0b_S
    albedo_sfc_cpl = csf.albedo
    beta_sfc_cpl = csf.beta

    thermo_params = CAP.thermodynamics_params(atmos_sim.integrator.p.params)

    T_land = get_land_temp(land_sim)
    z0m_land, z0b_land = get_land_roughness(land_sim)
    T_ocean = ocean_sim.integrator.u.T_sfc
    z0m_ocean = ocean_sim.integrator.p.params.z0m
    z0b_ocean = ocean_sim.integrator.p.params.z0b
    α_ocean = ocean_sim.integrator.p.params.α
    T_ice = ice_sim.integrator.u.T_sfc
    ice_mask = ice_sim.integrator.p.ice_mask
    z0m_ice = ice_sim.integrator.p.params.z0m
    z0b_ice = ice_sim.integrator.p.params.z0b

    update_masks!(cs)

    # combine models' surfaces onlo one coupler field
    combined_field = zeros(boundary_space)

    # surface temperature
    combine_surfaces!(combined_field, cs.surface_masks, (; land = T_land, ocean = T_ocean, ice = T_ice))
    dummmy_remap!(T_sfc_cpl, combined_field)

    # roughness length for momentum
    combine_surfaces!(combined_field, cs.surface_masks, (; land = z0m_land, ocean = z0m_ocean, ice = z0m_ice))
    dummmy_remap!(z0m_cpl, combined_field)

    # roughness length for tracers
    combine_surfaces!(combined_field, cs.surface_masks, (; land = z0b_land, ocean = z0b_ocean, ice = z0b_ice))
    dummmy_remap!(z0b_cpl, combined_field)

    # TODO replace functions that used ρ_sfc and q_sfc from atmos

    # albedo
    α_land = land_albedo(land_sim)

    α_ice = ice_sim.integrator.p.params.α
    combine_surfaces!(combined_field, cs.surface_masks, (; land = α_land, ocean = α_ocean, ice = α_ice))
    dummmy_remap!(albedo_sfc_cpl, combined_field)

    if radiation != nothing
        atmos_sim.integrator.p.radiation_model.diffuse_sw_surface_albedo .=
            reshape(RRTMGPI.field2array(albedo_sfc_cpl), 1, length(parent(albedo_sfc_cpl)))
        atmos_sim.integrator.p.radiation_model.direct_sw_surface_albedo .=
            reshape(RRTMGPI.field2array(albedo_sfc_cpl), 1, length(parent(albedo_sfc_cpl)))
        atmos_sim.integrator.p.radiation_model.surface_temperature .= RRTMGPI.field2array(T_sfc_cpl)
    end

    # beta factor for evaporation - should be 1 over ocean and ice
    β_ocean = ones(axes(combined_field))
    β_ice = ones(axes(combined_field))
    β_land = land_beta(land_sim)

    combine_surfaces!(combined_field, cs.surface_masks, (; land = β_land, ocean = β_ocean, ice = β_ice))
    dummmy_remap!(beta_sfc_cpl, combined_field)
    @assert all(i -> parent(beta_sfc_cpl)[i] == 1, findall(==(1), parent(ice_mask)))

    # calculate turbulent fluxes on atmos grid and save in atmos cache
    parent(atmos_sim.integrator.p.T_sfc) .= parent(T_sfc_cpl)
    if :z0b in propertynames(atmos_sim.integrator.p.surface_scheme)
        parent(atmos_sim.integrator.p.sfc_inputs.z0b) .= parent(z0b_cpl)
        parent(atmos_sim.integrator.p.sfc_inputs.z0m) .= parent(z0m_cpl)
    end

    Fields.bycolumn(axes(atmos_sim.integrator.p.ts_sfc)) do colidx
        ClimaAtmos.set_surface_thermo_state!(
            ClimaAtmos.Decoupled(),
            atmos_sim.integrator.p.surface_scheme.sfc_thermo_state_type,
            atmos_sim.integrator.p.ts_sfc[colidx],
            atmos_sim.integrator.p.T_sfc[colidx],
            Spaces.level(atmos_sim.integrator.p.ᶜts[colidx], 1),
            thermo_params,
            atmos_sim.integrator.t,
        )

        get_surface_fluxes!(
            atmos_sim.integrator.u,
            atmos_sim.integrator.p,
            atmos_sim.integrator.t,
            colidx,
            atmos_sim.integrator.p.atmos.vert_diff,
        )
    end

    # TODO correct for ice coverage
end

function atmos_pull!(cs, surfces)
    # placehoolder: add method to calculate fluxes above individual surfaces and then split fluxes (separate PR)
end

maximumfield(el1, el2) = maximum([el1, el2])
