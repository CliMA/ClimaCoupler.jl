# init coupler's boundary fields for regridding (TODO: technically this can be bypassed by directly rigridding on model grids)
T_S = ClimaCore.Fields.zeros(boundary_space) # temperature
z0m_S = ClimaCore.Fields.zeros(boundary_space)
z0b_S = ClimaCore.Fields.zeros(boundary_space)
ρ_sfc = ClimaCore.Fields.zeros(boundary_space)
q_sfc = ClimaCore.Fields.zeros(boundary_space)
albedo = ClimaCore.Fields.zeros(boundary_space)

F_A = ClimaCore.Fields.zeros(boundary_space) # aerodynamic turbulent fluxes
F_E = ClimaCore.Fields.zeros(boundary_space) # evaporation due to turbulent fluxes
F_R = ClimaCore.Fields.zeros(boundary_space) # radiative fluxes
univ_mask = parent(ClimaCore.Fields.zeros(boundary_space))

# coupling methods
function atmos_push!(atmos_sim, boundary_space, F_A, F_E, F_R, parsed_args)
    F_A .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_A, atmos_sim.integrator.p.dif_flux_energy)
    F_E .= ClimaCore.Fields.zeros(boundary_space)
    dummmy_remap!(F_E, atmos_sim.integrator.p.dif_flux_ρq_tot)
    F_R .= ClimaCore.Fields.zeros(boundary_space)
    parsed_args["rad"] == "gray" ? dummmy_remap!(F_R, level(atmos_sim.integrator.p.ᶠradiation_flux, half)) : nothing
end

function bucket_pull!(bucket_sim, F_A, F_E, F_R, ρ_sfc)
    @. bucket_sim.integrator.p.bucket.ρ_sfc = ρ_sfc
    @. bucket_sim.integrator.p.bucket.SHF = F_A
    @. bucket_sim.integrator.p.bucket.LHF = FT(0.0)
    @. bucket_sim.integrator.p.bucket.E = F_E ./ ρ_cloud_liq(bucket_sim.params.earth_param_set)
    @. bucket_sim.integrator.p.bucket.R_n = F_R
end

# Ocean does not apply water flux boundary conditions, does not need F_E
function ocean_pull!(slab_ocean_sim, F_A, F_R)
    @. slab_ocean_sim.integrator.p.F_aero = F_A
    @. slab_ocean_sim.integrator.p.F_rad = F_R
end

# We assume ice is not sublimating, and have zeroed out that contributed to the atmos fluxes.
function ice_pull!(slab_ice_sim, F_A, F_R)
    @. slab_ice_sim.integrator.p.F_aero = F_A
    @. slab_ice_sim.integrator.p.F_rad = F_R
end


function atmos_pull!(coupler_sim, atmos_sim, slab_ice_sim, bucket_sim, slab_ocean_sim, boundary_space, prescribed_sst,  z0m_S,  z0b_S, T_S, ocean_params, SST, univ_mask)
    combined_field = zeros(boundary_space)
    if prescribed_sst == true
        ocean_p = ocean_params
        sst = SST
    else
        ocean_p = slab_ocean_sim.integrator.p.params
        sst = slab_ocean_sim.integrator.u.T_sfc
    end
    # This needs to be updated if the ice mask has changed in the past step
    mask = coupler_sim.mask
    univ_mask .= parent(mask) .- parent(slab_ice_sim.integrator.p.ice_mask .* FT(2))
    parent(combined_field) .=
        combine_surface.(
            univ_mask,
            parent(bucket_sim.integrator.u.bucket.T_sfc),
            parent(sst),
            parent(slab_ice_sim.integrator.u.T_sfc),
        ) # prescribed SSTs
    dummmy_remap!(T_S, combined_field)
    
    # Assume ice and ocean have the same roughness lengths
    parent(combined_field) .=
        combine_surface.(
            univ_mask,
            bucket_sim.params.z_0m,
            ocean_p.z0m
        )
    dummmy_remap!(z0m_S, combined_field)
    parent(combined_field) .=
        combine_surface.(
            univ_mask,
            bucket_sim.params.z_0b,
            ocean_p.z0b
        )
    dummmy_remap!(z0b_S, combined_field)
    
    # Compute ρ_sfc based atmos properties at lowest level
    set_ρ_sfc!(ρ_sfc, T_S, atmos_sim.integrator)
    # Now compute ocean and sea ice q_sat
    ocean_q_sfc = TD.q_vap_saturation_generic.(atmos_sim.integrator.p.params, sst, ρ_sfc, TD.Liquid())
    sea_ice_q_sfc = TD.q_vap_saturation_generic.(atmos_sim.integrator.p.params, slab_ice_sim.integrator.u.T_sfc, ρ_sfc, TD.Ice())
    # Pull q_sfc from land, and set q_sfc on surface with it and the above computed values.
    parent(combined_field) .=
        combine_surface.(
            univ_mask,
            parent(bucket_sim.integrator.p.bucket.q_sfc),
            parent(ocean_q_sfc),
            parent(sea_ice_q_sfc),
        )
    dummmy_remap!(q_sfc, combined_field)

    atmos_sim.integrator.p.rrtmgp_model.surface_temperature .= field2array(T_S) # supplied to atmos for radiation

    coords = ClimaCore.Fields.coordinate_field(axes(bucket_sim.integrator.u.bucket.S))
    α_land = surface_albedo.(Ref(bucket_sim.params.albedo), coords, bucket_sim.integrator.u.bucket.S, bucket_sim.params.S_c)
    parent(combined_field) .=
        combine_surface.(
            univ_mask,
            parent(α_land),
            ocean_p.α,
            slab_ice_sim.params.α
        )
    dummmy_remap!(albedo, combined_field)
    atmos_sim.integrator.p.rrtmgp_model.diffuse_sw_surface_albedo .=reshape(field2array(albedo), 1, length(parent(albedo)))
    atmos_sim.integrator.p.rrtmgp_model.direct_sw_surface_albedo .=reshape(field2array(albedo), 1, length(parent(albedo)))
    
    # calculate turbulent fluxes on atmos grid and save in atmos cache
    info_sfc = (; T_sfc = T_S, ρ_sfc = ρ_sfc, q_sfc = q_sfc, z0m = z0m_S, z0b = z0b_S, ice_mask = slab_ice_sim.integrator.p.ice_mask)
    calculate_surface_fluxes_atmos_grid!(atmos_sim.integrator, info_sfc)
end
