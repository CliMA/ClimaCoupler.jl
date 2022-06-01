using ClimaCore.Geometry: ⊗
using ClimaCore.Utilities: half, PlusHalf

"""
calculate_surface_fluxes_atmos_grid!(integrator)

- TODO: generalize interface for regridding and take land state out of atmos's integrator.p
"""
function calculate_surface_fluxes_atmos_grid!(integrator, info_sfc)
    p = integrator.p
    (; ᶜts, dif_flux_energy, dif_flux_ρq_tot, dif_flux_uₕ, ∂F_aero∂T_sfc, params, Cd, Ch) = p

    (; T_sfc, z0m, z0b, ice_mask) = info_sfc
    Y = integrator.u

    # Turbulent surface flux calculation
    tsf =
        constant_T_saturated_surface_coefs_coupled.(
            Spaces.level(ᶜts, 1),
            Geometry.UVVector.(Spaces.level(Y.c.uₕ, 1)),
            Spaces.level(Fields.coordinate_field(Y.c).z, 1),
            FT(0), # TODO: get actual value of z_sfc
            swap_space!(T_sfc, axes(Spaces.level(Y.c, 1))), # remove when same instance issue is resolved
            params,
            swap_space!(z0m, axes(Spaces.level(Y.c, 1))), # TODO: get these roughness lengths from land
            swap_space!(z0b, axes(Spaces.level(Y.c, 1))),
            Cd,
            Ch,
        )

    # Total energy flux
    if :ρe in propertynames(Y.c)

        flux_energy = ones(axes(dif_flux_energy))
        parent(flux_energy) .= parent(tsf.shf .+ tsf.lhf .* swap_space!(abs.(ice_mask .- FT(1)), axes(tsf.shf)))  # only SHF above sea ice
        @. dif_flux_energy = Geometry.WVector(flux_energy) #Geometry.WVector.(swap_space!(tsf.shf .+ tsf.lhf, axes(dif_flux_energy)) )
    end

    # Moisture mass flux 
    if :ρq_tot in propertynames(Y.c)
        flux_mass = ones(axes(dif_flux_ρq_tot))
        parent(flux_mass) .= parent(tsf.E .* swap_space!(abs.(ice_mask .- FT(1)), axes(tsf.E)))
        @. dif_flux_ρq_tot = Geometry.WVector(flux_mass) # no E above sea ice
    end

    # Momentum flux
    u_space = axes(tsf.ρτxz) # TODO: delete when "space not the same instance" error is dealt with 
    normal = Geometry.WVector.(ones(u_space)) # TODO: this will need to change for topography
    ρ_1 = Fields.Field(Fields.field_values(Fields.level(Y.c.ρ, 1)), u_space) # TODO: delete when "space not the same instance" error is dealt with 
    if :uₕ in propertynames(Y.c)
        parent(dif_flux_uₕ) .=  # TODO: remove parent when "space not the same instance" error is dealt with 
            parent(
                Geometry.Contravariant3Vector.(normal) .⊗
                Geometry.Covariant12Vector.(Geometry.UVVector.(tsf.ρτxz ./ ρ_1, tsf.ρτyz ./ ρ_1)),
            )
    end

    # calculate gradient - TODO: make this just optional (only required by stub sea ice model)
    ΔT_sfc = FT(0.1) # following FMS
    tsf1 =
        constant_T_saturated_surface_coefs_coupled.(
            Spaces.level(ᶜts, 1),
            Geometry.UVVector.(Spaces.level(Y.c.uₕ, 1)),
            Spaces.level(Fields.coordinate_field(Y.c).z, 1),
            FT(0), # TODO: get actual value of z_sfc
            swap_space!(T_sfc .+ ΔT_sfc, axes(Spaces.level(Y.c, 1))), # remove when same instance issue is resolved
            params,
            swap_space!(z0m, axes(Spaces.level(Y.c, 1))), # TODO: get these roughness lengths from land
            swap_space!(z0b, axes(Spaces.level(Y.c, 1))),
            Cd,
            Ch,
        )

    p.∂F_aero∂T_sfc .= ((tsf1.shf .+ tsf1.lhf) .- (tsf.shf .+ tsf.lhf)) ./ ΔT_sfc

    return nothing
end

function variable_T_saturated_surface_coefs(ts_int, uₕ_int, z_int, z_sfc, T_sfc, params, z0m, z0b)

    # get the near-surface thermal state
    T_int = TD.air_temperature(params, ts_int)
    Rm_int = TD.gas_constant_air(params, ts_int)
    ρ_sfc = TD.air_density(params, ts_int) * (T_sfc / T_int)^(TD.cv_m(params, ts_int) / Rm_int) # use ideal gas law and hydrostatic balance to extrapolate for surface density

    q_sfc = TD.q_vap_saturation_generic(params, T_sfc, ρ_sfc, TD.Liquid()) # TODO: assumes all surface is water covered. Generalize!
    ts_sfc = TD.PhaseEquil_ρTq(params, ρ_sfc, T_sfc, q_sfc)

    # wrap state values
    sc = SF.ValuesOnly{FT}(;
        state_in = SF.InteriorValues(z_int, (uₕ_int.u, uₕ_int.v), ts_int),
        state_sfc = SF.SurfaceValues(z_sfc, (FT(0), FT(0)), ts_sfc),
        z0m = z0m,
        z0b = z0b,
    )

    # calculate all fluxes
    tsf = SF.surface_conditions(params, sc, SF.UniversalFunctions.Businger())

    E = SF.evaporation(sc, params, tsf.Ch)

    return (; shf = tsf.shf, lhf = tsf.lhf, E = E, ρτxz = tsf.ρτxz, ρτyz = tsf.ρτyz)
end

function constant_T_saturated_surface_coefs_coupled(ts_int, uₕ_int, z_int, z_sfc, T_sfc, params, z0m, z0b, Cd, Ch)

    # get the near-surface thermal state
    T_int = TD.air_temperature(params, ts_int)
    Rm_int = TD.gas_constant_air(params, ts_int)
    ρ_sfc = TD.air_density(params, ts_int) * (T_sfc / T_int)^(TD.cv_m(params, ts_int) / Rm_int) # use ideal gas law and hydrostatic balance to extrapolate for surface density

    q_sfc = TD.q_vap_saturation_generic(params, T_sfc, ρ_sfc, TD.Liquid()) # TODO: assumes all surface is water covered. Generalize!
    ts_sfc = TD.PhaseEquil_ρTq(params, ρ_sfc, T_sfc, q_sfc)

    # wrap state values
    sc = SF.Coefficients{FT}(;
        state_in = SF.InteriorValues(z_int, (uₕ_int.u, uₕ_int.v), ts_int),
        state_sfc = SF.SurfaceValues(z_sfc, (FT(0), FT(0)), ts_sfc),
        z0m = z0m, #FT(1e-3),
        z0b = z0b, #FT(1e-5),
        Cd = Cd, #FT(0.001),
        Ch = Ch, #FT(0.0001),
    )

    # calculate all fluxes
    tsf = SF.surface_conditions(params, sc, SF.UniversalFunctions.Businger()) # here can specify tol, maxiter

    E = SF.evaporation(sc, params, tsf.Ch)

    return (; shf = tsf.shf, lhf = tsf.lhf, E = E, ρτxz = tsf.ρτxz, ρτyz = tsf.ρτyz)
end
