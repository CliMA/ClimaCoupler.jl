using ClimaCore.Geometry: ⊗
"""
calculate_surface_fluxes_atmos_grid!(integrator)

- TODO: generalize interface for regridding and take land state out of atmos's integrator.p
"""
function calculate_surface_fluxes_atmos_grid!(integrator, T_sfc)
    p = integrator.p
    (; ᶜts, dif_flux_energy, dif_flux_ρq_tot, dif_flux_uₕ, params) = p

    Y = integrator.u

    # Turbulent surface flux calculation
    tsf =
        variable_T_saturated_surface_coefs.(
            Spaces.level(ᶜts, 1),
            Geometry.UVVector.(Spaces.level(Y.c.uₕ, 1)),
            Spaces.level(Fields.coordinate_field(Y.c).z, 1),
            FT(0), # TODO: get actual value of z_sfc
            T_sfc,
            params,
            FT(0.001), # TODO: get these roughness lengths from land
            FT(0.001),
        )

    # Radiation fluxes 
    t = integrator.t
    Rn = FT(10 * sin(t / 2π * 1000)) # TODO: link to the SW and LW fluxes from atmos

    # Total energy flux
    if :ρe in propertynames(Y.c)
        @. dif_flux_energy = Geometry.WVector(tsf.shf + tsf.shf + Rn)
    end

    # Moisture mass flux
    if :ρq_tot in propertynames(Y.c)
        @. dif_flux_ρq_tot = Geometry.WVector(tsf.E)
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

    _ρ_liq = FT(1e3)# TODO: use CLIMAParameters
    E = SF.evaporation(sc, params, tsf.Ch) / _ρ_liq
    return (; shf = tsf.shf, lhf = tsf.lhf, E = E, ρτxz = tsf.ρτxz, ρτyz = tsf.ρτyz)
end
