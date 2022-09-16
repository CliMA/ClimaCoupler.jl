using ClimaCore.Geometry: ⊗
using ClimaCore.Utilities: half, PlusHalf

"""
    set_ρ_sfc!(ρ_sfc, T_S, integrator)

sets the value of the ρ_sfc field based on the temperature of the surface, 
the temperature of the atmosphere at the lowest level, and the heigh
of the lowest level.
"""
function set_ρ_sfc!(ρ_sfc, T_S, integrator)
    ts = integrator.p.ᶜts
    thermo_params = CAP.thermodynamics_params(integrator.p.params)
    ts_int = Spaces.level(ts, 1)
    parent(ρ_sfc) .= parent(ρ_sfc_at_point.(thermo_params, ts_int, swap_space!(T_S, axes(ts_int))))
end

"""
    ρ_sfc_at_point(params, ts_int, T_sfc)

Computes the surface density at a point given the atmospheric state
at the lowest level, the surface temperature, and the assumption of
an ideal gas and hydrostatic balance.

Required because the surface models do not compute air density as a 
variable.
"""
function ρ_sfc_at_point(params, ts_int, T_sfc)
    T_int = TD.air_temperature(params, ts_int)
    Rm_int = TD.gas_constant_air(params, ts_int)
    ρ_air = TD.air_density(params, ts_int)
    ρ_sfc = ρ_air * (T_sfc / T_int)^(TD.cv_m(params, ts_int) / Rm_int)  # use ideal gas law and hydrostatic balance to extrapolate for surface density
    return ρ_sfc
end

"""
calculate_surface_fluxes_atmos_grid!(integrator)

- TODO: generalize interface for regridding and take land state out of atmos's integrator.p
"""
function calculate_surface_fluxes_atmos_grid!(integrator, info_sfc)
    p = integrator.p
    (; ᶜts, dif_flux_energy, dif_flux_ρq_tot, dif_flux_uₕ, params, Cd, Ch) = p

    (; T_sfc, ρ_sfc, q_sfc, z0m, z0b, ice_mask) = info_sfc
    Y = integrator.u
    FT = eltype(integrator.u.c.ρ)
    thermo_params = CAP.thermodynamics_params(integrator.p.params)
    surface_flux_params = CAP.surface_fluxes_params(integrator.p.params)
    # Turbulent surface flux calculation

    tsf =
        constant_T_saturated_surface_coefs_coupled.(
            Spaces.level(ᶜts, 1),
            Geometry.UVVector.(Spaces.level(Y.c.uₕ, 1)),
            Spaces.level(Fields.coordinate_field(Y.c).z, 1),
            FT(0), # TODO: get actual value of z_sfc
            swap_space!(T_sfc, axes(Spaces.level(Y.c, 1))), # remove when same instance issue is resolved
            swap_space!(ρ_sfc, axes(Spaces.level(Y.c, 1))), # remove when same instance issue is resolved
            swap_space!(q_sfc, axes(Spaces.level(Y.c, 1))), # remove when same instance issue is resolved
            thermo_params,
            surface_flux_params,
            swap_space!(z0m, axes(Spaces.level(Y.c, 1))), # TODO: get these roughness lengths from land
            swap_space!(z0b, axes(Spaces.level(Y.c, 1))),
            Cd,
            Ch,
        )

    # Total energy flux
    if :ρe_tot in propertynames(Y.c)

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

    return nothing
end

function constant_T_saturated_surface_coefs_coupled(
    ts_int,
    uₕ_int,
    z_int::FT,
    z_sfc::FT,
    T_sfc::FT,
    ρ_sfc::FT,
    q_sfc::FT,
    thermo_params,
    surface_flux_params,
    z0m::FT,
    z0b::FT,
    Cd::FT,
    Ch::FT,
) where {FT}
    # get the near-surface thermal state
    ts_sfc = TD.PhaseEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q_sfc)

    # wrap state values

    sc = SF.Coefficients{FT}(;
        state_in = SF.InteriorValues(z_int, (uₕ_int.u, uₕ_int.v), ts_int),
        state_sfc = SF.SurfaceValues(z_sfc, (FT(0), FT(0)), ts_sfc),
        z0m = z0m,
        z0b = z0b,
        Cd = Cd,
        Ch = Ch,
    )

    # calculate all fluxes
    tsf = SF.surface_conditions(surface_flux_params, sc)

    E = SF.evaporation(surface_flux_params, sc, tsf.Ch)

    return (; shf = tsf.shf, lhf = tsf.lhf, E = E, ρτxz = tsf.ρτxz, ρτyz = tsf.ρτyz)
end
