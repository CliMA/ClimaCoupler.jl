#import Diagnostics: get_var

# Atmos diagnostics
import ClimaAtmos.Parameters as CAP
import Thermodynamics as TD
using ClimaCoupler.Interfacer: CoupledSimulation, float_type

"""
    get_var(cs::CoupledSimulation, ::Val{:T})

Air temperature (K).
"""
function get_var(cs::CoupledSimulation, ::Val{:T})
    p = cs.model_sims.atmos_sim.integrator.p
    (; ᶜts, params) = p
    thermo_params = CAP.thermodynamics_params(params)
    @. TD.air_temperature(thermo_params, ᶜts)
end

"""
    get_var(cs::CoupledSimulation, ::Val{:u})

Zonal wind (m s⁻¹).
"""
get_var(cs::CoupledSimulation, ::Val{:u}) =
    ClimaCore.Geometry.UVVector.(cs.model_sims.atmos_sim.integrator.u.c.uₕ).components.data.:1

"""
    get_var(cs::CoupledSimulation, ::Val{:v})

Meridional wind (m s⁻¹).
"""
get_var(cs::CoupledSimulation, ::Val{:v}) =
    ClimaCore.Geometry.UVVector.(cs.model_sims.atmos_sim.integrator.u.c.uₕ).components.data.:2

"""
    get_var(cs::CoupledSimulation, ::Val{:w})

Vertical wind (m s⁻¹).
"""
get_var(cs::CoupledSimulation, ::Val{:w}) =
    ClimaCore.Geometry.WVector.(CA.ᶜinterp.(cs.model_sims.atmos_sim.integrator.u.f.u₃)).components.data.:1

"""
    get_var(cs::CoupledSimulation, ::Val{:moist_static_energy})

Moist static energy. [J / kg]
"""
function get_var(cs::CoupledSimulation, ::Val{:moist_static_energy})
    p = cs.model_sims.atmos_sim.integrator.p
    (; ᶜts, params) = p
    c_space = axes(cs.model_sims.atmos_sim.integrator.u.c)
    thermo_params = CAP.thermodynamics_params(params)
    e_pot = 9.81 .* Fields.coordinate_field(c_space).z
    TD.moist_static_energy.(thermo_params, ᶜts, e_pot)
end

"""
    get_var(cs::CoupledSimulation, ::Val{:lapse_rate})

Lapse rate (K/m).
"""
function get_var(cs::CoupledSimulation, ::Val{:lapse_rate})
    p = cs.model_sims.atmos_sim.integrator.p
    (; ᶜts, params) = p
    thermo_params = CAP.thermodynamics_params(params)
    ᶜT = @. TD.air_temperature(thermo_params, ᶜts)
    ClimaCore.Geometry.WVector.(CA.ᶜgradᵥ.(CA.ᶠinterp.(ᶜT))).components.data.:1
end

"""
    get_var(cs::CoupledSimulation, ::Val{:eddy_diffusivity})

Eddy diffusivity. [m2/s]
"""
function get_var(cs::CoupledSimulation, ::Val{:eddy_diffusivity})
    p = cs.model_sims.atmos_sim.integrator.p
    Y = cs.model_sims.atmos_sim.integrator.u
    (; ᶜp) = p # assume ᶜts and ᶜp have been updated
    (; C_E) = p.atmos.vert_diff

    interior_uₕ = Fields.level(Y.c.uₕ, 1)
    ᶠp = ᶠK_E = p.ᶠtemp_scalar
    Fields.bycolumn(axes(ᶜp)) do colidx
        @. ᶠp[colidx] = CA.ᶠinterp(ᶜp[colidx])
        ᶜΔz_surface = Fields.Δz_field(interior_uₕ)
        @. ᶠK_E[colidx] = CA.eddy_diffusivity_coefficient(
                C_E,
                CA.norm(interior_uₕ[colidx]),
                ᶜΔz_surface[colidx] / 2,
                ᶠp[colidx],
            )
    end
    return CA.ᶜinterp.(ᶠK_E)
end

"""
    get_var(cs::CoupledSimulation, ::Val{:q_tot})

Total specific humidity (g kg⁻¹).
"""
get_var(cs::CoupledSimulation, ::Val{:q_tot}) =
    cs.model_sims.atmos_sim.integrator.u.c.ρq_tot ./ cs.model_sims.atmos_sim.integrator.u.c.ρ .* float_type(cs)(1000)

"""
    get_var(cs::CoupledSimulation, ::Val{:q_liq_ice})

Cloud specific humidity (g kg⁻¹).
"""
function get_var(cs::CoupledSimulation, ::Val{:q_liq_ice})
    p = cs.model_sims.atmos_sim.integrator.p
    (; ᶜts, params) = p
    thermo_params = CAP.thermodynamics_params(params)
    TD.liquid_specific_humidity.(thermo_params, ᶜts) .* float_type(cs)(1000) .+     TD.ice_specific_humidity.(thermo_params, ᶜts) .* float_type(cs)(1000)
end

# need to integrate between z and z0
# """
#     get_var(cs::CoupledSimulation, ::Val{:mass_streamfunction})

# Mass streamfunction (kg s⁻¹).
# """
# function get_var(cs::CoupledSimulation, ::Val{:mass_streamfunction})
#     v = ClimaCore.Geometry.UVVector.(cs.model_sims.atmos_sim.integrator.u.c.uₕ).components.data.:2
#     ρ = cs.model_sims.atmos_sim.integrator.u.c.ρ
#     c_space = axes(ρ)
#     coslat = cos.(Fields.coordinate_field(c_space).lat .* π ./ 180)
#     strf = @. vert_int(2π * 6.371e6 * coslat * v * ρ)
#     return strf
# end

"""
    get_var(cs::CoupledSimulation, ::Val{:toa_fluxes})

Top of the atmosphere radiation fluxes (W m⁻²).
"""
function get_var(cs::CoupledSimulation, ::Val{:toa_fluxes})
    atmos_sim = cs.model_sims.atmos_sim
    face_space = axes(atmos_sim.integrator.u.f)
    z = parent(Fields.coordinate_field(face_space).z)
    Δz_top = round(FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1]))
    n_faces = length(z[:, 1, 1, 1, 1])

    LWd_TOA = Fields.level(
        RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_dn), face_space),
        n_faces - half,
    )
    LWu_TOA = Fields.level(
        RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_up), face_space),
        n_faces - half,
    )
    SWd_TOA = Fields.level(
        RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_dn), face_space),
        n_faces - half,
    )
    SWu_TOA = Fields.level(
        RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_up), face_space),
        n_faces - half,
    )

    radiation_sources = @. -(LWd_TOA + SWd_TOA - LWu_TOA - SWu_TOA)
    swap_space!(zeros(cs.boundary_space), radiation_sources)
end

"""
    get_var(cs::CoupledSimulation, ::Val{:longwave_down_sfc})

Downward shortwave radiation fluxes (W m⁻²).
"""
function get_var(cs::CoupledSimulation, ::Val{:longwave_down_sfc})
    atmos_sim = cs.model_sims.atmos_sim
    face_space = axes(atmos_sim.integrator.u.f)
    z = parent(Fields.coordinate_field(face_space).z)
    Δz_top = round(FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1]))
    n_faces = length(z[:, 1, 1, 1, 1])

    LWd_TOA = Fields.level(
        RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_dn), face_space),
        half,
    )
    swap_space!(zeros(cs.boundary_space), LWd_TOA)
end

"""
    get_var(cs::CoupledSimulation, ::Val{:longwave_up_sfc})

Upward shortwave radiation fluxes (W m⁻²).
"""
function get_var(cs::CoupledSimulation, ::Val{:longwave_up_sfc})
    atmos_sim = cs.model_sims.atmos_sim
    face_space = axes(atmos_sim.integrator.u.f)
    z = parent(Fields.coordinate_field(face_space).z)
    Δz_top = round(FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1]))
    n_faces = length(z[:, 1, 1, 1, 1])

    LWu_TOA = Fields.level(
        RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_up), face_space),
        half,
    )
    swap_space!(zeros(cs.boundary_space), LWu_TOA)
end

"""
    get_var(cs::CoupledSimulation, ::Val{:shortwave_down_sfc})

Downward shortwave radiation fluxes (W m⁻²).
"""
function get_var(cs::CoupledSimulation, ::Val{:shortwave_down_sfc})
    atmos_sim = cs.model_sims.atmos_sim
    face_space = axes(atmos_sim.integrator.u.f)
    z = parent(Fields.coordinate_field(face_space).z)
    Δz_top = round(FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1]))
    n_faces = length(z[:, 1, 1, 1, 1])

    SWd_TOA = Fields.level(
        RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_dn), face_space),
        half,
    )
    swap_space!(zeros(cs.boundary_space), SWd_TOA)
end

"""
    get_var(cs::CoupledSimulation, ::Val{:shortwave_up_sfc})

Upward shortwave radiation fluxes (W m⁻²).
"""
function get_var(cs::CoupledSimulation, ::Val{:shortwave_up_sfc})
    atmos_sim = cs.model_sims.atmos_sim
    face_space = axes(atmos_sim.integrator.u.f)
    z = parent(Fields.coordinate_field(face_space).z)
    Δz_top = round(FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1]))
    n_faces = length(z[:, 1, 1, 1, 1])

    SWu_TOA = Fields.level(
        RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_up), face_space),
        half,
    )
    swap_space!(zeros(cs.boundary_space), SWu_TOA)
end


"""
    get_var(cs::CoupledSimulation, ::Val{:precipitation_rate})

Precipitation rate (Kg m⁻² s⁻¹).
"""
get_var(cs::CoupledSimulation, ::Val{:precipitation_rate}) =
    .-swap_space!(
        zeros(cs.boundary_space),
        cs.model_sims.atmos_sim.integrator.p.col_integrated_rain .+
        cs.model_sims.atmos_sim.integrator.p.col_integrated_snow,
    )

# coupler diagnotics
"""
    get_var(cs::CoupledSimulation, ::Val{:T_sfc})

Combined surface temperature (K).
"""
get_var(cs::CoupledSimulation, ::Val{:T_sfc}) = swap_space!(zeros(cs.boundary_space), cs.fields.T_S)

"""
    get_var(cs::CoupledSimulation, ::Val{:tubulent_energy_fluxes})

Combined aerodynamic turbulent energy surface fluxes (W m⁻²).
"""
get_var(cs::CoupledSimulation, ::Val{:tubulent_energy_fluxes}) =
    swap_space!(zeros(cs.boundary_space), cs.fields.F_turb_energy)


"""
    get_var(cs::CoupledSimulation, ::Val{:evaporation})

Combined aerodynamic turbulent moisture surface fluxes (kg m⁻² s-1).
"""
get_var(cs::CoupledSimulation, ::Val{:evaporation}) =
    swap_space!(zeros(cs.boundary_space), cs.fields.F_turb_moisture)



