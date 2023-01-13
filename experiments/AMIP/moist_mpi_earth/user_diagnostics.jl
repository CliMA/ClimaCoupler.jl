#import Diagnostics: get_var

# Atmos diagnostics

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
    get_var(cs::CoupledSimulation, ::Val{:q_tot}) 

Total specific humidity (kg kg⁻¹).
"""
get_var(cs::CoupledSimulation, ::Val{:q_tot}) =
    cs.model_sims.atmos_sim.integrator.u.c.ρq_tot ./ cs.model_sims.atmos_sim.integrator.u.c.ρ .* float_type_cs(cs)(1000)

"""
    get_var(cs::CoupledSimulation, ::Val{:toa}) 

Top of the atmosphere radiation fluxes (W m⁻²).
"""
function get_var(cs::CoupledSimulation, ::Val{:toa})
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
    swap_space!(radiation_sources, cs.boundary_space)
end

"""
    get_var(cs::CoupledSimulation, ::Val{:precipitation}) 

Precipitation (m m⁻²).
"""
get_var(cs::CoupledSimulation, ::Val{:precipitation}) =
    .-swap_space!(
        cs.model_sims.atmos_sim.integrator.p.col_integrated_rain .+
        cs.model_sims.atmos_sim.integrator.p.col_integrated_snow,
        cs.boundary_space,
    )

# coupler diagnotics
"""
    get_var(cs::CoupledSimulation, ::Val{:T_sfc}) 

Combined surface temperature (K).
"""
get_var(cs::CoupledSimulation, ::Val{:T_sfc}) = swap_space!(cs.fields.T_S, cs.boundary_space)

# land diagnotics
