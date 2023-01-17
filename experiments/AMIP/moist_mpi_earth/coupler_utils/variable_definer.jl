# coupler_diagnostics_list.jl - a temporary prequel to ClimaDiagnostics's diagnostic table

"""
    get_var(cs::CouplerSimulation, ::Symbol)

Defines variable extraction from the coupler simulation.
"""
get_var(cs, ::Val{:T}) = air_temperature(cs)

get_var(cs, ::Val{:u}) = ClimaCore.Geometry.UVVector.(cs.model_sims.atmos_sim.integrator.u.c.uₕ).components.data.:1

get_var(cs, ::Val{:q_tot}) =
    cs.model_sims.atmos_sim.integrator.u.c.ρq_tot ./ cs.model_sims.atmos_sim.integrator.u.c.ρ .* float_type(cs)(1000)

get_var(cs, ::Val{:toa}) = swap_space!(toa_fluxes(cs), cs.boundary_space)

get_var(cs, ::Val{:precipitation}) =
    .-swap_space!(
        cs.model_sims.atmos_sim.integrator.p.col_integrated_rain .+
        cs.model_sims.atmos_sim.integrator.p.col_integrated_snow,
        cs.boundary_space,
    )

get_var(cs, ::Val{:T_sfc}) = swap_space!(cs.fields.T_S, cs.boundary_space)

# more complex calculations
function zonal_wind(cs)
    cs.model_sims.atmos_sim.integrator.u.c.uₕ
end
function air_temperature(cs)
    p = cs.model_sims.atmos_sim.integrator.p
    (; ᶜts, params) = p
    thermo_params = CAP.thermodynamics_params(params)
    ᶜT = @. TD.air_temperature(thermo_params, ᶜts)
end
function toa_fluxes(cs)
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
end
