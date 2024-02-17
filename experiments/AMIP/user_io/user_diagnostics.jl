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
    (; ᶜts) = p.precomputed
    (; params) = p
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
    (; ᶜts) = p.precomputed
    (; params) = p
    thermo_params = CAP.thermodynamics_params(params)
    TD.liquid_specific_humidity.(thermo_params, ᶜts) .* float_type(cs)(1000)
end

"""
    get_var(cs::CoupledSimulation, ::Val{:toa_fluxes})

Top of the atmosphere radiation fluxes (W m⁻²).
"""
function get_var(cs::CoupledSimulation, ::Val{:toa_fluxes})
    atmos_sim = cs.model_sims.atmos_sim
    face_space = axes(atmos_sim.integrator.u.f)
    nz_faces = length(face_space.grid.vertical_grid.topology.mesh.faces)

    LWd_TOA = Fields.level(
        CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation.radiation_model.face_lw_flux_dn), face_space),
        nz_faces - half,
    )
    LWu_TOA = Fields.level(
        CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation.radiation_model.face_lw_flux_up), face_space),
        nz_faces - half,
    )
    SWd_TOA = Fields.level(
        CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation.radiation_model.face_sw_flux_dn), face_space),
        nz_faces - half,
    )
    SWu_TOA = Fields.level(
        CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation.radiation_model.face_sw_flux_up), face_space),
        nz_faces - half,
    )

    radiation_sources = @. -(LWd_TOA + SWd_TOA - LWu_TOA - SWu_TOA)
    swap_space!(zeros(cs.boundary_space), radiation_sources)
end

"""
    get_var(cs::CoupledSimulation, ::Val{:precipitation_rate})

Precipitation rate (Kg m⁻² s⁻¹).
"""
get_var(cs::CoupledSimulation, ::Val{:precipitation_rate}) =
    .-swap_space!(
        zeros(cs.boundary_space),
        cs.model_sims.atmos_sim.integrator.p.precipitation.col_integrated_rain .+
        cs.model_sims.atmos_sim.integrator.p.precipitation.col_integrated_snow,
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

# land diagnotics
