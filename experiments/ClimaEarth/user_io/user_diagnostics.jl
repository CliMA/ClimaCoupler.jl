import ClimaCore as CC
import ClimaAtmos.Parameters as CAP
import Thermodynamics as TD
import ClimaCoupler: Diagnostics, Interfacer, Utilities

"""
    Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:T})

Air temperature (K).
"""
function Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:T})
    p = cs.model_sims.atmos_sim.integrator.p
    (; ᶜts) = p.precomputed
    (; params) = p
    thermo_params = CAP.thermodynamics_params(params)
    @. TD.air_temperature(thermo_params, ᶜts)
end

"""
    Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:u})

Zonal wind (m s⁻¹).
"""
Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:u}) =
    CC.Geometry.UVVector.(cs.model_sims.atmos_sim.integrator.u.c.uₕ).components.data.:1

"""
    Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:q_tot})

Total specific humidity (g kg⁻¹).
"""
Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:q_tot}) =
    cs.model_sims.atmos_sim.integrator.u.c.ρq_tot ./ cs.model_sims.atmos_sim.integrator.u.c.ρ .*
    Interfacer.float_type(cs)(1000)


"""
    Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:q_liq_ice})

Cloud specific humidity (g kg⁻¹).
"""
function Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:q_liq_ice})
    p = cs.model_sims.atmos_sim.integrator.p
    (; ᶜts) = p.precomputed
    (; params) = p
    thermo_params = CAP.thermodynamics_params(params)
    TD.liquid_specific_humidity.(thermo_params, ᶜts) .* Interfacer.float_type(cs)(1000)
end

"""
    Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:toa_fluxes})

Top of the atmosphere radiation fluxes (W m⁻²).
"""
function Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:toa_fluxes})
    atmos_sim = cs.model_sims.atmos_sim
    face_space = axes(atmos_sim.integrator.u.f)
    nz_faces = length(face_space.grid.vertical_grid.topology.mesh.faces)

    LWd_TOA = CC.Fields.level(
        CC.Fields.array2field(FT.(atmos_sim.integrator.p.radiation.rrtmgp_model.face_lw_flux_dn), face_space),
        nz_faces - CC.Utilities.half,
    )
    LWu_TOA = CC.Fields.level(
        CC.Fields.array2field(FT.(atmos_sim.integrator.p.radiation.rrtmgp_model.face_lw_flux_up), face_space),
        nz_faces - CC.Utilities.half,
    )
    SWd_TOA = CC.Fields.level(
        CC.Fields.array2field(FT.(atmos_sim.integrator.p.radiation.rrtmgp_model.face_sw_flux_dn), face_space),
        nz_faces - CC.Utilities.half,
    )
    SWu_TOA = CC.Fields.level(
        CC.Fields.array2field(FT.(atmos_sim.integrator.p.radiation.rrtmgp_model.face_sw_flux_up), face_space),
        nz_faces - CC.Utilities.half,
    )

    radiation_sources = @. -(LWd_TOA + SWd_TOA - LWu_TOA - SWu_TOA)
    Utilities.swap_space!(cs.boundary_space, radiation_sources)
end

"""
    Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:precipitation_rate})

Precipitation rate (Kg m⁻² s⁻¹).
"""
Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:precipitation_rate}) =
    .-Utilities.swap_space!(
        cs.boundary_space,
        cs.model_sims.atmos_sim.integrator.p.precipitation.surface_rain_flux .+
        cs.model_sims.atmos_sim.integrator.p.precipitation.surface_snow_flux,
    )

# coupler diagnotics
"""
    Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:T_sfc})

Combined surface temperature (K).
"""
Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:T_sfc}) =
    Utilities.swap_space!(cs.boundary_space, cs.fields.T_S)

"""
    Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:tubulent_energy_fluxes})

Combined aerodynamic turbulent energy surface fluxes (W m⁻²).
"""
Diagnostics.get_var(cs::Interfacer.CoupledSimulation, ::Val{:tubulent_energy_fluxes}) =
    Utilities.swap_space!(cs.boundary_space, cs.fields.F_turb_energy)

# land diagnotics
