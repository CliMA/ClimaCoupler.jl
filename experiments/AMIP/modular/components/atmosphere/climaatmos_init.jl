# atmos_init: for ClimaAtmos pre-AMIP interface
import ClimaAtmos
import ClimaCoupler.FluxCalculator: compute_atmos_turbulent_fluxes!

driver_file = joinpath(pkgdir(ClimaAtmos), "examples", "hybrid", "driver.jl")
ENV["CI_PERF_SKIP_RUN"] = true
try
    include(driver_file)
catch err
    if err.error !== :exit_profile
        rethrow(err.error)
    end
end
# the clima atmos `integrator` is now defined
struct ClimaAtmosSimulation{P, Y, D, I} <: AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function atmos_init(::Type{FT}, Y, integrator; params = nothing) where {FT}
    center_space = axes(Y.c.ρe_tot)
    face_space = axes(Y.f.w)
    spaces = (; center_space = center_space, face_space = face_space)
    if :ρe_int in propertynames(Y.c)
        @warn("Running with ρe_int in coupled mode is not tested yet.")
    end

    ClimaAtmosSimulation(params, Y, spaces, integrator)
end

# required by the Interfacer
get_field(sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux}) = level(sim.integrator.p.ᶠradiation_flux, half)
get_field(sim::ClimaAtmosSimulation, ::Val{:liquid_precipitation}) =
    sim.integrator.p.col_integrated_rain .+ sim.integrator.p.col_integrated_snow # all fallen snow melts for now
get_field(sim::ClimaAtmosSimulation, ::Val{:snow_precipitation}) = sim.integrator.p.col_integrated_snow .* FT(0)

get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_energy_flux}) = sim.integrator.p.ρ_dif_flux_h_tot
get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_moisture_flux}) = sim.integrator.p.ρ_dif_flux_q_tot

function update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_temperature}, field)
    sim.integrator.p.radiation_model.surface_temperature .= RRTMGPI.field2array(field)
    parent(sim.integrator.p.T_sfc) .= parent(field)

end
function update_field!(sim::ClimaAtmosSimulation, ::Val{:albedo}, field)
    sim.integrator.p.radiation_model.diffuse_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(field), 1, length(parent(field)))
    sim.integrator.p.radiation_model.direct_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(field), 1, length(parent(field)))
end
function update_field!(sim::ClimaAtmosSimulation, ::Val{:roughness_momentum}, field)
    parent(sim.integrator.p.sfc_inputs.z0m) .= parent(field)
end
function update_field!(sim::ClimaAtmosSimulation, ::Val{:roughness_buoyancy}, field)
    parent(sim.integrator.p.sfc_inputs.z0b) .= parent(field)
end

step!(sim::ClimaAtmosSimulation, t) = step!(sim.integrator, t - sim.integrator.t, true)

reinit!(sim::ClimaAtmosSimulation) = reinit!(sim.integrator)


# flux calculation borrowed from atmos
function compute_atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)
    # TODO: dependence on atmos will be removed

    thermo_params = CAP.thermodynamics_params(atmos_sim.integrator.p.params)

    # calculate turbulent fluxes on atmos grid and save in atmos cache
    update_field!(atmos_sim, Val(:surface_temperature), csf.T_S)

    if :z0b in propertynames(atmos_sim.integrator.p.surface_scheme)
        update_field!(atmos_sim, Val(:roughness_momentum), csf.z0m_S)
        update_field!(atmos_sim, Val(:roughness_buoyancy), csf.z0b_S)
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
end
