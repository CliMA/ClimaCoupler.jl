# atmos_init: for ClimaAtmos pre-AMIP interface
import ClimaAtmos
using ClimaAtmos: RRTMGPI
import ClimaCoupler.FluxCalculator: compute_atmos_turbulent_fluxes!
using ClimaCore: Fields.level, Geometry

# the clima atmos `integrator` is now defined
struct ClimaAtmosSimulation{P, Y, D, I} <: AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
function atmos_init(::Type{FT}, parsed_args::Dict) where {FT}

    atmos_config = ClimaAtmos.AtmosConfig(argparse_settings(); parsed_args)
    integrator = ClimaAtmos.get_integrator(atmos_config)
    Y = integrator.u
    center_space = axes(Y.c.ρe_tot)
    face_space = axes(Y.f.u₃)
    spaces = (; center_space = center_space, face_space = face_space)
    if :ρe_int in propertynames(Y.c)
        @warn("Running with ρe_int in coupled mode is not tested yet.")
    end

    # set initial fluxes to zero
    @. integrator.p.sfc_conditions.ρ_flux_h_tot = Geometry.Covariant3Vector(FT(0.0))
    @. integrator.p.sfc_conditions.ρ_flux_q_tot = Geometry.Covariant3Vector(FT(0.0))
    @. integrator.p.sfc_conditions.ρ_flux_uₕ.components = zeros(axes(integrator.p.sfc_conditions.ρ_flux_uₕ.components))
    parent(integrator.p.ᶠradiation_flux) .= parent(zeros(axes(integrator.p.ᶠradiation_flux)))

    ClimaAtmosSimulation(integrator.p.params, Y, spaces, integrator)
end

# required by the Interfacer
get_field(sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux}) =
    Fields.level(sim.integrator.p.ᶠradiation_flux, half)
get_field(sim::ClimaAtmosSimulation, ::Val{:liquid_precipitation}) =
    sim.integrator.p.col_integrated_rain .+ sim.integrator.p.col_integrated_snow # all fallen snow melts for now
get_field(sim::ClimaAtmosSimulation, ::Val{:snow_precipitation}) = sim.integrator.p.col_integrated_snow .* FT(0)

get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_energy_flux}) =
    Geometry.WVector.(sim.integrator.p.sfc_conditions.ρ_flux_h_tot)
get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_moisture_flux}) =
    Geometry.WVector.(sim.integrator.p.sfc_conditions.ρ_flux_q_tot)

function update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_temperature}, field)
    sim.integrator.p.radiation_model.surface_temperature .= RRTMGPI.field2array(field)
end
function update_field!(sim::ClimaAtmosSimulation, ::Val{:albedo}, field)
    sim.integrator.p.radiation_model.diffuse_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(field), 1, length(parent(field)))
    sim.integrator.p.radiation_model.direct_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(field), 1, length(parent(field)))
end

step!(sim::ClimaAtmosSimulation, t) = step!(sim.integrator, t - sim.integrator.t, true)

reinit!(sim::ClimaAtmosSimulation) = reinit!(sim.integrator)

"""
    update_sim!(atmos_sim::ClimaAtmosSimulation, csf)

Updates the surface fields for temperature and albedo.

# Arguments
- `atmos_sim`: [ClimaAtmosSimulation] containing an atmospheric model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(atmos_sim::ClimaAtmosSimulation, csf, turbulent_fluxes)
    update_field!(atmos_sim, Val(:albedo), csf.albedo)
    update_field!(atmos_sim, Val(:surface_temperature), csf.T_S)
end


# flux calculation borrowed from atmos
"""
    CoupledMoninObukhov()
A modified version of a Monin-Obukhov surface for the Coupler, see the link below for more information
https://clima.github.io/SurfaceFluxes.jl/dev/SurfaceFluxes/#Monin-Obukhov-Similarity-Theory-(MOST)
"""
struct CoupledMoninObukhov end
"""
    coupler_surface_setup(::CoupledMoninObukhov, p, csf_sfc = (; T = nothing, z0m = nothing, z0b = nothing, beta = nothing))

Sets up `surface_setup` as a `Fields.Field` of `SurfaceState`s.
"""
function coupler_surface_setup(
    ::CoupledMoninObukhov,
    p;
    csf_sfc = (; T = nothing, z0m = nothing, z0b = nothing, beta = nothing),
)

    surface_state(z0m, z0b, T, beta) = ClimaAtmos.SurfaceConditions.SurfaceState(;
        parameterization = ClimaAtmos.SurfaceConditions.MoninObukhov(; z0m, z0b),
        T,
        beta,
    )
    surface_state_field = @. surface_state(csf_sfc.z0m, csf_sfc.z0b, csf_sfc.T, csf_sfc.beta)
    return surface_state_field
end

"""
    compute_atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)

Computes turbulent surface fluxes using ClimaAtmos's `update_surface_conditions!`. This
requires that we define a new temporary parameter Tuple, `new_p`, and save the new surface state
in it. We do not want `new_p` to live in the atmospheric model permanently, because that would also
trigger flux calculation during Atmos `step!`. We only want to trigger this once per coupling
timestep from ClimaCoupler.
"""
function compute_atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)

    p = atmos_sim.integrator.p

    coupler_sfc_setup = coupler_surface_setup(
        CoupledMoninObukhov(),
        p;
        csf_sfc = (; T = csf.T_S, z0m = csf.z0m_S, z0b = csf.z0b_S, beta = csf.beta),
    )

    p_names = propertynames(p)
    p_values = map(x -> x == :sfc_setup ? coupler_sfc_setup : getproperty(p, x), p_names)

    new_p = (; zip(p_names, p_values)...)

    ClimaAtmos.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t)

    # p.sfc_conditions.ρ_flux_h_tot .= new_p.sfc_conditions.ρ_flux_h_tot
    # p.sfc_conditions.ρ_flux_q_tot .= new_p.sfc_conditions.ρ_flux_q_tot
    # p.sfc_conditions.ρ_flux_uₕ .= new_p.sfc_conditions.ρ_flux_uₕ
    # p.sfc_conditions.ts .= new_p.sfc_conditions.ts
    # p.sfc_conditions.buoyancy_flux .= new_p.sfc_conditions.buoyancy_flux
    # p.sfc_conditions.obukhov_length .= new_p.sfc_conditions.obukhov_length

    p.sfc_conditions .= new_p.sfc_conditions

end
