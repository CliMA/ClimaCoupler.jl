"""
    AbstractSurfaceStub

An abstract type representing a simple stand-in surface model.
Any concrete type extending this abstract type should have a `cache` field that
contains the necessary fields for the simulation, at minimum the following:
    - `T_sfc` (surface temperature [K])
    - `ρ_sfc` (surface air density [kg / m3])
    - `z0m` (roughness length for momentum [m])
    - `z0b` (roughness length for tracers [m])
    - `beta` (evaporation scaling factor)
    - `α_direct` (direct albedo)
    - `α_diffuse` (diffuse albedo)
    - `area_fraction` (fraction of the grid cell covered by the ocean)
    - `phase` (phase of the water used to calculate surface humidity)
    - `thermo_params` (thermodynamic parameters)
"""
abstract type AbstractSurfaceStub <: SurfaceModelSimulation end

"""
    SurfaceStub

On object containing simulation-like info, used as a stub or for prescribed data.
"""
struct SurfaceStub{I} <: AbstractSurfaceStub
    cache::I
end

## Extensions of Interfacer.jl functions

"""
    get_field(::AbstractSurfaceStub, ::Val)

A getter function, that should not allocate. If undefined, it returns a descriptive error.
"""
get_field(sim::AbstractSurfaceStub, ::Val{:area_fraction}) = sim.cache.area_fraction
get_field(sim::AbstractSurfaceStub, ::Val{:beta}) = sim.cache.beta
get_field(sim::AbstractSurfaceStub, ::Val{:energy}) = nothing
get_field(sim::AbstractSurfaceStub, ::Val{:roughness_buoyancy}) = sim.cache.z0b
get_field(sim::AbstractSurfaceStub, ::Val{:roughness_momentum}) = sim.cache.z0m
get_field(sim::AbstractSurfaceStub, ::Val{:surface_direct_albedo}) = sim.cache.α_direct
get_field(sim::AbstractSurfaceStub, ::Val{:surface_diffuse_albedo}) = sim.cache.α_diffuse
get_field(sim::AbstractSurfaceStub, ::Val{:surface_humidity}) =
    TD.q_vap_saturation_generic.(sim.cache.thermo_params, sim.cache.T_sfc, sim.cache.ρ_sfc, sim.cache.phase)
get_field(sim::AbstractSurfaceStub, ::Val{:surface_temperature}) = sim.cache.T_sfc
get_field(sim::AbstractSurfaceStub, ::Val{:water}) = nothing

"""
    update_field!(sim::AbstractSurfaceStub, ::Val{:area_fraction}, field::CC.Fields.Field)

Updates the specified value in the cache of `SurfaceStub`.
"""
function update_field!(sim::AbstractSurfaceStub, ::Val{:area_fraction}, field::CC.Fields.Field)
    sim.cache.area_fraction .= field
end
function update_field!(sim::AbstractSurfaceStub, ::Val{:surface_temperature}, field::CC.Fields.Field)
    sim.cache.T_sfc .= field
end
function update_field!(sim::AbstractSurfaceStub, ::Val{:air_density}, field)
    parent(sim.cache.ρ_sfc) .= parent(field)
end
function update_field!(sim::AbstractSurfaceStub, ::Val{:surface_direct_albedo}, field::CC.Fields.Field)
    sim.cache.α_direct .= field
end
function update_field!(sim::AbstractSurfaceStub, ::Val{:surface_diffuse_albedo}, field::CC.Fields.Field)
    sim.cache.α_diffuse .= field
end
update_field!(::AbstractSurfaceStub, ::Val{:liquid_precipitation}, field) = nothing
update_field!(::AbstractSurfaceStub, ::Val{:radiative_energy_flux_sfc}, field) = nothing
update_field!(::AbstractSurfaceStub, ::Val{:snow_precipitation}, field) = nothing
update_field!(::AbstractSurfaceStub, ::Val{:turbulent_energy_flux}, field) = nothing
update_field!(::AbstractSurfaceStub, ::Val{:turbulent_moisture_flux}, field) = nothing

"""
    name(::ComponentModelSimulation)

Returns simulation name, if defined, or `Unnamed` if not.
"""
name(::AbstractSurfaceStub) = "SurfaceStub"


## Extensions of FieldExchanger.jl functions

"""
    update_sim!(::AbstractSurfaceStub, csf, area_fraction)

The stub surface simulation only updates the air density (needed for the turbulent flux calculation).
"""
function update_sim!(sim::AbstractSurfaceStub, csf, area_fraction)
    update_field!(sim, Val(:air_density), csf.ρ_sfc)
end

"""
    reinit!(cs::AbstractSurfaceStub)

The stub surface simulation is not updated by this function. Extends `SciMLBase.reinit!`.
"""
reinit!(::AbstractSurfaceStub) = nothing

"""
Extend Interfacer.add_coupler_fields! to add the fields required for AbstractSurfaceStub.

The fields added are:
- `:ρ_sfc` (for humidity calculation)
"""
function Interfacer.add_coupler_fields!(coupler_field_names, ::AbstractSurfaceStub)
    surface_coupler_fields = [:ρ_sfc]
    push!(coupler_field_names, surface_coupler_fields...)
end

"""
    step!(::AbstractSurfaceStub, t)

The stub surface simulation is not updated by this function. Extends `SciMLBase.step!`.
"""
step!(::AbstractSurfaceStub, _) = nothing
