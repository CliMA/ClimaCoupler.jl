"""
    AbstractSurfaceStub

An abstract type representing a simple stand-in surface model.
Any concrete type extending this abstract type should have a `cache` field that
contains the necessary fields for the simulation, at minimum the following:
    - `T_sfc` (surface temperature [K])
    - `z0m` (roughness length for momentum [m])
    - `z0b` (roughness length for tracers [m])
    - `α_direct` (direct albedo)
    - `α_diffuse` (diffuse albedo)
    - `area_fraction` (fraction of the grid cell covered by the ocean)
    - `phase` (phase of the water used to calculate surface humidity)
    - `thermo_params` (thermodynamic parameters)
"""
abstract type AbstractSurfaceStub <: AbstractSurfaceSimulation end

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
get_field(sim::AbstractSurfaceStub, ::Val{:energy}) = nothing
get_field(sim::AbstractSurfaceStub, ::Val{:roughness_buoyancy}) = sim.cache.z0b
get_field(sim::AbstractSurfaceStub, ::Val{:roughness_momentum}) = sim.cache.z0m
get_field(sim::AbstractSurfaceStub, ::Val{:surface_direct_albedo}) = sim.cache.α_direct
get_field(sim::AbstractSurfaceStub, ::Val{:surface_diffuse_albedo}) = sim.cache.α_diffuse
get_field(sim::AbstractSurfaceStub, ::Val{:surface_temperature}) = sim.cache.T_sfc

"""
    update_field!(sim::AbstractSurfaceStub, ::Val{:area_fraction}, field::CC.Fields.Field)

Updates the specified value in the cache of `SurfaceStub`.
"""
function update_field!(
    sim::AbstractSurfaceStub,
    ::Val{:area_fraction},
    field::CC.Fields.Field,
)
    sim.cache.area_fraction .= field
    return nothing
end
function update_field!(
    sim::AbstractSurfaceStub,
    ::Val{:surface_temperature},
    field::CC.Fields.Field,
)
    Interfacer.remap!(sim.cache.T_sfc, field)
end
function update_field!(
    sim::AbstractSurfaceStub,
    ::Val{:surface_direct_albedo},
    field::CC.Fields.Field,
)
    Interfacer.remap!(sim.cache.α_direct, field)
end
function update_field!(
    sim::AbstractSurfaceStub,
    ::Val{:surface_diffuse_albedo},
    field::CC.Fields.Field,
)
    Interfacer.remap!(sim.cache.α_diffuse, field)
end
update_field!(::AbstractSurfaceStub, ::Val{:liquid_precipitation}, field) = nothing
update_field!(::AbstractSurfaceStub, ::Val{:SW_d}, field) = nothing
update_field!(::AbstractSurfaceStub, ::Val{:LW_d}, field) = nothing
update_field!(::AbstractSurfaceStub, ::Val{:snow_precipitation}, field) = nothing
update_field!(::AbstractSurfaceStub, ::Val{:turbulent_energy_flux}, field) = nothing
update_field!(::AbstractSurfaceStub, ::Val{:turbulent_moisture_flux}, field) = nothing

## Extensions of FieldExchanger.jl functions

"""
    step!(::AbstractSurfaceStub, t)

The stub surface simulation is not updated by this function. Extends `SciMLBase.step!`.
"""
step!(::AbstractSurfaceStub, _) = nothing
