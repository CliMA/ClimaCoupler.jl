"""
    Interfacer

This modules contains abstract types, interface templates and model stubs for coupling component models.
"""
module Interfacer
import Thermodynamics as TD

using ClimaCore: Fields
export ComponentModelSimulation,
    AtmosModelSimulation,
    SurfaceModelSimulation,
    SurfaceStub,
    SeaIceModelSimulation,
    LandModelSimulation,
    OceanModelSimulation,
    name,
    get_field,
    update_field!

"""
    ComponentModelSimulation

An abstract type encompassing all component model (and model stub) simulations.
"""
abstract type ComponentModelSimulation end

"""
    AtmosModelSimulation

An abstract type for an atmospheric model simulation.
"""
abstract type AtmosModelSimulation <: ComponentModelSimulation end

"""
    SurfaceModelSimulation

An abstract type for surface model simulations.
"""
abstract type SurfaceModelSimulation <: ComponentModelSimulation end

abstract type SeaIceModelSimulation <: SurfaceModelSimulation end
abstract type LandModelSimulation <: SurfaceModelSimulation end
abstract type OceanModelSimulation <: SurfaceModelSimulation end


"""
    SurfaceStub

On object containing simulation-like info, used as a stub or for prescribed data.
"""
struct SurfaceStub{I} <: SurfaceModelSimulation
    cache::I
end

"""
    get_field(::SurfaceStub, ::Val)

A getter function, that should not allocate. If undefined, it returns a descriptive error.
"""
get_field(sim::SurfaceStub, ::Val{:area_fraction}) = sim.cache.area_fraction
get_field(sim::SurfaceStub, ::Val{:surface_temperature}) = sim.cache.T_sfc
get_field(sim::SurfaceStub, ::Val{:albedo}) = sim.cache.α
get_field(sim::SurfaceStub, ::Val{:roughness_momentum}) = sim.cache.z0m
get_field(sim::SurfaceStub, ::Val{:roughness_buoyancy}) = sim.cache.z0b
get_field(sim::SurfaceStub, ::Val{:beta}) = sim.cache.beta
get_field(sim::SurfaceStub, ::Val{:surface_humidity}) =
    TD.q_vap_saturation_generic.(sim.cache.thermo_params, sim.cache.T_sfc, sim.cache.ρ_sfc, sim.cache.phase)

function get_field(sim::ComponentModelSimulation, val::Val)
    error("undefined get_field for $val for " * name(sim))
end

get_field(sim::SurfaceStub, ::Val{:energy}) = nothing
get_field(sim::SurfaceStub, ::Val{:water}) = nothing

"""
    get_field(::ComponentModelSimulation, ::Val, colidx::Fields.ColumnIndex)

Extension of `get_field(::ComponentModelSimulation, ::Val)`, indexing into the specified colum index.
"""
function get_field(sim::ComponentModelSimulation, val::Val, colidx::Fields.ColumnIndex)
    if get_field(sim, val) isa AbstractFloat
        get_field(sim, val)
    else
        get_field(sim, val)[colidx]
    end
end

"""
    update_field!(::ComponentModelSimulation, ::Val, _...)

No update in unspecified in the particular component model simulation.
"""
update_field!(sim::ComponentModelSimulation, val::Val, _...) = nothing

# TODO:
# function update_field!(sim::ComponentModelSimulation, val::Val, _...)
#     warning = Warning("undefined `update!` for $val in " * name(sim) * ": skipping")
#     @warn(warning.message, maxlog=10)
#     return warning
# end
# struct Warning
#     message::String
# end

"""
    update_field!(sim::SurfaceStub, ::Val{:area_fraction}, field::Fields.Field)

Updates the specified value in the cache of `SurfaceStub`.
"""
function update_field!(sim::SurfaceStub, ::Val{:area_fraction}, field::Fields.Field)
    sim.cache.area_fraction .= field
end
function update_field!(sim::SurfaceStub, ::Val{:surface_temperature}, field::Fields.Field)
    sim.cache.T_sfc .= field
end
function update_field!(sim::SurfaceStub, ::Val{:air_density}, field)
    parent(sim.cache.ρ_sfc) .= parent(field)
end

"""
    name(::ComponentModelSimulation)

Returns simulation name, if defined, or `Unnamed` if not.
"""
name(::ComponentModelSimulation) = @error("undefined name for name(sim)")
name(::SurfaceStub) = "SurfaceStub"

end # module
