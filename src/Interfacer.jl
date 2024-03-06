"""
    Interfacer

This modules contains abstract types, interface templates and model stubs for coupling component models.
"""
module Interfacer
import Thermodynamics as TD

using ClimaCore: Fields
export CoupledSimulation,
    float_type,
    ComponentModelSimulation,
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
    AbstractSimulation

An abstract super-type representing a simulation.
"""
abstract type AbstractSimulation{FT} end

"""
    CoupledSimulation
Stores information needed to run a simulation with the coupler.
"""
struct CoupledSimulation{
    FT <: Real,
    X,
    D,
    B,
    FV,
    P,
    E,
    TS,
    TI <: Real,
    DTI <: Real,
    NTSM <: NamedTuple,
    NTMS <: NamedTuple,
    NTM <: NamedTuple,
    TD <: Tuple,
    NTC <: NamedTuple,
    NTP <: NamedTuple,
}
    comms_ctx::X
    dates::D
    boundary_space::B
    fields::FV
    parsed_args::P
    conservation_checks::E
    tspan::TS
    t::TI
    Δt_cpl::DTI
    surface_fractions::NTSM
    model_sims::NTMS
    mode::NTM
    diagnostics::TD
    callbacks::NTC
    dirs::NTP
end

CoupledSimulation{FT}(args...) where {FT} = CoupledSimulation{FT, typeof.(args[1:end])...}(args...)

"""
    float_type(::CoupledSimulation)

Return the floating point type backing `T`: `T` can either be an object or a type.
"""
float_type(::CoupledSimulation{FT}) where {FT} = FT

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
    get_field(sim::AtmosModelSimulation, val::Val)

A getter function that should not allocate. Here we implement a default that
will raise an error if `get_field` isn't defined for all required fields of
an atmosphere component model.
"""
get_field(
    sim::AtmosModelSimulation,
    val::Union{
        Val{:air_density},
        Val{:air_temperature},
        Val{:energy},
        Val{:height_int},
        Val{:height_sfc},
        Val{:liquid_precipitation},
        Val{:radiative_energy_flux_sfc},
        Val{:radiative_energy_flux_toa},
        Val{:snow_precipitation},
        Val{:turblent_energy_flux},
        Val{:turbulent_moisture_flux},
        Val{:thermo_state_int},
        Val{:uv_int},
        Val{:water},
    },
) = get_field_error(sim, val)

"""
    get_field(sim::SurfaceModelSimulation, val::Val)

A getter function that should not allocate. Here we implement a default that
will raise an error if `get_field` isn't defined for all required fields of
a surface component model.
"""
get_field(
    sim::SurfaceModelSimulation,
    val::Union{
        Val{:air_density},
        Val{:area_fraction},
        Val{:beta},
        Val{:roughness_buoyancy},
        Val{:roughness_momentum},
        Val{:surface_albedo},
        Val{:surface_humidity},
        Val{:surface_temperature},
    },
) = get_field_error(sim, val)

"""
    get_field(sim::ComponentModelSimulation, val::Val)

Generic fallback for `get_field` that raises an error.
"""
get_field(sim::ComponentModelSimulation, val::Val) = get_field_error(sim, val)

get_field_error(sim, val::Val{X}) where {X} = error("undefined field `$X` for " * name(sim))


"""
    SurfaceStub

On object containing simulation-like info, used as a stub or for prescribed data.
"""
struct SurfaceStub{I} <: SurfaceModelSimulation
    cache::I
end

"""
    stub_init(cache)

Initialization function for SurfaceStub simulation type.
"""
stub_init(cache) = SurfaceStub(cache)

"""
    get_field(::SurfaceStub, ::Val)

A getter function, that should not allocate. If undefined, it returns a descriptive error.
"""
get_field(sim::SurfaceStub, ::Val{:air_density}) = sim.cache.ρ_sfc
get_field(sim::SurfaceStub, ::Val{:area_fraction}) = sim.cache.area_fraction
get_field(sim::SurfaceStub, ::Val{:beta}) = sim.cache.beta
get_field(sim::SurfaceStub, ::Val{:energy}) = nothing
get_field(sim::SurfaceStub, ::Val{:roughness_buoyancy}) = sim.cache.z0b
get_field(sim::SurfaceStub, ::Val{:roughness_momentum}) = sim.cache.z0m
get_field(sim::SurfaceStub, ::Val{:surface_albedo}) = sim.cache.α
get_field(sim::SurfaceStub, ::Val{:surface_humidity}) =
    TD.q_vap_saturation_generic.(sim.cache.thermo_params, sim.cache.T_sfc, sim.cache.ρ_sfc, sim.cache.phase)
get_field(sim::SurfaceStub, ::Val{:surface_temperature}) = sim.cache.T_sfc
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
    update_field!(::AtmosModelSimulation, ::Val, _...)

Default functions for updating fields at each timestep in an atmosphere
component model simulation. This should be extended by component models.
If it isn't extended, the field won't be updated and a warning will be raised.
"""
update_field!(
    sim::AtmosModelSimulation,
    val::Union{Val{:co2}, Val{:surface_albedo}, Val{:surface_temperature}, Val{:turbulent_fluxes}},
    _,
) = update_field_warning(sim, val)

"""
    update_field!(::SurfaceModelSimulation, ::Val, _...)

Default functions for updating fields at each timestep in an atmosphere
component model simulation. This should be extended by component models.
If it isn't extended, the field won't be updated and a warning will be raised.
"""
update_field!(
    sim::SurfaceModelSimulation,
    val::Union{
        Val{:air_density},
        Val{:area_fraction},
        Val{:liquid_precipitation},
        Val{:radiative_energy_flux_sfc},
        Val{:snow_precipitation},
        Val{:turbulent_energy_flux},
        Val{:turbulent_moisture_flux},
    },
    _,
) = update_field_warning(sim, val)

update_field_warning(sim, val::Val{X}) where {X} =
    @warn("`update_field!` is not extended for the `$X` field of " * name(sim) * ": skipping update.", maxlog = 1)

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
name(::ComponentModelSimulation) = "Unnamed"
name(::SurfaceStub) = "SurfaceStub"

end # module
