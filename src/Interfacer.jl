"""
    Interfacer

This modules contains abstract types, interface templates and model stubs for coupling component models.
"""
module Interfacer

import SciMLBase
import ClimaCore as CC
import Thermodynamics as TD
import SciMLBase: step!, reinit! # explicitly import to extend these functions

export CoupledSimulation,
    float_type,
    ComponentModelSimulation,
    AtmosModelSimulation,
    SurfaceModelSimulation,
    SeaIceModelSimulation,
    LandModelSimulation,
    OceanModelSimulation,
    name,
    get_field,
    update_field!,
    SurfaceStub,
    step!,
    reinit!,
    AbstractSlabplanetSimulationMode,
    AMIPMode,
    SlabplanetMode,
    SlabplanetAquaMode,
    SlabplanetTerraMode,
    SlabplanetEisenmanMode


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
    E,
    TS,
    DTI <: Real,
    NTMS <: NamedTuple,
    NTM <: NamedTuple,
    NTC <: NamedTuple,
    NTP <: NamedTuple,
    TF,
    TP,
    DH,
}
    comms_ctx::X
    dates::D
    boundary_space::B
    fields::FV
    conservation_checks::E
    tspan::TS
    Δt_cpl::DTI
    model_sims::NTMS
    mode::NTM
    callbacks::NTC
    dirs::NTP
    turbulent_fluxes::TF
    thermo_params::TP
    amip_diags_handler::DH
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
        Val{:surface_direct_albedo},
        Val{:surface_diffuse_albedo},
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
    update_field!(::AtmosModelSimulation, ::Val, _...)

Default functions for updating fields at each timestep in an atmosphere
component model simulation. This should be extended by component models.
If it isn't extended, the field won't be updated and a warning will be raised.
"""
update_field!(
    sim::AtmosModelSimulation,
    val::Union{
        Val{:co2},
        Val{:surface_direct_albedo},
        Val{:surface_diffuse_albedo},
        Val{:surface_temperature},
        Val{:turbulent_fluxes},
    },
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
        Val{:surface_direct_albedo},
        Val{:surface_diffuse_albedo},
    },
    _,
) = update_field_warning(sim, val)

update_field_warning(sim, val::Val{X}) where {X} =
    @warn("`update_field!` is not extended for the `$X` field of " * name(sim) * ": skipping update.", maxlog = 1)

"""
    name(::ComponentModelSimulation)

Returns simulation name, if defined, or `Unnamed` if not.
"""
name(::ComponentModelSimulation) = "Unnamed"

"""
    step!(sim::ComponentModelSimulation, t)

A function to update the simulation in-place with values calculate for time `t`.
For the models we currently have implemented, this is a simple wrapper around
the `step!` function implemented in SciMLBase.jl.

This must be extended for all component models - otherwise this default
function will be called and an error will be raised.
"""
step!(sim::ComponentModelSimulation, t) = error("undefined step! for " * name(sim))

"""
    reinit!(sim::ComponentModelSimulation)

A function to restart a simulation after solving of the simulation has been
paused or interrupted. Like `step!`, this is currently a simple wrapper
around the `reinit!` function of SciMLBase.jl.

This must be extended for all component models - otherwise this default
function will be called and an error will be raised.
"""
reinit!(sim::ComponentModelSimulation) = error("undefined reinit! for " * name(sim))


# Include file containing the surface stub simulation type.
include("surface_stub.jl")

"""
    AbstractModeType

An abstract type representing a simulation mode.
"""
abstract type AbstractSimulationMode end

"""
    AbstractSlabplanetSimulationMode

An abstract type representing a simulation mode for slabplanet models.
"""
abstract type AbstractSlabplanetSimulationMode <: AbstractSimulationMode end

struct AMIPMode <: AbstractSimulationMode end
struct SlabplanetMode <: AbstractSlabplanetSimulationMode end
struct SlabplanetAquaMode <: AbstractSlabplanetSimulationMode end
struct SlabplanetTerraMode <: AbstractSlabplanetSimulationMode end
struct SlabplanetEisenmanMode <: AbstractSlabplanetSimulationMode end

end # module
