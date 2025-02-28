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
    AbstractSurfaceStub,
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
    DTI,
    NTMS <: NamedTuple,
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
    callbacks::NTC
    dirs::NTP
    turbulent_fluxes::TF
    thermo_params::TP
    diags_handler::DH
end

CoupledSimulation{FT}(args...) where {FT} = CoupledSimulation{FT, typeof.(args[1:end])...}(args...)

"""
    float_type(::CoupledSimulation)

Return the floating point type backing `T`: `T` can either be an object or a type.
"""
float_type(::CoupledSimulation{FT}) where {FT} = FT

"""
    default_coupler_fields()

Return a list of default coupler fields needed to run a simulation.
"""
default_coupler_fields() = [
    # fields needed for flux calculations and exchange
    :z0m_sfc,
    :z0b_sfc,
    :beta,
    :F_turb_energy,
    :F_turb_moisture,
    :F_turb_ρτxz,
    :F_turb_ρτyz,
    # fields used for temporary storage during calculations
    :temp1,
    :temp2,
]

"""
    init_coupler_fields(FT, coupler_field_names, boundary_space)

Allocate a Field of NamedTuples on the provided boundary space to store
the provided coupler fields.
"""
function init_coupler_fields(FT, coupler_field_names, boundary_space)
    # First remove any duplicate field names
    unique!(coupler_field_names)

    key_types = (coupler_field_names...,)
    val_types = Tuple{(FT for _ in 1:length(coupler_field_names))...}
    nt_type = NamedTuple{key_types, val_types}
    coupler_fields = zeros(nt_type, boundary_space)
    return coupler_fields
end

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
        Val{:air_pressure},
        Val{:air_temperature},
        Val{:cos_zenith_angle},
        Val{:co2},
        Val{:diffuse_fraction},
        Val{:height_int},
        Val{:height_sfc},
        Val{:humidity},
        Val{:liquid_precipitation},
        Val{:LW_d},
        Val{:radiative_energy_flux_sfc},
        Val{:radiative_energy_flux_toa},
        Val{:snow_precipitation},
        Val{:SW_d},
        Val{:turblent_energy_flux},
        Val{:turbulent_moisture_flux},
        Val{:thermo_state_int},
        Val{:uv_int},
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
        Val{:area_fraction},
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

# Set default values for fields that are not defined in all component models
get_field(::ComponentModelSimulation, ::Val{:energy}) = nothing
get_field(::ComponentModelSimulation, ::Val{:water}) = nothing
get_field(sim::SurfaceModelSimulation, ::Val{:beta}) = convert(eltype(sim.integrator.u), 1.0)
get_field(sim::SurfaceModelSimulation, ::Val{:emissivity}) = convert(eltype(sim.integrator.u), 1.0)
get_field(sim::SurfaceModelSimulation, ::Val{:height_disp}) = convert(eltype(sim.integrator.u), 0.0)

"""
    update_field!(::AtmosModelSimulation, ::Val, _...)

Default functions for updating fields at each timestep in an atmosphere
component model simulation. This should be extended by component models.
If it isn't extended, the field won't be updated and a warning will be raised.
"""
update_field!(
    sim::AtmosModelSimulation,
    val::Union{
        Val{:emissivity},
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
    add_coupler_fields!(coupler_fields, sim::ComponentModelSimulation, fields)

A function to add fields to the set of coupler fields. This should be extended
by component models that require coupler fields beyond the defaults.

If this function isn't extended, no additional fields will be added.
"""
add_coupler_fields!(coupler_fields, sim::ComponentModelSimulation) = nothing

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
    AbstractSimulationMode

An abstract type representing a simulation mode.
"""
abstract type AbstractSimulationMode end

"""
    AbstractSlabplanetSimulationMode

An abstract type representing a simulation mode for slabplanet models. Slabplanet simulations
are more idealized than the AMIP configuration, but provide valuable insight about
conservation and individual model behavior.
"""
abstract type AbstractSlabplanetSimulationMode <: AbstractSimulationMode end

"""
    AMIPMode

An abstract type representing the AMIP simulation mode. AMIP is currently the most complex
configuration of the ClimaEarth model. It runs a ClimaAtmos.jl atmosphere model,
ClimaLand.jl bucket land model, a prescribed ocean model, and a simple thermal sea ice model.
"""
abstract type AMIPMode <: AbstractSimulationMode end

"""
    SlabplanetMode

An abstract type represeting the slabplanet simulation mode with a ClimaAtmos.jl atmosphere model,
a ClimaLand.jl bucket land model, a thermal slab ocean model, and no sea ice model. Instead
of using a sea ice model, the ocean is evaluated in areas that would be covered in ice.
"""
abstract type SlabplanetMode <: AbstractSlabplanetSimulationMode end

"""
    SlabplanetAquaMode

An abstract type representing the slabplanet simulation mode with a ClimaAtmos.jl atmosphere model,
and only once surface model, a thermal slab ocean model, which is evaluated over the entire
surface. There are no land or sea ice models.
"""
abstract type SlabplanetAquaMode <: AbstractSlabplanetSimulationMode end

"""
    SlabplanetTerraMode

An abstract type representing the slabplanet simulation mode with a ClimaAtmos.jl atmosphere model,
and only once surface model, a ClimaLand.jl bucket land model, which is evaluated over the
entire surface. There are no ocean or sea ice models.
"""
abstract type SlabplanetTerraMode <: AbstractSlabplanetSimulationMode end

"""
    SlabplanetEisenmanMode

An abstract type representing the slabplanet simulation mode with a ClimaAtmos.jl atmosphere model,
a ClimaLand.jl bucket land model, and Eisenman sea ice model. The ocean model
is included in the Eisenman sea ice model.
"""
abstract type SlabplanetEisenmanMode <: AbstractSlabplanetSimulationMode end

end # module
