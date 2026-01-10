"""
    Interfacer

This modules contains abstract types, interface templates and model stubs for coupling component models.
"""
module Interfacer

import SciMLBase
import ClimaComms
import ClimaCore as CC
import Dates
import Thermodynamics as TD
import SciMLBase: step!
import ClimaUtilities.TimeManager: ITime, date

export CoupledSimulation,
    LandSimulation,
    OceanSimulation,
    SeaIceSimulation,
    AtmosSimulation,
    AbstractComponentSimulation,
    AbstractAtmosSimulation,
    AbstractSurfaceSimulation,
    AbstractSeaIceSimulation,
    AbstractLandSimulation,
    AbstractOceanSimulation,
    get_field,
    update_field!,
    AbstractSurfaceStub,
    SurfaceStub,
    step!,
    set_cache!,
    remap,
    remap!,
    AbstractSlabplanetSimulationMode,
    AMIPMode,
    CMIPMode,
    SlabplanetMode,
    SlabplanetAquaMode,
    SlabplanetTerraMode

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
    D,
    FV,
    E,
    TS,
    DTI,
    TT,
    CTT,
    NTMS <: NamedTuple,
    CALLBACKS,
    NTP <: NamedTuple,
    TP,
    DH,
    SC <: Bool,
}
    start_date::D
    fields::FV
    conservation_checks::E
    tspan::TS
    Δt_cpl::DTI
    t::TT
    prev_checkpoint_t::CTT
    model_sims::NTMS
    callbacks::CALLBACKS
    dir_paths::NTP
    thermo_params::TP
    diags_handler::DH
    save_cache::SC
end

CoupledSimulation{FT}(args...) where {FT} = CoupledSimulation{FT, typeof.(args)...}(args...)

function Base.show(io::IO, sim::CoupledSimulation)
    device_type = nameof(typeof(ClimaComms.device(sim)))
    return print(
        io,
        "Coupled Simulation\n",
        "├── Running on: $(device_type)\n",
        "├── Output folder: $(sim.dir_paths.output_dir_root)\n",
        "└── Current date: $(current_date(sim))\n",
    )
end

"""
    current_date(cs::CoupledSimulation)

Return the model date at the current timestep.
# Arguments
- `cs`: [CoupledSimulation] containing info about the simulation
"""
current_date(cs::CoupledSimulation) =
    cs.t[] isa ITime ? date(cs.t[]) : cs.start_date[] + Dates.Second(cs.t[])

"""
    default_coupler_fields()

Return a list of default coupler fields needed to run a simulation.
"""
default_coupler_fields() = [
    # fields used to compute turbulent fluxes
    :T_atmos,
    :q_tot_atmos,
    :q_liq_atmos,
    :q_ice_atmos,
    :ρ_atmos,
    :height_int,
    :height_sfc,
    :height_delta,
    :u_int,
    :v_int,
    # fields used for flux exchange
    :F_lh,
    :F_sh,
    :F_turb_moisture,
    :F_turb_ρτxz,
    :F_turb_ρτyz,
    :SW_d,
    :LW_d,
    # fields used to compute radiation in the atmosphere
    :emissivity,
    :T_sfc,
    # fields used for water fluxes and tracking water conservation
    :P_liq,
    :P_snow,
    # fields used for temporary storage during calculations
    :scalar_temp1,
    :scalar_temp2,
    :scalar_temp3,
    :scalar_temp4,
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
    AbstractComponentSimulation

An abstract type encompassing all component model (and model stub) simulations.
"""
abstract type AbstractComponentSimulation end

"""
    AbstractAtmosSimulation

An abstract type for an atmospheric model simulation.
"""
abstract type AbstractAtmosSimulation <: AbstractComponentSimulation end

"""
    AbstractSurfaceSimulation

An abstract type for surface model simulations.
"""
abstract type AbstractSurfaceSimulation <: AbstractComponentSimulation end

abstract type AbstractSeaIceSimulation <: AbstractSurfaceSimulation end
abstract type AbstractLandSimulation <: AbstractSurfaceSimulation end
abstract type AbstractOceanSimulation <: AbstractSurfaceSimulation end

"""
    AbstractImplicitFluxSimulation

An abstract type for surface model simulations that compute fluxes implicitly,
rather than explicitly. At the moment, this means the fluxes are computed in the
component model's `step!` function, rather than in the coupler's `compute_surface_fluxes!`
function.

Currently, the only implicit flux simulation is the integrated land model.
"""
abstract type AbstractImplicitFluxSimulation <: AbstractSurfaceSimulation end

# Simulation objects tend to be very big, so it is best to make sure they are not printed in the REPL
function Base.show(io::IO, @nospecialize(sim::AbstractComponentSimulation))
    return println(io, "$(nameof(sim)) without a specialized `Base.show` method")
end

"""
    get_field(sim::AbstractAtmosSimulation, val::Val)

A getter function that should not allocate. Here we implement a default that
will raise an error if `get_field` isn't defined for all required fields of
an atmosphere component model.
"""
get_field(
    sim::AbstractAtmosSimulation,
    val::Union{
        Val{:height_int},
        Val{:height_sfc},
        Val{:liquid_precipitation},
        Val{:SW_d},
        Val{:LW_d},
        Val{:snow_precipitation},
        Val{:turblent_energy_flux},
        Val{:turbulent_moisture_flux},
        Val{:u_int},
        Val{:v_int},
    },
) = get_field_error(sim, val)

"""
    get_field(sim::AbstractSurfaceSimulation, val::Val)

A getter function that should not allocate. Here we implement a default that
will raise an error if `get_field` isn't defined for all required fields of
a surface component model.
"""
get_field(
    sim::AbstractSurfaceSimulation,
    val::Union{
        Val{:area_fraction},
        Val{:roughness_buoyancy},
        Val{:roughness_momentum},
        Val{:surface_direct_albedo},
        Val{:surface_diffuse_albedo},
        Val{:surface_temperature},
    },
) = get_field_error(sim, val)

# Sea ice models need to provide ice concentration
get_field(sim::AbstractSeaIceSimulation, val::Val{:ice_concentration}) =
    get_field_error(sim, val)

"""
    get_field(sim::AbstractComponentSimulation, val::Val)

Generic fallback for `get_field` that raises an error.
"""
get_field(sim::AbstractComponentSimulation, val::Val) = get_field_error(sim, val)

get_field_error(sim, val::Val{X}) where {X} =
    error("undefined field `$X` for $(nameof(sim))")

# Set default values for fields that are not defined in all component models
get_field(::AbstractComponentSimulation, ::Val{:energy}) = nothing
get_field(::AbstractComponentSimulation, ::Val{:water}) = nothing
get_field(sim::AbstractSurfaceSimulation, ::Val{:emissivity}) =
    convert(eltype(sim.integrator.u), 1.0)
get_field(sim::AbstractSurfaceSimulation, ::Val{:height_disp}) =
    convert(eltype(sim.integrator.u), 0.0)


"""
    get_field(target_space, sim, quantity)

Return `quantity` in `sim` remapped onto the `target_space`

This is equivalent to calling `get_field`, and then `remap`.
"""
function get_field(target_space, sim, quantity)
    return remap(target_space, get_field(sim, quantity))
end

"""
    get_field!(target_field, sim, quantity)

Remap `quantity` in `sim` remapped onto the `target_field`.
"""
function get_field!(target_field, sim, quantity)
    remap!(target_field, get_field(sim, quantity))
    return nothing
end

"""
    update_field!(::AbstractAtmosSimulation, ::Val, _...)

Default functions for updating fields at each timestep in an atmosphere
component model simulation. This should be extended by component models.
If it isn't extended, the field won't be updated and a warning will be raised.
"""
update_field!(
    sim::AbstractAtmosSimulation,
    val::Union{
        Val{:emissivity},
        Val{:surface_direct_albedo},
        Val{:surface_diffuse_albedo},
        Val{:surface_temperature},
    },
    _,
) = update_field_warning(sim, val)

"""
    update_field!(::AbstractSurfaceSimulation, ::Val, _...)

Default functions for updating fields at each timestep in an atmosphere
component model simulation. This should be extended by component models.
If it isn't extended, the field won't be updated and a warning will be raised.
"""
update_field!(
    sim::AbstractSurfaceSimulation,
    val::Union{
        Val{:area_fraction},
        Val{:liquid_precipitation},
        Val{:SW_d},
        Val{:LW_d},
        Val{:snow_precipitation},
        Val{:turbulent_energy_flux},
        Val{:turbulent_moisture_flux},
    },
    _,
) = update_field_warning(sim, val)

update_field_warning(sim, val::Val{X}) where {X} = @warn(
    "`update_field!` is not extended for the `$X` field of $(nameof(sim)): skipping update.",
    maxlog = 1
)


"""
    add_coupler_fields!(coupler_fields, sim::AbstractComponentSimulation, fields)

A function to add fields to the set of coupler fields. This should be extended
by component models that require coupler fields beyond the defaults.

If this function isn't extended, no additional fields will be added.
"""
add_coupler_fields!(coupler_fields, sim::AbstractComponentSimulation) = nothing

"""
    Base.nameof(::AbstractComponentSimulation)

Return the simulation name, if defined, or the type name if not.
"""
Base.nameof(sim::AbstractComponentSimulation) = string(nameof(typeof(sim)))

"""
    step!(sim::AbstractComponentSimulation, t)

A function to update the simulation in-place with values calculate for time `t`.
For the models we currently have implemented, this is a simple wrapper around
the `step!` function implemented in SciMLBase.jl.

This must be extended for all component models - otherwise this default
function will be called and an error will be raised.
"""
step!(sim::AbstractComponentSimulation, t) = error("undefined step! for $(nameof(sim))")

"""
    close_output_writers(sim::AbstractComponentSimulation)

A function to close all output writers associated with the given
component model, at the end of a simulation.

This should be extended for any component model that uses
an output writer.
"""
close_output_writers(sim::AbstractComponentSimulation) = nothing

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

An abstract type representing the AMIP simulation mode. It runs a ClimaAtmos.jl atmosphere model,
ClimaLand.jl bucket land model, a prescribed ocean model, and a simple thermal sea ice model.
"""
abstract type AMIPMode <: AbstractSimulationMode end


"""
    SubseasonalMode

An abstract type representing the subseasonal simulation mode. This mode is similar to AMIP
but uses different data sources and initialization pathways tailored for subseasonal runs.

Inputs are ERA5-derived netcdfs with initial conditions produced by `https://github.com/CliMA/WeatherQuest`. Given `start_date` (YYYYMMDD)
and directory `era5_initial_condition_dir`, filenames containing the initial conditions are inferred as:
- `sst_processed_YYYYMMDD_0000.nc` (variable `SST`)
- `sic_processed_YYYYMMDD_0000.nc` (variable `SEAICE`)
- Land IC (integrated land): `era5_land_processed_YYYYMMDD_0000.nc`, with fields
  - `skt` (K), `tsn` (K),`swe` (m), `swvl` (m^3/m^3), `si` (m^3/m^3), `sie` (J/m^3), `stl` (K)
- Land IC (bucket land): `era5_bucket_processed_YYYYMMDD_0000.nc`, with fields
  - `W` (m), `Ws` (m), `S` (m), `T` (K), `tsn` (K), `skt` (K); dims `(lat, lon)`
- Albedo (optional, when `bucket_albedo_type: "era5"`): `albedo_processed_YYYYMMDD_0000.nc`, with fields
  - `sw_alb_clr` (clear-sky surface albedo, fraction 0-1); dims `(time, lat, lon)`
and are used to initialize the coupler components.

"""
abstract type SubseasonalMode <: AbstractSimulationMode end

"""
    CMIPMode

An abstract type representing the CMIP simulation mode. CMIP is currently the most complex
configuration of the ClimaEarth model. It runs a ClimaAtmos.jl atmosphere model,
ClimaLand.jl bucket land model, a ClimaOcean ocean model, and a simple thermal sea ice model.
"""
abstract type CMIPMode <: AbstractSimulationMode end

"""
    SlabplanetMode

An abstract type representing the slabplanet simulation mode with a ClimaAtmos.jl atmosphere model,
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
    remap(target_space, source_field)

Remap the given `source_field` onto the `target_space`. Note that if the field is already
on the target space or a compatible one (e.g. another instance of the same space),
it is returned unchanged. Users should use caution in modifying the returned field
in this case.

This is a convenience wrapper around `remap!` that allocates the output field.

Non-ClimaCore fields should provide a method to this function.
"""
function remap end

function remap(target_space::CC.Spaces.AbstractSpace, source_field::CC.Fields.Field)
    source_space = axes(source_field)

    # Check if the source and target spaces are compatible
    spaces_are_compatible =
        source_space == target_space ||
        CC.Spaces.issubspace(source_space, target_space) ||
        CC.Spaces.issubspace(target_space, source_space)

    # TODO: Handle remapping of Vectors correctly
    if hasproperty(source_field, :components)
        @assert length(source_field.components) == 1 "Can only work with simple vectors"
        source_field = source_field.components.data.:1
    end

    # If the spaces are the same or one is a subspace of the other, we can just return the input field
    spaces_are_compatible && return source_field

    # Allocate target field and call remap!
    target_field = CC.Fields.zeros(target_space)
    remap!(target_field, source_field)
    return target_field
end

function remap(target_space::CC.Spaces.AbstractSpace, source_field::Number)
    # Allocate target field and call remap!
    target_field = CC.Fields.zeros(target_space)
    remap!(target_field, source_field)
    return target_field
end

"""
    remap!(target_field, source_field)

Remap the given `source_field` onto the `target_field`. This is the core non-allocating
implementation.

Non-ClimaCore fields should provide a method to this function.

Note that this method has a lot of allocations and is not efficient.
"""
function remap! end

function remap!(target_field::CC.Fields.Field, source_field::CC.Fields.Field)
    source_space = axes(source_field)
    target_space = axes(target_field)
    comms_ctx = ClimaComms.context(source_space)

    # Check if the source and target spaces are compatible
    spaces_are_compatible =
        source_space == target_space ||
        CC.Spaces.issubspace(source_space, target_space) ||
        CC.Spaces.issubspace(target_space, source_space)

    # TODO: Handle remapping of Vectors correctly
    if hasproperty(source_field, :components)
        @assert length(source_field.components) == 1 "Can only work with simple vectors"
        source_field = source_field.components.data.:1
    end

    # If the spaces are the same or one is a subspace of the other, we can just copy
    if spaces_are_compatible
        target_field .= source_field
        return nothing
    end

    # Get vector of LatLongPoints for the target space to get the hcoords
    # Copy target coordinates to CPU if they are on GPU
    coords = CC.to_cpu(CC.Fields.coordinate_field(target_space))
    lats = CC.Fields.field2array(coords.lat)
    lons = CC.Fields.field2array(coords.long)
    hcoords = CC.Geometry.LatLongPoint.(lats, lons)

    # Remap the field, using MPI if applicable
    if comms_ctx isa ClimaComms.SingletonCommsContext
        # Remap source field to target space as an array
        remapped_array = CC.Remapping.interpolate(source_field, hcoords, [])

        # Write directly to target field's underlying array to avoid temporary field allocation
        CC.Fields.field2array(target_field) .= remapped_array
    else
        # Gather then broadcast the global hcoords and offsets
        offset = [length(hcoords)]
        all_hcoords = ClimaComms.bcast(comms_ctx, ClimaComms.gather(comms_ctx, hcoords))
        all_offsets = ClimaComms.bcast(comms_ctx, ClimaComms.gather(comms_ctx, offset))

        # Interpolate on root and broadcast to all processes
        remapper = CC.Remapping.Remapper(source_space; target_hcoords = all_hcoords)
        remapped_array =
            ClimaComms.bcast(comms_ctx, CC.Remapping.interpolate(remapper, source_field))

        my_ending_offset = sum(all_offsets[1:ClimaComms.mypid(comms_ctx)])
        my_starting_offset = my_ending_offset - offset[]

        # Write directly to target field's underlying array to avoid temporary field allocation
        CC.Fields.field2array(target_field) .=
            remapped_array[(1 + my_starting_offset):my_ending_offset]
    end
    return nothing
end

function remap!(target_field::CC.Fields.Field, source::Number)
    fill!(target_field, source)
    return nothing
end


"""
    set_cache!(sim::AbstractComponentSimulation, csf)

Perform any initialization of the component model cache that must be done
after the initial exchange.
This is not required to be extended, but may be necessary for some models.
"""
set_cache!(sim::AbstractComponentSimulation, csf) = nothing

"""
    LandSimulation(::Type{FT}, ::Val{model_type}; kwargs...)

Generic constructor for land model simulations. Dispatches to specific implementations
based on `model_type` (`:bucket`, `:integrated`, or `:nothing`).

If `model_type` is `:nothing`, returns `nothing`.
"""
LandSimulation(::Type{FT}, ::Val{:nothing}; kwargs...) where {FT} = nothing
function LandSimulation(::Type{FT}, ::Val{model_type}; kwargs...) where {FT, model_type}
    error(
        "Unknown land model type: $model_type. Valid options are: :bucket, :integrated, or :nothing",
    )
end

"""
    OceanSimulation(::Type{FT}, ::Val{model_type}; kwargs...)

Generic constructor for ocean model simulations. Dispatches to specific implementations
based on `model_type` (`:oceananigans`, `:slab`, `:prescribed`, or `:nothing`).

If `model_type` is `:nothing`, returns `nothing`.
Some ocean models (like `:oceananigans`) don't require FT as the first argument.
"""
OceanSimulation(::Type{FT}, ::Val{:nothing}; kwargs...) where {FT} = nothing
function OceanSimulation(::Type{FT}, ::Val{model_type}; kwargs...) where {FT, model_type}
    error(
        "Unknown ocean model type: $model_type. Valid options are: :oceananigans, :slab, :prescribed, or :nothing",
    )
end

"""
    SeaIceSimulation(::Type{FT}, ::Val{model_type}; kwargs...)

Generic constructor for sea ice model simulations. Dispatches to specific implementations
based on `model_type` (`:clima_seaice`, `:prescribed`, or `:nothing`).

If `model_type` is `:nothing`, returns `nothing`.
FT is passed to all sea ice models, though some may ignore it.
"""
SeaIceSimulation(::Type{FT}, ::Val{:nothing}; kwargs...) where {FT} = nothing
function SeaIceSimulation(::Type{FT}, ::Val{model_type}; kwargs...) where {FT, model_type}
    error(
        "Unknown sea ice model type: $model_type. Valid options are: :clima_seaice, :prescribed, or :nothing",
    )
end

"""
    AtmosSimulation(::Val{model_type}; kwargs...)

Generic constructor for atmosphere model simulations. Dispatches to specific implementations
based on `model_type` (`:climaatmos`).

Note that the atmosphere model cannot be nothing.
"""
AtmosSimulation(::Val{model_type}; kwargs...) where {model_type} =
    error("Unknown atmosphere model type: $model_type. Valid options are: :climaatmos")

"""
    boundary_space(sim::CoupledSimulation)

Return the `ClimaCore.Field` over which the exchange fields are defined.
"""
function boundary_space(sim::CoupledSimulation)
    return axes(sim.fields)
end

"""
    ClimaComms.context(sim::CoupledSimulation)

Return the `ClimaComms.context` associated to the simulation.
"""
function ClimaComms.context(sim::CoupledSimulation)
    return ClimaComms.context(sim.fields)
end

"""
    ClimaComms.device(sim::CoupledSimulation)

Return the `ClimaComms.device` associated to the simulation.
"""
function ClimaComms.device(sim::CoupledSimulation)
    return ClimaComms.device(sim.fields)
end

"""
    get_atmos_height_delta(height_int, height_sfc)

Return a Field of the height delta between the atmosphere bottom cell center
and bottom face, defined on the boundary space.
This is used to compute turbulent fluxes.

Since the atmospheric height is defined on centers, we need to copy the values onto
the boundary space to be able to subtract the surface elevation.
This pattern is not reliable and should not be reused.

Note this function allocates a new field, and the atmosphere heights won't change
during a simulation, so it should only be called at initialization.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.

# Returns
- [CC.Fields.Field] defined on the boundary space containing the height delta.
"""
function get_atmos_height_delta(height_int, height_sfc)
    return CC.Fields.Field(
        CC.Fields.field_values(height_int),
        axes(height_sfc), # boundary space
    ) .- height_sfc
end

end # module
