"""
    Interfacer

This modules contains abstract types, interface templates and model stubs for coupling component models.
"""
module Interfacer

import SciMLBase
import ClimaComms
import ClimaCore as CC
import NVTX
import Dates
import Thermodynamics as TD
import SciMLBase: step!
import ClimaUtilities.TimeManager: ITime, date

export CoupledSimulation,
    ComponentModelSimulation,
    AtmosModelSimulation,
    SurfaceModelSimulation,
    SeaIceModelSimulation,
    LandModelSimulation,
    OceanModelSimulation,
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
    :q_atmos,
    :ρ_atmos,
    :z_int,
    :z_sfc,
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

# Simulation objects tend to be very big, so it is best to make sure they are not printed in the REPL
function Base.show(io::IO, @nospecialize(sim::ComponentModelSimulation))
    return println(io, "$(nameof(sim)) without a specialized `Base.show` method")
end

"""
    get_field(sim::AtmosModelSimulation, val::Val)

A getter function that should not allocate. Here we implement a default that
will raise an error if `get_field` isn't defined for all required fields of
an atmosphere component model.
"""
get_field(
    sim::AtmosModelSimulation,
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
        Val{:surface_temperature},
    },
) = get_field_error(sim, val)

# Sea ice models need to provide ice concentration
get_field(sim::SeaIceModelSimulation, val::Val{:ice_concentration}) =
    get_field_error(sim, val)

"""
    get_field(sim::ComponentModelSimulation, val::Val)

Generic fallback for `get_field` that raises an error.
"""
get_field(sim::ComponentModelSimulation, val::Val) = get_field_error(sim, val)

get_field_error(sim, val::Val{X}) where {X} =
    error("undefined field `$X` for $(nameof(sim))")

# Set default values for fields that are not defined in all component models
get_field(::ComponentModelSimulation, ::Val{:energy}) = nothing
get_field(::ComponentModelSimulation, ::Val{:water}) = nothing
get_field(sim::SurfaceModelSimulation, ::Val{:beta}) =
    convert(eltype(sim.integrator.u), 1.0)
get_field(sim::SurfaceModelSimulation, ::Val{:emissivity}) =
    convert(eltype(sim.integrator.u), 1.0)
get_field(sim::SurfaceModelSimulation, ::Val{:height_disp}) =
    convert(eltype(sim.integrator.u), 0.0)


"""
    get_field(sim, what, target_space)

Return `quantity` in `sim` remapped onto the `target_space`

This is equivalent to calling `get_field`, and then `remap`.
"""
function get_field(sim, quantity, target_space)
    return remap(get_field(sim, quantity), target_space)
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
    add_coupler_fields!(coupler_fields, sim::ComponentModelSimulation, fields)

A function to add fields to the set of coupler fields. This should be extended
by component models that require coupler fields beyond the defaults.

If this function isn't extended, no additional fields will be added.
"""
add_coupler_fields!(coupler_fields, sim::ComponentModelSimulation) = nothing

"""
    Base.nameof(::ComponentModelSimulation)

Return the simulation name, if defined, or the type name if not.
"""
Base.nameof(sim::ComponentModelSimulation) = string(nameof(typeof(sim)))

"""
    step!(sim::ComponentModelSimulation, t)

A function to update the simulation in-place with values calculate for time `t`.
For the models we currently have implemented, this is a simple wrapper around
the `step!` function implemented in SciMLBase.jl.

This must be extended for all component models - otherwise this default
function will be called and an error will be raised.
"""
step!(sim::ComponentModelSimulation, t) = error("undefined step! for $(nameof(sim))")

"""
    close_output_writers(sim::ComponentModelSimulation)

A function to close all output writers associated with the given
component model, at the end of a simulation.

This should be extended for any component model that uses
an output writer.
"""
close_output_writers(sim::ComponentModelSimulation) = nothing

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
    remap(field, target_space)

Remap the given `field` onto the `target_space`. If the field is already
on the target space or a compatible one, it is returned unchanged.

Note that this method has a lot of allocations and is not efficient.

Non-ClimaCore fields should provide a method to this function.
"""
function remap end

# Pretty descriptions for debugging remaps
function describe_space(space)
    topo = try
        CC.Spaces.topology(space)
    catch
        nothing
    end
    mesh = isnothing(topo) ? nothing : (
        try
            topo.mesh
        catch
            nothing
        end
    )
    quad = try
        CC.Spaces.quadrature_style(space)
    catch
        nothing
    end
    deg = try
        isnothing(quad) ? nothing : CC.Spaces.Quadratures.polynomial_degree(quad)
    catch
        nothing
    end
    return (
        type = typeof(space),
        id = objectid(space),
        quad_degree = deg,
        topo_type = isnothing(topo) ? nothing : typeof(topo),
        mesh_type = isnothing(mesh) ? nothing : typeof(mesh),
    )
end

function describe_field(f)
    vals = CC.Fields.field_values(f)
    storage_t = typeof(vals)
    # Avoid hard dependency on CUDA
    on_gpu = occursin("CuArray", string(storage_t))
    sz = try
        size(vals)
    catch
        ()
    end
    return (eltype = eltype(vals), size = sz, on_gpu = on_gpu, storage_type = storage_t)
end

# A short, filtered callsite snapshot (skip Interfacer frames)
function short_callsite(; max_frames = 3)
    bt = Base.backtrace()
    frames = Base.stacktrace(bt)
    keep = filter(fr -> !occursin("Interfacer.jl", String(fr.file)), frames)
    n = min(length(keep), max_frames)
    return [
        (; file = String(keep[i].file), line = keep[i].line, func = keep[i].func) for
        i in 1:n
    ]
end

function remap(
    field::CC.Fields.Field,
    target_space::CC.Spaces.AbstractSpace;
    label::Union{Nothing, String, Symbol} = nothing,
)
    source_space = axes(field)
    comms_ctx = ClimaComms.context(source_space)

    # Check if the source and target spaces are compatible
    spaces_are_compatible =
        source_space == target_space ||
        CC.Spaces.issubspace(source_space, target_space) ||
        CC.Spaces.issubspace(target_space, source_space)
    if !spaces_are_compatible
        @info (
            "Remapping field from source space:\n$source_space\n" *
            "to incompatible target space:\n$target_space"
        )
        cs = short_callsite(max_frames = 5)
        @info "Interfacer.remap: incompatible spaces" label callsite = cs
    end

    # TODO: Handle remapping of Vectors correctly
    if hasproperty(field, :components)
        println("Field has components")
        @assert length(field.components) == 1 "Can only work with simple vectors"
        field = field.components.data.:1
    end

    # If the spaces are the same or one is a subspace of the other, we can just return the input field
    spaces_are_compatible && return field

    # Get vector of LatLongPoints for the target space to get the hcoords
    # Copy target coordinates to CPU if they are on GPU
    coords = CC.to_cpu(CC.Fields.coordinate_field(target_space))
    lats = CC.Fields.field2array(coords.lat)
    lons = CC.Fields.field2array(coords.long)
    hcoords = CC.Geometry.LatLongPoint.(lats, lons)

    # Remap the field, using MPI if applicable, NVTX-labeled
    range_name =
        isnothing(label) ? "Interfacer.remap:interpolate" :
        "Interfacer.remap:interpolate (" * String(label) * ")"
    NVTX.range_push(; message = range_name)
    out_field = begin
        if comms_ctx isa ClimaComms.SingletonCommsContext
            # Remap source field to target space as an array
            remapped_array = CC.Remapping.interpolate(field, hcoords, [])
            # Convert remapped array to a field in the target space
            CC.Fields.array2field(remapped_array, target_space)
        else
            # Gather then broadcast the global hcoords and offsets
            offset = [length(hcoords)]
            all_hcoords =
                ClimaComms.bcast(comms_ctx, ClimaComms.gather(comms_ctx, hcoords))
            all_offsets =
                ClimaComms.bcast(comms_ctx, ClimaComms.gather(comms_ctx, offset))
            # Interpolate on root and broadcast to all processes
            remapper = CC.Remapping.Remapper(source_space; target_hcoords = all_hcoords)
            remapped_array =
                ClimaComms.bcast(comms_ctx, CC.Remapping.interpolate(remapper, field))
            my_ending_offset = sum(all_offsets[1:ClimaComms.mypid(comms_ctx)])
            my_starting_offset = my_ending_offset - offset[]
            # Convert remapped array to a field in the target space on each process
            CC.Fields.array2field(
                remapped_array[(1 + my_starting_offset):my_ending_offset],
                target_space,
            )
        end
    end
    NVTX.range_pop()
    return out_field
end

function remap(
    num::Number,
    target_space::CC.Spaces.AbstractSpace;
    label::Union{Nothing, String, Symbol} = nothing,
)
    return num
end

"""
    remap!(target_field, source)

Remap the given `source` onto the `target_field`.

Non-ClimaCore fields should provide a method to [`Interfacer.remap`](@ref), or directly to this
function.
"""
function remap!(target_field, source; label::Union{Nothing, String, Symbol} = nothing)
    range_name =
        isnothing(label) ? "Interfacer.remap!" : "Interfacer.remap! (" * String(label) * ")"
    NVTX.@range range_name begin
        target_field .= remap(source, axes(target_field); label)
    end
    return nothing
end


"""
    set_cache!(sim::ComponentModelSimulation)

Perform any initialization of the component model cache that must be done
after the initial exchange.
This is not required to be extended, but may be necessary for some models.
"""
set_cache!(sim::ComponentModelSimulation) = nothing

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


end # module
