"""
    Checkpointer

This module contains template functions for checkpointing the model states and restarting the simulations from the Coupler.
"""
module Checkpointer

import ClimaComms
import ClimaCore as CC
import ClimaUtilities.Utils: sort_by_creation_time
import ClimaUtilities.TimeManager: ITime, seconds
import ClimaUtilities.TimeVaryingInputs: AbstractTimeVaryingInput
import ..Interfacer
import Dates
import StaticArrays

import JLD2

export get_model_prog_state, checkpoint_model_state, checkpoint_sims, restore!

"""
    get_model_prog_state(sim::Interfacer.ComponentModelSimulation)

Returns the model state of a simulation as a `ClimaCore.FieldVector`.
This is a template function that should be implemented for each component model.
"""
get_model_prog_state(sim::Interfacer.ComponentModelSimulation) = nothing

"""
    get_model_cache(sim::Interfacer.ComponentModelSimulation)

Returns the model cache of a simulation.
This is a template function that should be implemented for each component model.
"""
get_model_cache(sim::Interfacer.ComponentModelSimulation) = nothing

"""
    checkpoint_model_state(
        sim::Interfacer.ComponentModelSimulation,
        comms_ctx::ClimaComms.AbstractCommsContext,
        t::Int,
        prev_checkpoint_t::Int;
        output_dir = "output")

Checkpoint the model state of a simulation to a HDF5 file at a given time, t (in seconds).

If a previous checkpoint exists, it is removed. This is to avoid accumulating
many checkpoint files in the output directory. A value of -1 for `prev_checkpoint_t`
is used to indicate that there is no previous checkpoint to remove.
"""
function checkpoint_model_state(
    sim::Interfacer.ComponentModelSimulation,
    comms_ctx::ClimaComms.AbstractCommsContext,
    t::Int,
    prev_checkpoint_t::Int;
    output_dir = "output",
)
    Y = get_model_prog_state(sim)
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    @info "Saving checkpoint $(nameof(sim)) model state to HDF5 on day $day second $sec"
    output_file = joinpath(output_dir, "checkpoint_$(nameof(sim))_$t.hdf5")
    checkpoint_writer = CC.InputOutput.HDF5Writer(output_file, comms_ctx)
    CC.InputOutput.HDF5.write_attribute(checkpoint_writer.file, "time", t)
    CC.InputOutput.write!(checkpoint_writer, Y, "model_state")
    Base.close(checkpoint_writer)

    # Remove previous checkpoint if it exists
    prev_checkpoint_file =
        joinpath(output_dir, "checkpoint_$(nameof(sim))_$(prev_checkpoint_t).hdf5")
    remove_checkpoint(prev_checkpoint_file, prev_checkpoint_t, comms_ctx)
    return nothing
end

"""
    checkpoint_model_cache(
        sim::Interfacer.ComponentModelSimulation,
        comms_ctx::ClimaComms.AbstractCommsContext,
        t::Int,
        prev_checkpoint_t::Int;
        output_dir = "output")

Checkpoint the model cache to N JLD2 files at a given time, t (in seconds),
where N is the number of MPI ranks.

Objects are saved to JLD2 files because caches are generally not ClimaCore
objects (and ClimaCore.InputOutput can only save `Field`s or `FieldVector`s).

If a previous checkpoint exists, it is removed. This is to avoid accumulating
many checkpoint files in the output directory. A value of -1 for `prev_checkpoint_t`
is used to indicate that there is no previous checkpoint to remove.
"""
function checkpoint_model_cache(
    sim::Interfacer.ComponentModelSimulation,
    comms_ctx::ClimaComms.AbstractCommsContext,
    t::Int,
    prev_checkpoint_t::Int;
    output_dir = "output",
)
    # Move p to CPU (because we cannot save CUArrays)
    p = get_model_cache_to_checkpoint(sim)
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    @info "Saving checkpoint $(nameof(sim)) model cache to JLD2 on day $day second $sec"
    pid = ClimaComms.mypid(comms_ctx)
    output_file = joinpath(output_dir, "checkpoint_cache_$(pid)_$(nameof(sim))_$t.jld2")
    JLD2.jldsave(output_file, cache = p)

    # Remove previous checkpoint if it exists
    prev_checkpoint_file = joinpath(
        output_dir,
        "checkpoint_cache_$(pid)_$(nameof(sim))_$(prev_checkpoint_t).jld2",
    )
    remove_checkpoint(prev_checkpoint_file, prev_checkpoint_t, comms_ctx)
    return nothing
end

"""
    get_model_cache_to_checkpoint(sim::Interfacer.ComponentModelSimulation)

Prepare the cache for checkpointing by moving the entire cache to CPU.
"""
function get_model_cache_to_checkpoint(sim::Interfacer.ComponentModelSimulation)
    return CC.Adapt.adapt(Array, get_model_cache(sim))
end

"""
    get_model_cache_to_checkpoint(sim::Interfacer.AtmosModelSimulation)

Prepare the atmos cache for checkpoint by selectively moving parts of the atmos
cache to CPU instead of moving the entire atmos cache to CPU, resulting in a
much smaller saved file.

# Implementation Details

When moving the cache from GPU to CPU, calling `adapt` on the entire cache
creates unnecessary duplicate objects because `adapt` is not properly defined
for the entire cache structure. This function addresses three key issues:

1. **Individual adaptation**: On GPU, `adapt` is called on each object
   separately rather than on the entire cache at once.

2. **Deduplication**: Objects sharing the same object ID are not duplicated.
   Instead, references to already-processed objects are reused.

3. **Selective saving**: Only the parts of the cache needed for restoration are
   saved.

4. **Recreated views**: Views are preserved by recreating the views with the
   parent arrays on CPU. Otherwise, the arrays are unnecessarily duplicated when
   `adapt` is called since checks using `objectid` fail for views, resulting in
   a larger saved cache on GPU.

# Returns
A vector of objects from the cache. Elements may reference the same underlying
data if they share object IDs. The order of the objects in the vector is
determined by `CacheIterator`.
"""
function get_model_cache_to_checkpoint(sim::Interfacer.AtmosModelSimulation)
    atmos_cache_itr = CacheIterator(sim)
    cache_vec = [] # Elements of the vector can be any type
    # Keep track of the object ID and its index in the vector while iterating
    # through the fields of the cache
    obj_id_to_idx_map = Dict{UInt64, Int}()
    for (i, obj_saved) in enumerate(atmos_cache_itr)
        obj_id = objectid(obj_saved)
        if obj_id in keys(obj_id_to_idx_map)
            # Reuse the same object if we encountered the same object ID while
            # traversing the object
            push!(cache_vec, cache_vec[obj_id_to_idx_map[obj_id]])
        else
            # Call adapt on each individual object
            push!(cache_vec, CC.Adapt.adapt(Array, obj_saved))
            obj_id_to_idx_map[obj_id] = i
        end
    end

    # Recreate views if needed
    atmos_cache_itr = CacheIterator(sim)
    for (i, (obj_cache, obj_saved)) in enumerate(zip(cache_vec, atmos_cache_itr))
        obj_cache isa SubArray || continue
        # When the SubArray is adapted to CPU, its parent array is lost. Here we
        # look up the parent array's object ID in the cache, find the
        # corresponding parent array that is on CPU, and reconstruct the view
        # with the parent array.
        parent_obj_saved = parent(obj_saved)
        parent_obj_saved_obj_id = objectid(parent_obj_saved)
        if parent_obj_saved_obj_id in keys(obj_id_to_idx_map)
            parent_obj_cache = cache_vec[obj_id_to_idx_map[parent_obj_saved_obj_id]]
            recreated_obj_saved = view(parent_obj_cache, parentindices(obj_cache)...)
            if !isequal(obj_cache, recreated_obj_saved)
                @warn("Object at index $i in the cache vector is not identical to the view")
                continue
            else
                cache_vec[i] = recreated_obj_saved
            end
        end
    end

    return cache_vec
end

"""
    restore_cache!(sim::Interfacer.ComponentModelSimulation, new_cache)

Replace the cache in `sim` with `new_cache`.

Component models can define new methods for this to change how cache is restored.
"""
function restore_cache!(sim::Interfacer.ComponentModelSimulation, new_cache)
    return nothing
end

"""
    checkpoint_sims(cs::CoupledSimulation)

This is a callback function that checkpoints all simulations defined in the
current coupled simulation.
"""
function checkpoint_sims(cs::Interfacer.CoupledSimulation)
    time = Int(round(float(cs.t[])))
    day = floor(Int, time / (60 * 60 * 24))
    sec = floor(Int, time % (60 * 60 * 24))
    output_dir = cs.dir_paths.checkpoints_dir
    prev_checkpoint_t = cs.prev_checkpoint_t[]
    comms_ctx = ClimaComms.context(cs)

    for sim in cs.model_sims
        if !isnothing(Checkpointer.get_model_prog_state(sim))
            Checkpointer.checkpoint_model_state(
                sim,
                comms_ctx,
                time,
                prev_checkpoint_t;
                output_dir,
            )
        end
        if !isnothing(Checkpointer.get_model_cache(sim)) && cs.save_cache
            Checkpointer.checkpoint_model_cache(
                sim,
                comms_ctx,
                time,
                prev_checkpoint_t;
                output_dir,
            )
        end
    end

    # Checkpoint the Coupler fields
    pid = ClimaComms.mypid(comms_ctx)
    @info "Saving coupler fields to JLD2 on day $day second $sec"
    output_file = joinpath(output_dir, "checkpoint_coupler_fields_$(pid)_$time.jld2")
    # Adapt to Array move fields to the CPU
    JLD2.jldsave(output_file, coupler_fields = CC.Adapt.adapt(Array, cs.fields))

    # Remove previous Coupler fields checkpoint if it exists
    prev_checkpoint_file =
        joinpath(output_dir, "checkpoint_coupler_fields_$(pid)_$(prev_checkpoint_t).jld2")
    remove_checkpoint(prev_checkpoint_file, prev_checkpoint_t, comms_ctx)

    # Update previous checkpoint time stored in the coupled simulation
    cs.prev_checkpoint_t[] = time
end

"""
    restart!(cs::CoupledSimulation, checkpoint_dir, checkpoint_t, restart_cache)

Overwrite the content of `cs` with checkpoints in `checkpoint_dir` at time `checkpoint_t`.

If `restart_cache` is true, the cache will be read from the restart file using `restore_cache!`.
Otherwise, the cache will be left unchanged.

Return a true if the simulation was restarted.
"""
function restart!(cs, checkpoint_dir, checkpoint_t, restart_cache)
    @info "Restarting from time $(checkpoint_t) and directory $(checkpoint_dir)"
    pid = ClimaComms.mypid(ClimaComms.context(cs))
    for sim in cs.model_sims
        if !isnothing(Checkpointer.get_model_prog_state(sim))
            input_file_state =
                output_file = joinpath(
                    checkpoint_dir,
                    "checkpoint_$(nameof(sim))_$(checkpoint_t).hdf5",
                )
            restart_model_state!(sim, input_file_state, ClimaComms.context(cs))
        end
        if !isnothing(Checkpointer.get_model_cache(sim)) && restart_cache
            input_file_cache = joinpath(
                checkpoint_dir,
                "checkpoint_cache_$(pid)_$(nameof(sim))_$(checkpoint_t).jld2",
            )
            restart_model_cache!(sim, input_file_cache)
        end
    end
    input_file_coupler_fields =
        joinpath(checkpoint_dir, "checkpoint_coupler_fields_$(pid)_$(checkpoint_t).jld2")
    restart_coupler_fields!(cs, input_file_coupler_fields)
    return true
end

"""
    restart_model_cache!(sim, input_file)

Overwrite the content of `sim` with the cache from the `input_file`.

It relies on `restore_cache!(sim, old_cache)`, which has to be implemented by
the component models that have a cache.
"""
function restart_model_cache!(sim, input_file)
    ispath(input_file) || error("File $(input_file) not found")
    # Component models are responsible for defining a method for this
    JLD2.jldopen(input_file) do file
        restore_cache!(sim, file["cache"])
    end
end

"""
    restart_model_state!(sim, input_file, comms_ctx)

Overwrite the content of `sim` with the state from the `input_file`.
"""
function restart_model_state!(sim, input_file, comms_ctx)
    ispath(input_file) || error("File $(input_file) not found")
    Y = get_model_prog_state(sim)
    # open file and read
    CC.InputOutput.HDF5Reader(input_file, comms_ctx) do restart_reader
        Y_new = CC.InputOutput.read_field(restart_reader, "model_state")
        # set new state
        Y .= Y_new
    end
    return nothing
end

"""
    restart_coupler_fields!(cs, input_file)

Overwrite the content of the coupled simulation `cs` with the coupler fields
read from `input_file`.
"""
function restart_coupler_fields!(cs, input_file)
    ispath(input_file) || error("File $(input_file) not found")
    JLD2.jldopen(input_file) do file
        fields_read = file["coupler_fields"]
        for name in propertynames(cs.fields)
            ArrayType = ClimaComms.array_type(ClimaComms.device(cs))
            parent(getproperty(cs.fields, name)) .=
                ArrayType(parent(getproperty(fields_read, name)))
        end
    end
end

"""
    t_start_from_checkpoint(checkpoint_dir)

Look for restart files in `checkpoint_dir`, if found, return the time of the latest.
If not found, return `nothing`.
"""
function t_start_from_checkpoint(checkpoint_dir)
    isdir(checkpoint_dir) || return nothing
    restart_file_rx = r"checkpoint_(\w+)_(\d+).hdf5"
    restarts = filter(f -> !isnothing(match(restart_file_rx, f)), readdir(checkpoint_dir))
    isempty(restarts) && return nothing
    latest_restart = last(sort_by_creation_time(restarts))
    return parse(Int, match(restart_file_rx, latest_restart)[2])
end

"""
    remove_checkpoint(prev_checkpoint_file, prev_checkpoint_t, comms_ctx)

Delete the provided checkpoint file on the root process and print a helpful
info message. This can be used to remove intermediate checkpoints, to prevent
saving excessively large amounts of output.
"""
function remove_checkpoint(prev_checkpoint_file, prev_checkpoint_t, comms_ctx)
    if ClimaComms.iamroot(comms_ctx) &&
       prev_checkpoint_t != -1 &&
       isfile(prev_checkpoint_file)
        @info "Removing previous checkpoint file: $prev_checkpoint_file"
        rm(prev_checkpoint_file)
    end
    return nothing
end

"""
    restore!(v1, v2, comms_ctx; name = "", ignore = Set())

Recursively traverse `v1` and `v2`, setting each field of `v1` with the
corresponding field in `v2`. In this, ignore all the properties that have name
within the `ignore` iterable.

This is intended to be used when restarting a simulation's cache object
from a checkpoint.

`ignore` is useful when there are stateful properties, such as live pointers.
"""
function restore!(v1::T1, v2::T2, comms_ctx; name = "", ignore = Set()) where {T1, T2}
    # We pick fieldnames(T2) because v2 tend to be simpler (Array as opposed
    # to CuArray)
    fields = filter(x -> !(x in ignore), fieldnames(T2))
    # If there are no fields to restore, we check for consistency
    if isempty(fields)
        v1 == v2 || error("$v1 != $v2")
    else
        # Recursive case: restore each field
        for p in fields
            restore!(
                getfield(v1, p),
                getfield(v2, p),
                comms_ctx;
                name = "$(name).$(p)",
                ignore,
            )
        end
    end
    return nothing
end

"""
    restore!(
        v1::Union{
            AbstractTimeVaryingInput,
            ClimaComms.AbstractCommsContext,
            ClimaComms.AbstractDevice,
            UnionAll,
            DataType,
        },
        v2::Union{
            AbstractTimeVaryingInput,
            ClimaComms.AbstractCommsContext,
            ClimaComms.AbstractDevice,
            UnionAll,
            DataType,
        },
        _comms_ctx;
        name = "",
        ignore = Set(),
    )

Ignore certain types that don't need to be restored.
`UnionAll` and `DataType` are infinitely recursive, so we also ignore those.
"""
function restore!(
    v1::Union{
        AbstractTimeVaryingInput,
        ClimaComms.AbstractCommsContext,
        ClimaComms.AbstractDevice,
        UnionAll,
        DataType,
    },
    v2::Union{
        AbstractTimeVaryingInput,
        ClimaComms.AbstractCommsContext,
        ClimaComms.AbstractDevice,
        UnionAll,
        DataType,
    },
    _comms_ctx;
    name = "",
    ignore = Set(),
)
    return nothing
end

"""
    restore!(
        v1::Union{CC.DataLayouts.AbstractData, AbstractArray},
        v2::Union{CC.DataLayouts.AbstractData, AbstractArray},
        comms_ctx;
        name = "",
        ignore = Set(),
    )

For array-like objects, we move the original data (v2) to the
device of the new data (v1). Then we copy the original data to
the new object.
"""
function restore!(
    v1::Union{CC.DataLayouts.AbstractData, AbstractArray},
    v2::Union{CC.DataLayouts.AbstractData, AbstractArray},
    comms_ctx;
    name = "",
    ignore = Set(),
)
    ArrayType =
        parent(v1) isa Array ? Array : ClimaComms.array_type(ClimaComms.device(comms_ctx))
    moved_to_device = ArrayType(parent(v2))

    parent(v1) .= moved_to_device
    return nothing
end

"""
    restore!(v1::LinearIndices, v2::AbstractArray, comms_ctx; name = "", ignore = Set())

Special case to compare LinearIndices to AbstractArray, which is needed for ClimaAtmos v0.32.
"""
function restore!(
    v1::LinearIndices,
    v2::AbstractArray,
    comms_ctx;
    name = "",
    ignore = Set(),
)
    v1 == v2 || error("$name is immutable but it inconsistent ($(v1) != $(v2))")
    return nothing
end

"""
    restore!(
        v1::Union{StaticArrays.StaticArray, Number, UnitRange, LinRange, Symbol},
        v2::Union{StaticArrays.StaticArray, Number, UnitRange, LinRange, Symbol},
        comms_ctx;
        name = "",
        ignore = Set(),
    )

Ensure that immutable objects have been initialized correctly,
as they cannot be restored from a checkpoint.
"""
function restore!(
    v1::Union{
        StaticArrays.StaticArray,
        Number,
        UnitRange,
        LinRange,
        Symbol,
        CartesianIndices,
    },
    v2::Union{
        StaticArrays.StaticArray,
        Number,
        UnitRange,
        LinRange,
        Symbol,
        CartesianIndices,
    },
    comms_ctx;
    name = "",
    ignore = Set(),
)
    v1 == v2 || error("$name is immutable but it inconsistent ($(v1) != $(v2))")
    return nothing
end

"""
    restore!(v1::Dict, v2::Dict, comms_ctx; name = "", ignore = Set())

RRTMGP has some internal dictionaries, which we check for consistency.
"""
function restore!(v1::Dict, v2::Dict, comms_ctx; name = "", ignore = Set())
    v1 == v2 || error("$name is inconsistent")
    return nothing
end

"""
    restore!(
        v1::T1,
        v2::T2,
        comms_ctx;
        name = "",
        ignore = Set(),
    ) where {
        T1 <: Union{Dates.DateTime, Dates.UTInstant, Dates.Millisecond},
        T2 <: Union{Dates.DateTime, Dates.UTInstant, Dates.Millisecond},
    }

Special case to compare time-related types to allow different timestamps during restore.
"""
function restore!(
    v1::T1,
    v2::T2,
    comms_ctx;
    name = "",
    ignore = Set(),
) where {
    T1 <: Union{Dates.DateTime, Dates.UTInstant, Dates.Millisecond},
    T2 <: Union{Dates.DateTime, Dates.UTInstant, Dates.Millisecond},
}
    if v1 != v2
        @warn "Time value differs in restart" field = name original = v2 new = v1
    end
    return nothing
end

# Needed for CacheIterator, because calling adapt on a StaticArray makes it into
# an array
function restore!(v1::StaticArrays.StaticArray, v2::Array, comms_ctx; name, ignore)
    v1 == v2 || error("$name is a immutable but it inconsistent ($(v1) != $(v2))")
    return nothing
end

"""
    CacheIterator{F, G}

An iterator that recursively iterate through the object's field hierarchy.
"""
struct CacheIterator
    """A vector to hold the objects during depth first search of the object's
    field hierarchy"""
    stack::Vector

    """A set of field names (as `Symbols`s) to exclude from traversal."""
    ignore::Set{Symbol}
end

"""
    CacheIterator(cache, ignore::Set{Symbol})

A stateful iterator for recursively traversing the fields of `cache` by reverse
preorder traversal.

# Arguments
- `cache `: The cache object to traverse.
- `ignore`: A set of field names (as `Symbol`s) to skip during traversal.

# Returns
An iterator that yields objects that should be saved.

# Notes
- Objects returned by the iterator may reference the same memory.
- Fields listed in `ignore` are not traversed.
"""
function CacheIterator(cache, ignore::Set{Symbol})
    stack = []
    push!(stack, cache)
    return CacheIterator(stack, ignore)
end

"""
    iterate(itr::CacheIterator)

Advance the iterator to obtain the next field in the field hierarchy of the
object. If no fields remain, `nothing` is returned. Otherwise, a 2-tuple of the
next element and the new iteration state is returned.
"""
function Base.iterate(itr::CacheIterator)
    stack = itr.stack
    ignore = itr.ignore
    while !isempty(stack)
        obj = pop!(stack)
        stop(obj) && continue
        is_leaf(obj) && return (obj, nothing)
        T = typeof(obj)
        fields = filter(x -> !(x in ignore), fieldnames(T))
        for p in fields
            push!(stack, getfield(obj, p))
        end
    end
    return nothing
end

# Functions needed to implement the iterator interface
# see https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration
Base.iterate(iter::CacheIterator, _) = iterate(iter)
Base.IteratorEltype(::CacheIterator) = Base.EltypeUnknown()
Base.IteratorSize(::CacheIterator) = Base.SizeUnknown()
Base.isdone(itr::CacheIterator, state...) = isempty(itr.stack)

"""
    stop(::Type)

Returns `true` to stop recursing at the given object and `false` to continue
recursing through the fields of the object.

This function is used by `CacheIterator`. This is useful to stop recursing and
not return the object.
"""
function stop(
    ::Union{
        AbstractTimeVaryingInput,
        ClimaComms.AbstractCommsContext,
        ClimaComms.AbstractDevice,
        UnionAll,
        DataType,
    },
)
    return true # Stop recursing here
end

# Catch all for everything else
stop(obj) = false

"""
    is_leaf(::Type)

Returns `true` if the object should be saved and `false` otherwise.

This function is used by `CacheIterator`.
"""
is_leaf(::Union{CC.DataLayouts.AbstractData, AbstractArray}) = true # Needed for saving data

# Needed for error handling
is_leaf(
    ::Union{
        StaticArrays.StaticArray,
        Number,
        UnitRange,
        LinRange,
        Symbol,
        Dict,
        Dates.DateTime,
        Dates.UTInstant,
        Dates.Millisecond,
    },
) = true

# Catch all for everything else
is_leaf(::Any) = false

"""
    CacheIterator(sim::Interfacer.ComponentModelSimulation)

Create an iterator over the cache in `sim` of objects that should be
saved.

To use this constructor, the function `get_cache_ignore` must be implemented
for `sim` which takes in `sim` and returns a set of symbols to ignore while
iterating through the fields of the cache.
"""
function CacheIterator(sim::Interfacer.ComponentModelSimulation)
    cache = get_model_cache(sim)
    return CacheIterator(cache, get_cache_ignore(sim))
end

"""
    get_cache_ignore(::Interfacer.AtmosModelSimulation)

Return a set of symbols that should be ignored when iterating the fields of the
atmos cache. These fields will not be saved when checkpointing.
"""
function get_cache_ignore(::Interfacer.AtmosModelSimulation)
    return Set([
        :rc,
        :params,
        :ghost_buffer,
        :hyperdiffusion_ghost_buffer,
        :data_handler,
        :graph_context,
        :dt,
    ])
end

function get_cache_ignore(::Interfacer.ComponentModelSimulation)
    return error("List of symbols to ignore when iterating the fields of the
    cache is not supported for this model. Checkpoint by saving the entire
    cache instead or define get_cache_ignore to use CacheIterator.")
end

end # module
