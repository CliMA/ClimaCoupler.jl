"""
    Checkpointer

This module contains template functions for checkpointing the model states and restarting the simulations from the Coupler.
"""
module Checkpointer

import ClimaComms
import ClimaCore as CC
import CC.DataLayouts: AbstractData
import ClimaUtilities.Utils: sort_by_creation_time
import ClimaUtilities.TimeManager: ITime, seconds
import ClimaUtilities.TimeVaryingInputs: AbstractTimeVaryingInput
import Dates
import JLD2
import StaticArrays
import ..Interfacer

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
    p = CC.Adapt.adapt(Array, get_model_cache(sim))
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
    restore!(v1, v2, comms_ctx; ignore)

Recursively traverse `v1` and `v2`, setting each field of `v1` with the
corresponding field in `v2`. In this, ignore all the properties that have name
within the `ignore` iterable.

`ignore` is useful when there are stateful properties, such as live pointers.
"""
function restore!(v1::T1, v2::T2, comms_ctx; name = "", ignore) where {T1, T2}
    # We pick fieldnames(T2) because v2 tend to be simpler (Array as opposed
    # to CuArray)
    fields = filter(x -> !(x in ignore), fieldnames(T2))
    if isempty(fields)
        v1 == v2 || error("$v1 != $v2")
    else
        # Recursive case
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

# Ignoring certain types that don't need to be restored
# UnionAll and DataType are infinitely recursive, so we also ignore those
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
    name,
    ignore,
)
    return nothing
end

function restore!(
    v1::Union{AbstractData, AbstractArray},
    v2::Union{AbstractData, AbstractArray},
    comms_ctx;
    name,
    ignore,
)
    ArrayType =
        parent(v1) isa Array ? Array : ClimaComms.array_type(ClimaComms.device(comms_ctx))
    moved_to_device = ArrayType(parent(v2))

    parent(v1) .= moved_to_device
    return nothing
end

function restore!(
    v1::Union{StaticArrays.StaticArray, Number, UnitRange, LinRange, Symbol},
    v2::Union{StaticArrays.StaticArray, Number, UnitRange, LinRange, Symbol},
    comms_ctx;
    name,
    ignore,
)
    v1 == v2 || error("$name is a immutable but it inconsistent ($(v1) != $(v2))")
    return nothing
end

function restore!(v1::Dict, v2::Dict, comms_ctx; name, ignore)
    # RRTGMP has some internal dictionaries
    v1 == v2 || error("$name is inconsistent")
    return nothing
end

"""
    restore!(v1::T1, v2::T2, comms_ctx; name = "", ignore) where {T1 <: Union{Dates.DateTime, Dates.UTInstant, Dates.Millisecond}, T2 <: Union{Dates.DateTime, Dates.UTInstant, Dates.Millisecond}}

Special case for time-related types to allow different timestamps during restore.
"""
function restore!(
    v1::T1,
    v2::T2,
    comms_ctx;
    name,
    ignore,
) where {
    T1 <: Union{Dates.DateTime, Dates.UTInstant, Dates.Millisecond},
    T2 <: Union{Dates.DateTime, Dates.UTInstant, Dates.Millisecond},
}
    if v1 != v2
        @warn "Time value differs in restart" field = name original = v2 new = v1
    end
    return nothing
end

end # module
