"""
    Checkpointer

This module contains template functions for checkpointing the model states and restarting the simulations from the Coupler.
"""
module Checkpointer

import ClimaComms
import ClimaCore as CC
import ClimaUtilities.Utils: sort_by_creation_time
import ..Interfacer
import Dates

import JLD2

export get_model_prog_state, checkpoint_model_state, checkpoint_sims

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
    checkpoint_model_state(sim::Interfacer.ComponentModelSimulation, comms_ctx::ClimaComms.AbstractCommsContext, t::Int; output_dir = "output")

Checkpoint the model state of a simulation to a HDF5 file at a given time, t (in seconds).
"""
function checkpoint_model_state(
    sim::Interfacer.ComponentModelSimulation,
    comms_ctx::ClimaComms.AbstractCommsContext,
    t::Int;
    output_dir = "output",
)
    Y = get_model_prog_state(sim)
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    @info "Saving checkpoint " * Interfacer.name(sim) * " model state to HDF5 on day $day second $sec"
    output_file = joinpath(output_dir, "checkpoint_" * Interfacer.name(sim) * "_$t.hdf5")
    checkpoint_writer = CC.InputOutput.HDF5Writer(output_file, comms_ctx)
    CC.InputOutput.HDF5.write_attribute(checkpoint_writer.file, "time", t)
    CC.InputOutput.write!(checkpoint_writer, Y, "model_state")
    Base.close(checkpoint_writer)
    return nothing

end

"""
    checkpoint_model_cache(sim::Interfacer.ComponentModelSimulation, comms_ctx::ClimaComms.AbstractCommsContext, t::Int; output_dir = "output")

Checkpoint the model cache to N JLD2 files at a given time, t (in seconds),
where N is the number of MPI ranks.

Objects are saved to JLD2 files because caches are generally not ClimaCore
objects (and ClimaCore.InputOutput can only save `Field`s or `FieldVector`s).
"""
function checkpoint_model_cache(
    sim::Interfacer.ComponentModelSimulation,
    comms_ctx::ClimaComms.AbstractCommsContext,
    t::Int;
    output_dir = "output",
)
    # Move p to CPU (because we cannot save CUArrays)
    p = CC.Adapt.adapt(Array, get_model_cache(sim))
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    @info "Saving checkpoint " * Interfacer.name(sim) * " model cache to JLD2 on day $day second $sec"
    pid = ClimaComms.mypid(comms_ctx)
    output_file = joinpath(output_dir, "checkpoint_cache_$(pid)_" * Interfacer.name(sim) * "_$t.jld2")
    JLD2.jldsave(output_file, cache = p)
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

This is a callback function that checkpoints all simulations defined in the current coupled simulation.
"""
function checkpoint_sims(cs::Interfacer.CoupledSimulation)
    t = Dates.datetime2epochms(cs.dates.date[1])
    t0 = Dates.datetime2epochms(cs.dates.date0[1])
    time = Int((t - t0) / 1e3)
    day = floor(Int, time / (60 * 60 * 24))
    sec = floor(Int, time % (60 * 60 * 24))
    output_dir = cs.dirs.checkpoints
    comms_ctx = cs.comms_ctx
    for sim in cs.model_sims
        if !isnothing(Checkpointer.get_model_prog_state(sim))
            Checkpointer.checkpoint_model_state(sim, comms_ctx, time; output_dir)
        end
        if !isnothing(Checkpointer.get_model_cache(sim))
            Checkpointer.checkpoint_model_cache(sim, comms_ctx, time; output_dir)
        end
    end
    # Checkpoint the Coupler fields
    pid = ClimaComms.mypid(comms_ctx)
    @info "Saving coupler fields to JLD2 on day $day second $sec"
    output_file = joinpath(output_dir, "checkpoint_coupler_fields_$(pid)_$time.jld2")
    # Adapt to Array move fields to the CPU
    JLD2.jldsave(output_file, coupler_fields = CC.Adapt.adapt(Array, cs.fields))
end

"""
    restart!(cs::CoupledSimulation, checkpoint_dir, checkpoint_t)

Overwrite the content of `cs` with checkpoints in `checkpoint_dir` at time `checkpoint_t`.

Return a true if the simulation was restarted.
"""
function restart!(cs, checkpoint_dir, checkpoint_t)
    @info "Restarting from time $(checkpoint_t) and directory $(checkpoint_dir)"
    pid = ClimaComms.mypid(cs.comms_ctx)
    for sim in cs.model_sims
        if !isnothing(Checkpointer.get_model_prog_state(sim))
            input_file_state =
                output_file = joinpath(checkpoint_dir, "checkpoint_$(Interfacer.name(sim))_$(checkpoint_t).hdf5")
            restart_model_state!(sim, input_file_state, cs.comms_ctx)
        end
        if !isnothing(Checkpointer.get_model_cache(sim))
            input_file_cache =
                joinpath(checkpoint_dir, "checkpoint_cache_$(pid)_$(Interfacer.name(sim))_$(checkpoint_t).jld2")
            restart_model_cache!(sim, input_file_cache)
        end
    end
    input_file_coupler_fields = joinpath(checkpoint_dir, "checkpoint_coupler_fields_$(pid)_$(checkpoint_t).jld2")
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
            ArrayType = ClimaComms.array_type(ClimaComms.device(cs.comms_ctx))
            parent(getproperty(cs.fields, name)) .= ArrayType(parent(getproperty(fields_read, name)))
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

end # module
