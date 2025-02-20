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

export get_model_prog_state, checkpoint_model_state, checkpoint_sims, auto_restart!

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

Checkpoints the model state of a simulation to a HDF5 file at a given time, t (in seconds).
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

Checkpoints the model cache to a N JLD2 files at a given time, t (in seconds),
where N is the number of MPI ranks.
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
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            Checkpointer.checkpoint_model_state(
                sim,
                cs.comms_ctx,
                Int((t - t0) / 1e3),
                output_dir = cs.dirs.checkpoints,
            )
            Checkpointer.checkpoint_model_cache(
                sim,
                cs.comms_ctx,
                Int((t - t0) / 1e3),
                output_dir = cs.dirs.checkpoints,
            )
        end
    end

    # Checkpoint the Coupler fields
    output_dir = cs.dirs.checkpoints
    comms_ctx = cs.comms_ctx
    time = Int((t - t0) / 1e3)
    day = floor(Int, time / (60 * 60 * 24))
    sec = floor(Int, time % (60 * 60 * 24))
    pid = ClimaComms.mypid(comms_ctx)
    @info "Saving coupler fields to JLD2 on day $day second $sec"
    output_file = joinpath(output_dir, "checkpoint_coupler_fields_$(pid)_$time.jld2")
    JLD2.jldsave(output_file, coupler_fields = CC.Adapt.adapt(Array, cs.fields))
end

"""
    auto_restart!(cs::CoupledSimulation, checkpoint_dir)::Bool

Try finding restart files for `cs` in `checkpoint_dir` and use them to reset the
content of `cs`.

Return a true if the simulation was restarted.
"""
function auto_restart!(cs, checkpoint_dir)
    comms_ctx = cs.comms_ctx
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            auto_restart_model_state!(sim, comms_ctx, checkpoint_dir)
            auto_restart_model_cache!(sim, comms_ctx, checkpoint_dir)
        end
    end
    auto_restart_coupler_fields!(cs, checkpoint_dir)
end

"""
    restart!(cs::CoupledSimulation, checkpoint_dir, checkpoint_t)

Overwrite the content of `cs` with checkpoints in `checkpoint_dir` at time `checkpoint_t`.

Return a true if the simulation was restarted.
"""
function restart!(cs, checkpoint_dir, checkpoint_t)
    pid = ClimaComms.mypid(cs.comms_ctx)
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            input_file_state =
                output_file = joinpath(checkpoint_dir, "checkpoint_$(Interfacer.name(sim))_$(checkpoint_t).hdf5")
            restart_model_state!(sim, input_file_state, cs.comms_ctx)
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
    auto_restart_coupler_fields!(cs, checkpoint_dir)

Find checkpoints in `checkpoint_dir` for the coupler fields and write them in `cs`.
"""
function auto_restart_coupler_fields!(cs, checkpoint_dir)
    isdir(checkpoint_dir) || return false
    pid = ClimaComms.mypid(cs.comms_ctx)
    restart_file_rx = Regex("checkpoint_coupler_fields_$(pid)_\\d+.jld2")
    restarts = filter(f -> !isnothing(match(restart_file_rx, f)), readdir(checkpoint_dir))
    isempty(restarts) && return false
    input_file = joinpath(checkpoint_dir, last(sort_by_creation_time(restarts)))
    @info "Restoring coupler fields from $(input_file)"
    restart_coupler_fields!(cs, input_file)
    return true
end

"""
    restart_model_cache!(sim::Interfacer.ComponentModelSimulation, comms_ctx::ClimaComms.AbstractCommsContext, t::Int; input_dir = "input")

Sets the model cache of a simulation from a collection of JLD2 files from a given time, t (in seconds).
"""
function auto_restart_model_cache!(
    sim::Interfacer.ComponentModelSimulation,
    comms_ctx::ClimaComms.AbstractCommsContext,
    checkpoint_dir,
)
    isdir(checkpoint_dir) || return false
    pid = ClimaComms.mypid(comms_ctx)
    restart_file_rx = Regex("checkpoint_cache_$(pid)_$(Interfacer.name(sim))_\\d+.jld2")
    restarts = filter(f -> !isnothing(match(restart_file_rx, f)), readdir(checkpoint_dir))
    isempty(restarts) && return false
    input_file = joinpath(checkpoint_dir, last(sort_by_creation_time(restarts)))
    @info "Restarting $(Interfacer.name(sim)) cache from $input_file"
    return true
end

function restart_model_cache!(sim, input_file)
    ispath(input_file) || error("File $(input_file) not found")
    restore_cache!(sim, JLD2.jldopen(input_file)["cache"])
end

function restart_model_state!(sim, input_file, comms_ctx)
    ispath(input_file) || error("File $(input_file) not found")
    Y = get_model_prog_state(sim)
    # open file and read
    restart_reader = CC.InputOutput.HDF5Reader(input_file, comms_ctx)
    Y_new = CC.InputOutput.read_field(restart_reader, "model_state")
    close(restart_reader)

    # set new state
    Y .= Y_new
end

function restart_coupler_fields!(cs, input_file)
    ispath(input_file) || error("File $(input_file) not found")
    fields_read = JLD2.jldopen(input_file)["coupler_fields"]
    for name in propertynames(cs.fields)
        ArrayType = ClimaComms.array_type(ClimaComms.device(cs.comms_ctx))
        parent(getproperty(cs.fields, name)) .= ArrayType(parent(getproperty(fields_read, name)))
    end
end

"""
    auto_restart_model_state!(sim::Interfacer.ComponentModelSimulation,
                                    comms_ctx::ClimaComms.AbstractCommsContext;
                                    input_dir = "input")

Look for restart files for `sim` in `input_dir`. If anything exists, reset the
state in `sim` to what is in the file.

Detecting restart files is done with
`ClimaUtilities.OutputPathGenerator.detect_restart_file`.
"""
function auto_restart_model_state!(
    sim::Interfacer.ComponentModelSimulation,
    comms_ctx::ClimaComms.AbstractCommsContext,
    checkpoint_dir,
)
    isdir(checkpoint_dir) || return false
    restart_file_rx = Regex("checkpoint_$(Interfacer.name(sim))_\\d+.hdf5")
    restarts = filter(f -> !isnothing(match(restart_file_rx, f)), readdir(checkpoint_dir))
    isempty(restarts) && return false
    input_file = joinpath(checkpoint_dir, last(sort_by_creation_time(restarts)))

    # TODO: Add check that all t_restart for all the component_models are the same
    @info "Restarting $(Interfacer.name(sim)) state from $input_file"
    restart_model_state!(sim, input_file, comms_ctx)
    return true
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
