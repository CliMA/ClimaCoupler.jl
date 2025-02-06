"""
    Checkpointer

This module contains template functions for checkpointing the model states and restarting the simulations from the Coupler.
"""
module Checkpointer

import ClimaComms
import ClimaCore as CC
import ..Interfacer
import Dates

import JLD2

export get_model_prog_state, checkpoint_model_state, restart_model_state!, checkpoint_sims

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
    mkpath(joinpath(output_dir, "checkpoint"))
    output_file = joinpath(output_dir, "checkpoint", "checkpoint_" * Interfacer.name(sim) * "_$t.hdf5")
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
    p = get_model_cache(sim)
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    @info "Saving checkpoint " * Interfacer.name(sim) * " model cache to JLD2 on day $day second $sec"
    mkpath(joinpath(output_dir, "checkpoint"))
    pid = ClimaComms.mypid(comms_ctx)
    output_file = joinpath(output_dir, "checkpoint", "checkpoint_cache_$(pid)_" * Interfacer.name(sim) * "_$t.jld2")
    JLD2.jldsave(output_file, cache = p)
    return nothing
end


"""
    restart_model_state!(sim::Interfacer.ComponentModelSimulation, comms_ctx::ClimaComms.AbstractCommsContext, t::Int; input_dir = "input")

Sets the model state of a simulation from a HDF5 file from a given time, t (in seconds).
"""
function restart_model_state!(
    sim::Interfacer.ComponentModelSimulation,
    comms_ctx::ClimaComms.AbstractCommsContext,
    t::Int;
    input_dir = "input",
)
    Y = get_model_prog_state(sim)
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    input_file = joinpath(input_dir, "checkpoint", "checkpoint_" * Interfacer.name(sim) * "_$t.hdf5")

    @info "Reading " Interfacer.name(sim) " state from checkpoint: $input_file, corresponding to day $day second $sec"

    # open file and read
    restart_reader = CC.InputOutput.HDF5Reader(input_file, comms_ctx)
    Y_new = CC.InputOutput.read_field(restart_reader, "model_state")
    Base.close(restart_reader)

    # set new state
    Y .= Y_new
end

"""
    restart_model_state!(sim::Interfacer.ComponentModelSimulation, comms_ctx::ClimaComms.AbstractCommsContext, t::Int; input_dir = "input")

Sets the model cache of a simulation from a collection of JLD2 files from a given time, t (in seconds).
"""
function restart_model_cache!(
    sim::Interfacer.ComponentModelSimulation,
    comms_ctx::ClimaComms.AbstractCommsContext,
    t::Int;
    input_dir = "input",
)
    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))
    pid = ClimaComms.mypid(comms_ctx)
    input_file = joinpath(input_dir, "checkpoint", "checkpoint_cache_$(pid)_" * Interfacer.name(sim) * "_$t.jld2")

    @info "Setting " Interfacer.name(sim) " cache to checkpoint: $input_file, corresponding to day $day second $sec"

    restore_cache!(sim, JLD2.jldopen(input_file)["cache"])
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
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            t = Dates.datetime2epochms(cs.dates.date[1])
            t0 = Dates.datetime2epochms(cs.dates.date0[1])
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
end

end # module
