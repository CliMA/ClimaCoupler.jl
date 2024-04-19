"""
    Checkpointer

This module contains template functions for checkpointing the model states and restarting the simulations from the Coupler.
"""
module Checkpointer

import ClimaComms
import ClimaCore as CC
import ..Interfacer

export get_model_prog_state, checkpoint_model_state, restart_model_state!

"""
    get_model_prog_state(sim::Interfacer.ComponentModelSimulation)

Returns the model state of a simulation as a `ClimaCore.FieldVector`.
This is a template function that should be implemented for each component model.
"""
get_model_prog_state(sim::Interfacer.ComponentModelSimulation) = nothing

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

    @info "Setting " Interfacer.name(sim) " state to checkpoint: $input_file, corresponding to day $day second $sec"

    # open file and read
    restart_reader = CC.InputOutput.HDF5Reader(input_file, comms_ctx)
    Y_new = CC.InputOutput.read_field(restart_reader, "model_state")
    Base.close(restart_reader)

    # set new state
    Y .= Y_new
end

end # module
