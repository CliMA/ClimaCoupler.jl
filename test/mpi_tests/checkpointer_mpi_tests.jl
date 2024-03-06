#=
    Unit tests for ClimaCoupler Checkpointer module functions to exercise MPI

These are in a separate testing file from the other Checkpointer unit tests so
that MPI can be enabled for testing of these functions.
=#

using ClimaCore: Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaCoupler: TestHelper
using ClimaComms
using Test
import ClimaCoupler: Interfacer
import ClimaCoupler.Checkpointer: get_model_prog_state, restart_model_state!, checkpoint_model_state

# set up MPI communications context
const comms_ctx = ClimaComms.context(ClimaComms.CPUSingleThreaded())
const pid, nprocs = ClimaComms.init(comms_ctx)
@info pid
ClimaComms.barrier(comms_ctx)

FT = Float64
struct DummySimulation{S} <: Interfacer.AtmosModelSimulation
    state::S
end
get_model_prog_state(sim::DummySimulation) = sim.state
@testset "checkpoint_model_state, restart_model_state!" begin
    boundary_space = TestHelper.create_space(FT, comms_ctx = comms_ctx)
    t = 1

    # old sim run
    sim = DummySimulation(Fields.FieldVector(T = ones(boundary_space)))
    checkpoint_model_state(sim, comms_ctx, t, output_dir = "test_checkpoint")
    ClimaComms.barrier(comms_ctx)

    # new sim run
    sim_new = DummySimulation(Fields.FieldVector(T = zeros(boundary_space)))
    restart_model_state!(sim_new, comms_ctx, t, input_dir = "test_checkpoint")
    @test sim_new.state.T == sim.state.T

    # remove checkpoint directory
    ClimaComms.barrier(comms_ctx)
    if ClimaComms.iamroot(comms_ctx)
        rm("./test_checkpoint/", force = true, recursive = true)
    end
end
