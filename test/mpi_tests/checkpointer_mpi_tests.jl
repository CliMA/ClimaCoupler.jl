#=
    Unit tests for ClimaCoupler Checkpointer module functions to exercise MPI

These are in a separate testing file from the other Checkpointer unit tests so
that MPI can be enabled for testing of these functions.
=#
import Test: @test, @testset
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaCore as CC
import ClimaCoupler
import ClimaCoupler: Checkpointer, Interfacer, TestHelper

# set up MPI communications context
const comms_ctx = ClimaComms.context(ClimaComms.CPUSingleThreaded())
const pid, nprocs = ClimaComms.init(comms_ctx)
@info pid
ClimaComms.barrier(comms_ctx)

FT = Float64
struct DummySimulation{S} <: Interfacer.AtmosModelSimulation
    state::S
end
Checkpointer.get_model_prog_state(sim::DummySimulation) = sim.state
@testset "checkpoint_model_state, restart_model_state!" begin
    boundary_space = TestHelper.create_space(FT, comms_ctx = comms_ctx)
    t = 1

    # old sim run
    sim = DummySimulation(CC.Fields.FieldVector(T = ones(boundary_space)))
    Checkpointer.checkpoint_model_state(sim, comms_ctx, t, output_dir = "test_checkpoint")
    ClimaComms.barrier(comms_ctx)

    # new sim run
    sim_new = DummySimulation(CC.Fields.FieldVector(T = zeros(boundary_space)))
    Checkpointer.restart_model_state!(sim_new, comms_ctx, t, input_dir = "test_checkpoint")
    @test sim_new.state.T == sim.state.T

    # remove checkpoint directory
    ClimaComms.barrier(comms_ctx)
    if ClimaComms.iamroot(comms_ctx)
        rm("./test_checkpoint/", force = true, recursive = true)
    end
end
