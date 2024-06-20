import Test: @test, @testset
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaCore as CC
import ClimaCoupler: Checkpointer, Interfacer, TestHelper

FT = Float64

struct DummySimulation{S} <: Interfacer.AtmosModelSimulation
    state::S
end
Checkpointer.get_model_prog_state(sim::DummySimulation) = sim.state

@testset "get_model_prog_state" begin
    boundary_space = TestHelper.create_space(FT)
    sim = DummySimulation((; T = ones(boundary_space)))
    @test Checkpointer.get_model_prog_state(sim) == sim.state

    sim2 = Interfacer.SurfaceStub([])
    @test Checkpointer.get_model_prog_state(sim2) === nothing
end

@testset "checkpoint_model_state, restart_model_state!" begin
    comms_ctx = ClimaComms.context(ClimaComms.CPUSingleThreaded())
    boundary_space = TestHelper.create_space(FT; comms_ctx = comms_ctx)
    t = 1
    # old sim run
    sim = DummySimulation(CC.Fields.FieldVector(T = ones(boundary_space)))
    Checkpointer.checkpoint_model_state(sim, comms_ctx, t, output_dir = "test_checkpoint")

    # new sim run
    sim_new = DummySimulation(CC.Fields.FieldVector(T = zeros(boundary_space)))
    Checkpointer.restart_model_state!(sim_new, comms_ctx, t, input_dir = "test_checkpoint")
    @test sim_new.state.T == sim.state.T

    # remove checkpoint directory
    rm("./test_checkpoint/", force = true, recursive = true)
end
