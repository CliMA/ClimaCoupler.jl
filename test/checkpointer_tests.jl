using ClimaCore: Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaCoupler: TestHelper
using ClimaComms
using Test
import ClimaCoupler: Interfacer
import ClimaCoupler.Checkpointer: get_model_prog_state, restart_model_state!, checkpoint_model_state

FT = Float64

struct DummySimulation{S} <: Interfacer.AtmosModelSimulation
    state::S
end
get_model_prog_state(sim::DummySimulation) = sim.state

@testset "get_model_prog_state" begin
    boundary_space = TestHelper.create_space(FT)
    sim = DummySimulation((; T = ones(boundary_space)))
    @test get_model_prog_state(sim) == sim.state

    sim2 = Interfacer.SurfaceStub([])
    @test get_model_prog_state(sim2) == nothing
end

@testset "checkpoint_model_state, restart_model_state!" begin
    boundary_space = TestHelper.create_space(FT)
    t = 1
    comms_ctx = ClimaComms.context(ClimaComms.CPUSingleThreaded())
    # old sim run
    sim = DummySimulation(Fields.FieldVector(T = ones(boundary_space)))
    checkpoint_model_state(sim, comms_ctx, t, output_dir = "test_checkpoint")

    # new sim run
    sim_new = DummySimulation(Fields.FieldVector(T = zeros(boundary_space)))
    restart_model_state!(sim_new, comms_ctx, t, input_dir = "test_checkpoint")
    @test sim_new.state.T == sim.state.T

    # remove checkpoint directory
    rm("./test_checkpoint/", force = true, recursive = true)
end
