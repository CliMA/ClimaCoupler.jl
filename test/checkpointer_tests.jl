using Test
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC
import ClimaCoupler: Checkpointer, Interfacer
import StaticArrays
import Dates

const FT = Float64
const space_checkpointer = CC.CommonSpaces.CubedSphereSpace(
    FT;
    radius = FT(6.371e6), # in meters
    n_quad_points = 4,
    h_elem = 4,
)

struct DummySimulation{S} <: Interfacer.AbstractAtmosSimulation
    state::S
end
Checkpointer.get_model_prog_state(sim::DummySimulation) = sim.state

@testset "get_model_prog_state" begin
    sim = DummySimulation((; T = ones(space_checkpointer)))
    @test Checkpointer.get_model_prog_state(sim) == sim.state

    sim2 = Interfacer.SurfaceStub([])
    @test Checkpointer.get_model_prog_state(sim2) === nothing
end

@testset "checkpoint_model_state, restart_model_state!" begin
    comms_ctx = ClimaComms.context(ClimaComms.CPUSingleThreaded())
    t = 1
    prev_checkpoint_t = -1
    # old sim run
    sim = DummySimulation(CC.Fields.FieldVector(T = ones(space_checkpointer)))

    dir = mktempdir(; prefix = "test_checkpoint_")
    Checkpointer.checkpoint_model_state(
        sim,
        comms_ctx,
        t,
        prev_checkpoint_t,
        output_dir = dir,
    )

    # new sim run
    sim_new = DummySimulation(CC.Fields.FieldVector(T = zeros(space_checkpointer)))
    input_file = joinpath(dir, "checkpoint_$(nameof(sim_new))_$t.hdf5")
    Checkpointer.restart_model_state!(sim_new, input_file, comms_ctx)
    @test sim_new.state.T == sim.state.T

    # remove checkpoint directory
    rm(dir, force = true, recursive = true)
end

@testset "restore! for different types" begin
    comms_ctx = ClimaComms.context()

    # Test restore! for arrays
    v1 = [1.0, 2.0, 3.0]
    v2 = [4.0, 5.0, 6.0]
    Checkpointer.restore!(v1, v2, comms_ctx)
    @test v1 == v2

    # Test restore! for ClimaCore data layouts
    v1 = CC.Fields.field_values(ones(space_checkpointer))
    v2 = CC.Fields.field_values(zeros(space_checkpointer))
    Checkpointer.restore!(v1, v2, comms_ctx)
    @test v1 == v2

    # Test restore! for StaticArrays
    v1 = StaticArrays.SVector{3, Float64}(1.0, 2.0, 3.0)
    v2 = StaticArrays.SVector{3, Float64}(1.0, 2.0, 3.0)
    Checkpointer.restore!(v1, v2, comms_ctx)
    @test v1 == v2

    v3 = StaticArrays.SVector{3, Float64}(4.0, 5.0, 6.0)
    @test_throws ErrorException Checkpointer.restore!(v1, v3, comms_ctx)

    # Test restore! for Numbers
    v1 = 42
    v2 = 42
    Checkpointer.restore!(v1, v2, comms_ctx)
    @test v1 == v2

    v3 = 43
    @test_throws ErrorException Checkpointer.restore!(v1, v3, comms_ctx)

    # Test restore! for UnitRange
    v1 = 1:5
    v2 = 1:5
    Checkpointer.restore!(v1, v2, comms_ctx)
    @test v1 == v2

    v3 = 1:6
    @test_throws ErrorException Checkpointer.restore!(v1, v3, comms_ctx)

    # Test restore! for Symbol
    v1 = :test
    v2 = :test
    Checkpointer.restore!(v1, v2, comms_ctx)
    @test v1 == v2

    v3 = :other
    @test_throws ErrorException Checkpointer.restore!(v1, v3, comms_ctx)

    # Test restore! for Dict
    v1 = Dict(:a => 1, :b => 2)
    v2 = Dict(:a => 1, :b => 2)
    Checkpointer.restore!(v1, v2, comms_ctx)
    @test v1 == v2

    v3 = Dict(:a => 1, :b => 3)
    @test_throws ErrorException Checkpointer.restore!(v1, v3, comms_ctx)

    # Test restore! for Dates types
    v1 = Dates.DateTime(2000, 1, 1)
    v2 = Dates.DateTime(2000, 1, 1)
    Checkpointer.restore!(v1, v2, comms_ctx)
    @test v1 == v2

    v3 = Dates.DateTime(2000, 1, 2)
    @test_logs (:warn, "Time value differs in restart") Checkpointer.restore!(
        v1,
        v3,
        comms_ctx,
    )

    # Test restore! for Dates.UTInstant
    v1 = Dates.UTInstant(Dates.Minute(0))
    v2 = Dates.UTInstant(Dates.Minute(0))
    Checkpointer.restore!(v1, v3, comms_ctx)
    @test v1 == v2

    v3 = Dates.UTInstant(Dates.Minute(1))
    @test_logs (:warn, "Time value differs in restart") Checkpointer.restore!(
        v1,
        v3,
        comms_ctx,
    )

    # Test restore! for Dates.Millisecond
    v1 = Dates.Millisecond(1000)
    v2 = Dates.Millisecond(1000)
    Checkpointer.restore!(v1, v2, comms_ctx)
    @test v1 == v2

    v3 = Dates.Millisecond(2000)
    @test_logs (:warn, "Time value differs in restart") Checkpointer.restore!(
        v1,
        v3,
        comms_ctx,
    )

    # Test restore! does nothing for comms contexts, data types
    v1 = ClimaComms.context()
    v2 = ClimaComms.context()
    @test isnothing(Checkpointer.restore!(v1, v2, comms_ctx))

    v1 = Float64
    v2 = Float32
    @test isnothing(Checkpointer.restore!(v1, v2, comms_ctx))

    # Test restore! for structs with fields
    struct TestStructError
        a::Float64
    end
    v1 = TestStructError(1.0)
    v2 = TestStructError(5.0)
    @test_throws ErrorException Checkpointer.restore!(v1, v2, comms_ctx)
    @test v1.a != v2.a

    # Test restore! with ignore parameter
    struct TestStructWithIgnore
        a::Ref{Float64}
        b::Int
        c::Vector{Float64}
    end

    v1 = TestStructWithIgnore(Ref(1.0), 2, [3.0, 4.0])
    v2 = TestStructWithIgnore(Ref(1.0), 6, [7.0, 8.0])
    Checkpointer.restore!(v1, v2, comms_ctx; ignore = Set([:b]))
    @test v1.a[] == v2.a[] # Ref doesn't get restored
    @test v1.b != v2.b # Int doesn't get restored
    @test v1.c == v2.c
end

@testset "Field iterator" begin
    struct A
        B::Any
        C::Any
        D::Any
        ignore_this::Any
    end
    cache = A(A([4], [3], [2], [42]), [1], [0], [42])
    itr = Checkpointer.FieldIterator(cache, ignore = Set([:ignore_this]))
    @test !isempty(itr)
    cache_vec = collect(itr)
    @test cache_vec == [[0], [1], [2], [3], [4]]
    @test isempty(itr)
end

@testset "Atmos cache iterator" begin
    struct IgnoreFields
        rc::Any
        params::Any
        ghost_buffer::Any
        hyperdiffusion_ghost_buffer::Any
        data_handler::Any
        graph_context::Any
        dt::Any
    end

    ignore_cache = IgnoreFields([1], [1], [1], [1], [1], [1], [1])
    itr = Checkpointer.CacheIterator(
        cache,
        ignore = Set([
            :rc,
            :params,
            :ghost_buffer,
            :hyperdiffusion_ghost_buffer,
            :data_handler,
            :graph_context,
            :dt,
        ]),
    )
    cache_vec = collect(itr)
    @test isempty(cache_vec)
end
