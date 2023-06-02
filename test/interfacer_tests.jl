using ClimaCore: Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaCoupler: Utilities, Regridder, TestHelper
using Test
import ClimaCoupler.Interfacer:
    get_field,
    name,
    SurfaceModelSimulation,
    LandModelSimulation,
    OceanModelSimulation,
    SeaIceModelSimulation,
    SurfaceStub,
    update_field!
FT = Float64

# test for a simple generic surface model
struct DummySimulation <: SeaIceModelSimulation end
struct DummySimulation2 <: OceanModelSimulation end
struct DummySimulation3 <: LandModelSimulation end

boundary_space = TestHelper.create_space(FT);
get_field(::SurfaceModelSimulation, ::Val{:var}) = ones(boundary_space)
get_field(::SurfaceModelSimulation, ::Val{:var_float}) = FT(2)
get_field(::SurfaceModelSimulation, ::Val{:surface_temperature}) = ones(boundary_space) .* FT(300)

@testset "get_field indexing" begin
    for sim in (DummySimulation(), DummySimulation2(), DummySimulation3())
        # field
        colidx = Fields.ColumnIndex{2}((1, 1), 73)
        @test parent(get_field(sim, Val(:var), colidx))[1] == FT(1)
        # float
        @test get_field(sim, Val(:var_float), colidx) == FT(2)
    end
end

@testset "undefined get_field" begin
    sim = DummySimulation()
    val = Val(:v)
    @test_throws ErrorException("undefined field $val for " * name(sim)) get_field(sim, val)
end

# test for a simple generic surface model
@testset "get_field for a SurfaceStub" begin
    stub = SurfaceStub((; area_fraction = 1, T_sfc = 2, α = 3, z0m = 4, z0b = 5, beta = 6))
    @test get_field(stub, Val(:area_fraction)) == 1
    @test get_field(stub, Val(:surface_temperature)) == 2
    @test get_field(stub, Val(:albedo)) == 3
    @test get_field(stub, Val(:roughness_momentum)) == 4
    @test get_field(stub, Val(:roughness_buoyancy)) == 5
    @test get_field(stub, Val(:beta)) == 6
end

@testset "update_field! the SurfaceStub area_fraction" begin
    boundary_space = TestHelper.create_space(FT)

    stub = SurfaceStub((; area_fraction = zeros(boundary_space), T_sfc = zeros(boundary_space)))

    update_field!(stub, Val(:area_fraction), ones(boundary_space))
    update_field!(stub, Val(:surface_temperature), ones(boundary_space) .* 2)

    @test parent(get_field(stub, Val(:area_fraction)))[1] == FT(1)
    @test parent(get_field(stub, Val(:surface_temperature)))[1] == FT(2)
end
