using ClimaCore: Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaCoupler: Regridder, TestHelper
using Test
import Thermodynamics as TD
import CLIMAParameters as CP
import Thermodynamics.Parameters as TDP
import ClimaCoupler.Interfacer:
    CoupledSimulation,
    float_type,
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

for FT in (Float32, Float64)
    @testset "test CoupledSim construction, float_type for FT=$FT" begin
        cs = CoupledSimulation{FT}(
            nothing, # comms_ctx
            nothing, # dates
            nothing, # boundary_space
            nothing, # fields
            nothing, # parsed_args
            nothing, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # t
            Int(200), # Δt_cpl
            (;), # surface_masks
            (;), # model_sims
            (;), # mode
            (), # diagnostics
            (;), # callbacks
            (;), # dirs
        )

        @test float_type(cs) == FT
    end

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

        toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
        aliases = string.(fieldnames(TDP.ThermodynamicsParameters))
        param_pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
        thermo_params = TDP.ThermodynamicsParameters{FT}(; param_pairs...)

        stub = SurfaceStub((;
            area_fraction = FT(1),
            T_sfc = FT(280),
            α = 3,
            z0m = 4,
            z0b = 5,
            beta = 6,
            ρ_sfc = FT(1),
            phase = TD.Liquid(),
            thermo_params = thermo_params,
        ))
        @test get_field(stub, Val(:area_fraction)) == FT(1)
        @test get_field(stub, Val(:surface_temperature)) == FT(280)
        @test get_field(stub, Val(:albedo)) == 3
        @test get_field(stub, Val(:roughness_momentum)) == 4
        @test get_field(stub, Val(:roughness_buoyancy)) == 5
        @test get_field(stub, Val(:beta)) == 6
        @test ≈(get_field(stub, Val(:surface_humidity))[1], FT(0.0076), atol = FT(1e-4))
    end

    @testset "name(::SurfaceStub)" begin
        stub = SurfaceStub((;))
        @test name(stub) == "SurfaceStub"
    end

    @testset "update_field! the SurfaceStub area_fraction" begin
        boundary_space = TestHelper.create_space(FT)

        stub = SurfaceStub((; area_fraction = zeros(boundary_space), T_sfc = zeros(boundary_space)))

        update_field!(stub, Val(:area_fraction), ones(boundary_space))
        update_field!(stub, Val(:surface_temperature), ones(boundary_space) .* 2)

        @test parent(get_field(stub, Val(:area_fraction)))[1] == FT(1)
        @test parent(get_field(stub, Val(:surface_temperature)))[1] == FT(2)
    end
end
