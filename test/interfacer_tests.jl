using ClimaCore: Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaCoupler: Regridder, TestHelper
using Test
import Thermodynamics as TD
import ClimaParams as CP
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
    AtmosModelSimulation,
    SurfaceStub,
    update_field!,
    reinit!,
    step!,
    update_turbulent_fluxes_point!

# test for a simple generic surface model
struct DummySimulation{S} <: SeaIceModelSimulation
    space::S
end
struct DummySimulation2{S} <: OceanModelSimulation
    space::S
end
struct DummySimulation3{S} <: LandModelSimulation
    space::S
end
name(::DummySimulation3) = "DummySimulation3"
struct DummySimulation4{S} <: AtmosModelSimulation
    space::S
end
name(::DummySimulation4) = "DummySimulation4"

get_field(sim::SurfaceModelSimulation, ::Val{:var}) = ones(sim.space)
function get_field(sim::SurfaceModelSimulation, ::Val{:var_float})
    FT = Domains.float_type(Meshes.domain(sim.space.grid.topology.mesh))
    return FT(2)
end

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
            nothing, # turbulent_fluxes
            nothing, # thermo_params
        )

        @test float_type(cs) == FT
    end

    @testset "get_field indexing for FT=$FT" begin
        space = TestHelper.create_space(FT)
        for sim in (DummySimulation(space), DummySimulation2(space), DummySimulation3(space))
            # field
            colidx = Fields.ColumnIndex{2}((1, 1), 73)
            @test parent(get_field(sim, Val(:var), colidx))[1] == FT(1)
            # float
            @test get_field(sim, Val(:var_float), colidx) == FT(2)
        end
    end

    # test for a simple generic surface model
    @testset "get_field for a SurfaceStub for FT=$FT" begin
        thermo_params = TDP.ThermodynamicsParameters(FT)

        stub = SurfaceStub((;
            area_fraction = FT(1),
            T_sfc = FT(280),
            α_direct = 3,
            α_diffuse = 3,
            z0m = 4,
            z0b = 5,
            beta = 6,
            ρ_sfc = FT(1),
            phase = TD.Liquid(),
            thermo_params = thermo_params,
        ))
        @test get_field(stub, Val(:area_fraction)) == FT(1)
        @test get_field(stub, Val(:surface_temperature)) == FT(280)
        @test get_field(stub, Val(:surface_direct_albedo)) == 3
        @test get_field(stub, Val(:surface_diffuse_albedo)) == 3
        @test get_field(stub, Val(:roughness_momentum)) == 4
        @test get_field(stub, Val(:roughness_buoyancy)) == 5
        @test get_field(stub, Val(:beta)) == 6
        @test get_field(stub, Val(:air_density)) == FT(1)
        @test ≈(get_field(stub, Val(:surface_humidity))[1], FT(0.0076), atol = FT(1e-4))
    end

    @testset "update_field! the SurfaceStub area_fraction for FT=$FT" begin
        boundary_space = TestHelper.create_space(FT)

        stub = SurfaceStub((;
            area_fraction = zeros(boundary_space),
            T_sfc = zeros(boundary_space),
            α_direct = zeros(boundary_space),
            α_diffuse = zeros(boundary_space),
            z0m = zeros(boundary_space),
            z0b = zeros(boundary_space),
            beta = zeros(boundary_space),
        ))

        update_field!(stub, Val(:area_fraction), ones(boundary_space))
        update_field!(stub, Val(:surface_temperature), ones(boundary_space) .* 2)
        update_field!(stub, Val(:surface_direct_albedo), ones(boundary_space) .* 3)
        update_field!(stub, Val(:surface_diffuse_albedo), ones(boundary_space) .* 4)

        @test parent(get_field(stub, Val(:area_fraction)))[1] == FT(1)
        @test parent(get_field(stub, Val(:surface_temperature)))[1] == FT(2)
        @test parent(get_field(stub, Val(:surface_direct_albedo)))[1] == FT(3)
        @test parent(get_field(stub, Val(:surface_diffuse_albedo)))[1] == FT(4)
    end
end

@testset "name(::SurfaceStub)" begin
    stub = SurfaceStub((;))
    @test name(stub) == "SurfaceStub"
end

@testset "undefined get_field for generic val" begin
    FT = Float32
    space = TestHelper.create_space(FT)
    sim = DummySimulation(space)
    val = Val(:v)
    @test_throws ErrorException("undefined field `v` for " * name(sim)) get_field(sim, val)
end

@testset "undefined get_field for SurfaceModelSimulation" begin
    FT = Float32
    space = TestHelper.create_space(FT)
    sim = DummySimulation3(space)

    # Test that get_field gives correct warnings for unextended fields
    for value in (
        :air_density,
        :area_fraction,
        :beta,
        :roughness_buoyancy,
        :roughness_momentum,
        :surface_direct_albedo,
        :surface_diffuse_albedo,
        :surface_humidity,
        :surface_temperature,
    )
        val = Val(value)
        @test_throws ErrorException("undefined field `$value` for " * name(sim)) get_field(sim, val)
    end
end

@testset "undefined get_field for AtmosModelSimulation" begin
    FT = Float32
    space = TestHelper.create_space(FT)
    sim = DummySimulation4(space)

    # Test that get_field gives correct warnings for unextended fields
    for value in (
        :air_density,
        :air_temperature,
        :energy,
        :height_int,
        :height_sfc,
        :liquid_precipitation,
        :radiative_energy_flux_sfc,
        :radiative_energy_flux_toa,
        :snow_precipitation,
        :turbulent_energy_flux,
        :turbulent_moisture_flux,
        :thermo_state_int,
        :uv_int,
        :water,
    )
        val = Val(value)
        @test_throws ErrorException("undefined field `$value` for " * name(sim)) get_field(sim, val)
    end
end

@testset "update_field! warnings for SurfaceModelSimulation" begin
    FT = Float32
    space = TestHelper.create_space(FT)
    dummy_field = Fields.ones(space)
    sim = DummySimulation3(space)

    # Test that update_field! gives correct warnings for unextended fields
    for value in (
        :air_density,
        :area_fraction,
        :liquid_precipitation,
        :radiative_energy_flux_sfc,
        :snow_precipitation,
        :turbulent_energy_flux,
        :turbulent_moisture_flux,
    )
        val = Val(value)
        @test_logs (
            :warn,
            "`update_field!` is not extended for the `$value` field of " * name(sim) * ": skipping update.",
        ) update_field!(sim, val, dummy_field)
        @test_throws ErrorException("undefined field `$value` for " * name(sim)) get_field(sim, val)
    end
end

@testset "undefined update_field! warnings for AtmosModelSimulation" begin
    FT = Float32
    space = TestHelper.create_space(FT)
    dummy_field = Fields.ones(space)
    sim = DummySimulation4(space)

    # Test that update_field! gives correct warnings for unextended fields
    for value in (:co2, :surface_direct_albedo, :surface_diffuse_albedo, :surface_temperature, :turbulent_fluxes)
        val = Val(value)
        @test_logs (
            :warn,
            "`update_field!` is not extended for the `$value` field of " * name(sim) * ": skipping update.",
        ) update_field!(sim, val, dummy_field)
        @test_throws ErrorException("undefined field `$value` for " * name(sim)) get_field(sim, val)
    end
end

@testset "undefined step! error" begin
    FT = Float32
    sim = DummySimulation3(nothing)
    @test_throws ErrorException("undefined step! for " * name(sim)) step!(sim, 1)
end

@testset "undefined reinit! error" begin
    FT = Float32
    sim = DummySimulation3(nothing)
    @test_throws ErrorException("undefined reinit! for " * name(sim)) reinit!(sim)
end

@testset "SurfaceStub step!" begin
    FT = Float32
    @test step!(SurfaceStub(FT(0)), 1) == nothing
end

@testset "SurfaceStub reinit!" begin
    FT = Float32
    @test reinit!(SurfaceStub(FT(0))) == nothing
end

@testset "SurfaceStub update_turbulent_fluxes_point!" begin
    FT = Float32
    colidx = Fields.ColumnIndex{2}((1, 1), 73) # arbitrary index
    @test update_turbulent_fluxes_point!(SurfaceStub(FT(0)), (;), colidx) == nothing
end

# # Test that update_field! gives correct warnings for unextended fields
# for value in (
#     :air_density,
#     :air_temperature,
#     :energy,
#     :height_int,
#     :height_sfc,
#     :liquid_precipitation,
#     :radiative_energy_flux_sfc,
#     :radiative_energy_flux_toa,
#     :snow_precipitation,
#     :turbulent_energy_flux,
#     :turbulent_moisture_flux,
#     :thermo_state_int,
#     :uv_int,
#     :water,
# )
#     val = Val(value)
#     @test_throws ErrorException("undefined field `$value` for " * name(sim)) get_field(sim, val)
# end
