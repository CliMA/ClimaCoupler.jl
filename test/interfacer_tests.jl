import Test: @test, @testset, @test_throws, @test_logs
import ClimaCore as CC
import ClimaParams as CP
import ClimaComms
ClimaComms.@import_required_backends
import Dates
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
import ClimaCoupler: Interfacer

function Interfacer.remap(field, ::Nothing)
    return field
end

# test for a simple generic surface model
struct DummySimulation{S} <: Interfacer.SeaIceModelSimulation
    space::S
end
struct DummySimulation2{S} <: Interfacer.OceanModelSimulation
    space::S
end
struct DummySimulation3{S} <: Interfacer.LandModelSimulation
    space::S
end
struct DummySimulation4{S} <: Interfacer.AtmosModelSimulation
    space::S
end

Interfacer.get_field(sim::Interfacer.SurfaceModelSimulation, ::Val{:var}) = ones(sim.space)
Interfacer.get_field(sim::Interfacer.SurfaceModelSimulation, ::Val{:var_float}) = CC.Spaces.undertype(sim.space)(2)

context = ClimaComms.context()
ClimaComms.init(context)

for FT in (Float32, Float64)
    @testset "test CoupledSim construction, float_type for FT=$FT" begin
        boundary_space =
            CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4, context)

        cs = Interfacer.CoupledSimulation{FT}(
            context,
            nothing, # dates
            boundary_space,
            nothing, # fields
            nothing, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # Δt_cpl
            Ref(Int(0)), # t
            (;), # model_sims
            (;), # callbacks
            (;), # dirs
            nothing, # thermo_params
            nothing, # diags_handler
        )
        @test CC.Spaces.undertype(cs.boundary_space) == FT

    end

    @testset "get_field indexing for FT=$FT" begin
        space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4, context)
        for sim in (DummySimulation(space), DummySimulation2(space), DummySimulation3(space))
            # field
            @test Array(parent(Interfacer.get_field(sim, Val(:var))))[1] == FT(1)
            # float
            @test Interfacer.get_field(sim, Val(:var_float)) == FT(2)
        end
    end

    # test for a simple generic surface model
    @testset "get_field for a SurfaceStub for FT=$FT" begin
        thermo_params = TDP.ThermodynamicsParameters(FT)

        stub = Interfacer.SurfaceStub((;
            area_fraction = FT(1),
            T_sfc = FT(280),
            α_direct = 3,
            α_diffuse = 3,
            z0m = 4,
            z0b = 5,
            beta = 6,
            phase = TD.Liquid(),
            thermo_params = thermo_params,
        ))
        @test Interfacer.get_field(stub, Val(:area_fraction)) == FT(1)
        @test Interfacer.get_field(stub, Val(:surface_temperature)) == FT(280)
        @test Interfacer.get_field(stub, Val(:surface_direct_albedo)) == 3
        @test Interfacer.get_field(stub, Val(:surface_diffuse_albedo)) == 3
        @test Interfacer.get_field(stub, Val(:roughness_momentum)) == 4
        @test Interfacer.get_field(stub, Val(:roughness_buoyancy)) == 5
        @test Interfacer.get_field(stub, Val(:beta)) == 6
    end

    @testset "update_field! the SurfaceStub area_fraction for FT=$FT" begin
        boundary_space =
            CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4, context)

        stub = Interfacer.SurfaceStub((;
            area_fraction = zeros(boundary_space),
            T_sfc = zeros(boundary_space),
            α_direct = zeros(boundary_space),
            α_diffuse = zeros(boundary_space),
            z0m = zeros(boundary_space),
            z0b = zeros(boundary_space),
            beta = zeros(boundary_space),
        ))

        Interfacer.update_field!(stub, Val(:area_fraction), ones(boundary_space))
        Interfacer.update_field!(stub, Val(:surface_temperature), ones(boundary_space) .* 2)
        Interfacer.update_field!(stub, Val(:surface_direct_albedo), ones(boundary_space) .* 3)
        Interfacer.update_field!(stub, Val(:surface_diffuse_albedo), ones(boundary_space) .* 4)

        @test Array(parent(Interfacer.get_field(stub, Val(:area_fraction))))[1] == FT(1)
        @test Array(parent(Interfacer.get_field(stub, Val(:surface_temperature))))[1] == FT(2)
        @test Array(parent(Interfacer.get_field(stub, Val(:surface_direct_albedo))))[1] == FT(3)
        @test Array(parent(Interfacer.get_field(stub, Val(:surface_diffuse_albedo))))[1] == FT(4)
    end
end

@testset "nameof(::SurfaceStub)" begin
    stub = Interfacer.SurfaceStub((;))
    @test nameof(stub) == "SurfaceStub"
end

@testset "undefined get_field for generic val" begin
    FT = Float32
    space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4, context)
    sim = DummySimulation(space)
    val = Val(:v)
    @test_throws ErrorException("undefined field `v` for $(nameof(sim))") Interfacer.get_field(sim, val)
end

@testset "undefined get_field for SurfaceModelSimulation" begin
    FT = Float32
    space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4, context)
    sim = DummySimulation3(space)

    # Test that get_field gives correct warnings for unextended fields
    for value in (
        :area_fraction,
        :roughness_buoyancy,
        :roughness_momentum,
        :surface_direct_albedo,
        :surface_diffuse_albedo,
        :surface_temperature,
    )
        val = Val(value)
        @test_throws ErrorException("undefined field `$value` for $(nameof(sim))") Interfacer.get_field(sim, val)
    end
end

@testset "undefined get_field for AtmosModelSimulation" begin
    FT = Float32
    space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4, context)
    sim = DummySimulation4(space)

    # Test that get_field gives correct warnings for unextended fields
    for value in (
        :height_int,
        :height_sfc,
        :liquid_precipitation,
        :radiative_energy_flux_sfc,
        :radiative_energy_flux_toa,
        :snow_precipitation,
        :u_int,
        :v_int,
    )
        val = Val(value)
        @test_throws ErrorException("undefined field `$value` for $(nameof(sim))") Interfacer.get_field(sim, val)
    end
end

@testset "update_field! warnings for SurfaceModelSimulation" begin
    FT = Float32
    space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4, context)
    dummy_field = CC.Fields.ones(space)
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
        @test_logs (:warn, "`update_field!` is not extended for the `$value` field of $(nameof(sim)): skipping update.") Interfacer.update_field!(
            sim,
            val,
            dummy_field,
        )
        @test_throws ErrorException("undefined field `$value` for $(nameof(sim))") Interfacer.get_field(sim, val)
    end
end

@testset "undefined update_field! warnings for AtmosModelSimulation" begin
    FT = Float32
    space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4, context)
    dummy_field = CC.Fields.ones(space)
    sim = DummySimulation4(space)

    # Test that update_field! gives correct warnings for unextended fields
    for value in (:emissivity, :surface_direct_albedo, :surface_diffuse_albedo, :surface_temperature, :turbulent_fluxes)
        val = Val(value)
        @test_logs (:warn, "`update_field!` is not extended for the `$value` field of $(nameof(sim)): skipping update.") Interfacer.update_field!(
            sim,
            val,
            dummy_field,
        )
        @test_throws ErrorException("undefined field `$value` for $(nameof(sim))") Interfacer.get_field(sim, val)
    end
end

@testset "undefined step! error" begin
    FT = Float32
    sim = DummySimulation3(nothing)
    @test_throws ErrorException("undefined step! for $(nameof(sim))") Interfacer.step!(sim, 1)
end

@testset "SurfaceStub step!" begin
    FT = Float32
    @test isnothing(Interfacer.step!(Interfacer.SurfaceStub(FT(0)), 1))
end

@testset "remap" begin
    FT = Float32
    source_space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4, context)
    field = CC.Fields.coordinate_field(source_space).lat

    # Remap field to target space
    target_space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 6, context)
    field_target_space = Interfacer.remap(field, target_space)

    # remap back to source space
    field_source_space = Interfacer.remap(field_target_space, source_space)

    # The ClimaCore DataLayout underlying the remapped field uses
    # Array instead of SubArray, so we can't compare the fields directly without Copy
    @test field_source_space ≈ copy(field)
end
