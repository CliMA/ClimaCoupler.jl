# testing functions used to produce user-defined debugging plots for AMIP experiments
import Test: @test, @testset, @test_logs
import ClimaCore as CC
import ClimaCoupler: Interfacer
import ClimaComms
ClimaComms.@import_required_backends

# Prevent GKS headless operation mode warning
ENV["GKSwstype"] = "nul"
FT = Float64

struct ClimaAtmosSimulation{C} <: Interfacer.AtmosModelSimulation
    cache::C
end
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:atmos_field}) = sim.cache.atmos_field

struct BucketSimulation{C} <: Interfacer.SurfaceModelSimulation
    cache::C
end

struct ClimaLandSimulation{C} <: Interfacer.SurfaceModelSimulation
    cache::C
end

include("../user_io/debug_plots.jl")

Interfacer.get_field(sim::BucketSimulation, ::Val{:surface_field}) = sim.cache.surface_field
Interfacer.get_field(sim::ClimaLandSimulation, ::Val{:surface_field}) =
    sim.cache.surface_field
Interfacer.get_field(sim::Interfacer.SurfaceStub, ::Val{:stub_field}) = sim.cache.stub_field

plot_field_names(sim::ClimaAtmosSimulation) = (:atmos_field,)
plot_field_names(sim::BucketSimulation) = (:surface_field,)
plot_field_names(sim::ClimaLandSimulation) = (:surface_field,)
plot_field_names(sim::Interfacer.SurfaceStub) = (:stub_field,)

@testset "import_atmos_fields!" begin

    boundary_space = CC.CommonSpaces.CubedSphereSpace(
        FT;
        radius = FT(6371e3),
        n_quad_points = 4,
        h_elem = 4,
    )
    coupler_names = [
        :surface_direct_albedo,
        :surface_diffuse_albedo,
        :F_radiative,
        :F_lh,
        :F_sh,
        :F_turb_moisture,
        :F_turb_ρτxz,
        :F_turb_ρτyz,
        :P_liq,
        :P_snow,
        :T_sfc,
    ]
    atmos_names = (:atmos_field,)
    surface_names = (:surface_field,)
    stub_names = (:stub_field,)

    atmos_fields = NamedTuple{atmos_names}(
        ntuple(i -> CC.Fields.zeros(boundary_space), length(atmos_names)),
    )
    surface_fields = NamedTuple{surface_names}(
        ntuple(i -> CC.Fields.zeros(boundary_space), length(surface_names)),
    )
    stub_fields = NamedTuple{stub_names}(
        ntuple(i -> CC.Fields.zeros(boundary_space), length(stub_names)),
    )
    coupler_fields = Interfacer.init_coupler_fields(FT, coupler_names, boundary_space)

    model_sims = (;
        atmos_sim = ClimaAtmosSimulation(atmos_fields),
        surface_sim = BucketSimulation(surface_fields),
        surface_sim2 = ClimaLandSimulation(surface_fields),
        ice_sim = Interfacer.SurfaceStub(stub_fields),
    )
    cs = Interfacer.CoupledSimulation{FT}(
        nothing, # dates
        coupler_fields, # fields
        nothing, # conservation_checks
        (Int(0), Int(1)), # tspan
        Int(200), # Δt_cpl
        Ref(Int(0)), # t
        Ref(-1), # prev_checkpoint_t
        model_sims, # model_sims
        (;), # callbacks
        (;), # dirs
        nothing, # thermo_params
        nothing, # diags_handler
    )

    output_plots = "test_debug"
    @test_logs (:info, "plotting debug in test_debug") match_mode = :any debug(
        cs,
        output_plots,
    )
    @test isfile("test_debug/debug_ClimaAtmosSimulation.png")
    @test isfile("test_debug/debug_BucketSimulation.png")
    @test isfile("test_debug/debug_ClimaLandSimulation.png")
    @test isfile("test_debug/debug_SurfaceStub.png")
    @test isfile("test_debug/debug_coupler.png")

    # remove output
    rm(output_plots; recursive = true)
end
