# testing functions used to produce user-defined debugging plots for AMIP experiments

using Test
using ClimaCore
using ClimaCoupler: TestHelper
import ClimaCoupler.Interfacer:
    update_field!, AtmosModelSimulation, SurfaceModelSimulation, SurfaceStub, get_field, update_field!, name
using ClimaCoupler.Utilities: CoupledSimulation, CoupledSimulation

FT = Float64

struct ClimaAtmosSimulation{C} <: AtmosModelSimulation
    cache::C
end
name(sim::ClimaAtmosSimulation) = "ClimaAtmosSimulation"
get_field(sim::AtmosModelSimulation, ::Val{:atmos_field}) = sim.cache.atmos_field

struct BucketSimulation{C} <: SurfaceModelSimulation
    cache::C
end
name(sim::BucketSimulation) = "BucketSimulation"

include("../../experiments/AMIP/modular/user_io/debug_plots.jl")

get_field(sim::BucketSimulation, ::Val{:surface_field}) = sim.cache.surface_field
get_field(sim::SurfaceStub, ::Val{:stub_field}) = sim.cache.stub_field

plot_field_names(sim::ClimaAtmosSimulation) = (:atmos_field,)
plot_field_names(sim::BucketSimulation) = (:surface_field,)
plot_field_names(sim::SurfaceStub) = (:stub_field,)

@testset "import_atmos_fields!" begin

    boundary_space = TestHelper.create_space(FT)
    coupler_names = (:F_turb_energy, :F_turb_moisture, :P_liq, :T_S, :ρ_sfc, :q_sfc)
    atmos_names = (:atmos_field,)
    surface_names = (:surface_field,)
    stub_names = (:stub_field,)

    atmos_fields = NamedTuple{atmos_names}(ntuple(i -> ClimaCore.Fields.ones(boundary_space), length(atmos_names)))
    surface_fields =
        NamedTuple{surface_names}(ntuple(i -> ClimaCore.Fields.ones(boundary_space), length(surface_names)))
    stub_fields = NamedTuple{stub_names}(ntuple(i -> ClimaCore.Fields.ones(boundary_space), length(stub_names)))
    coupler_fields =
        NamedTuple{coupler_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_names)))

    model_sims = (;
        atmos_sim = ClimaAtmosSimulation(atmos_fields),
        surface_sim = BucketSimulation(surface_fields),
        ice_sim = SurfaceStub(stub_fields),
    )
    cs = CoupledSimulation{FT}(
        nothing, # comms_ctx
        nothing, # dates
        nothing, # boundary_space
        coupler_fields, # fields
        nothing, # parsed_args
        nothing, # conservation_checks
        (Int(0), Int(1)), # tspan
        Int(200), # t
        Int(200), # Δt_cpl
        (;), # surface_masks
        model_sims, # model_sims
        (;), # mode
        (), # diagnostics
    )

    output_plots = "test_debug"
    debug(cs, output_plots)
    @test isfile("test_debug/debug_ClimaAtmosSimulation.png")
    @test isfile("test_debug/debug_BucketSimulation.png")
    @test isfile("test_debug/debug_SurfaceStub.png")
    @test isfile("test_debug/debug_coupler.png")

    # remove output
    rm(output_plots; recursive = true)

end
