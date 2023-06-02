
using ClimaCore: Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaCoupler: Utilities, Regridder, TestHelper
using Test
import ClimaCoupler.FluxCalculator:
    compute_atmos_turbulent_fluxes!,
    compute_combined_turbulent_fluxes!,
    CombinedAtmosGrid,
    PartitionedComponentModelGrid
import ClimaCoupler.Interfacer: AtmosModelSimulation
FT = Float64

# test for a simple generic atmos model
struct DummySimulation{S, C} <: AtmosModelSimulation
    state::S
    cache::C
end
function compute_atmos_turbulent_fluxes!(sim::DummySimulation, csf)
    sim.cache.flux .= (csf.T_sfc .- sim.state.T) .* sim.cache.κ ./ sim.cache.dz # Eq. 1
end
struct DummySimulation2{C} <: AtmosModelSimulation
    cache::C
end

@testset "compute_combined_turbulent_fluxes!" begin
    boundary_space = TestHelper.create_space(FT)
    coupler_fields = (; T_sfc = 310 .* ones(boundary_space))
    sim =
        DummySimulation((; T = 300 .* ones(boundary_space)), (; κ = FT(0.01), dz = FT(1), flux = zeros(boundary_space)))
    model_sims = (; atmos_sim = sim)
    flux_types = (CombinedAtmosGrid(), PartitionedComponentModelGrid())
    # the result of Eq 1, given the states above is 0.1 W/m2, but under PartitionedComponentModelGrid() turbulent fluxes are
    # not calculated using this method (using combined surface properties), so the fluxes remain 0.
    results = [FT(0.1), FT(0.0)]
    for (i, t) in enumerate(flux_types)
        sim.cache.flux .= FT(0)
        compute_combined_turbulent_fluxes!(model_sims, coupler_fields, t)
        @test parent(sim.cache.flux)[1] == results[i]
    end
    sim2 = DummySimulation2((; cache = (; flux = zeros(boundary_space))))
    model_sims = (; atmos_sim = sim2)
    @test_throws ErrorException compute_atmos_turbulent_fluxes!(sim2, coupler_fields)

end
