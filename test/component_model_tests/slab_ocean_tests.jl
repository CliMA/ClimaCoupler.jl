using Test
import ClimaCoupler
using ClimaCoupler.Interfacer: OceanModelSimulation
using ClimaCoupler.TestHelper: create_space
using ClimaCore
using ClimaCore: Fields, Spaces
import CLIMAParameters as CP
import Thermodynamics.Parameters as TDP

include(pkgdir(ClimaCoupler, "experiments/AMIP/components/ocean/slab_ocean.jl"))

for FT in (Float32, Float64)
    @testset "dss_state! SlabOceanSimulation for FT=$FT" begin
        # use TestHelper to create space
        boundary_space = create_space(FT)

        # construct dss buffer to put in cache
        dss_buffer = Spaces.create_dss_buffer(Fields.zeros(boundary_space))

        # set up objects for test
        integrator = (;
            u = (; state_field1 = FT.(Fields.ones(boundary_space)), state_field2 = FT.(Fields.zeros(boundary_space))),
            p = (; cache_field = FT.(Fields.zeros(boundary_space)), dss_buffer = dss_buffer),
        )
        integrator_copy = deepcopy(integrator)
        sim = SlabOceanSimulation(nothing, nothing, nothing, integrator)

        # make field non-constant to check the impact of the dss step
        for i in eachindex(parent(sim.integrator.u.state_field2))
            parent(sim.integrator.u.state_field2)[i] = FT(sin(i))
        end

        # apply DSS
        dss_state!(sim)

        # test that uniform field and cache are unchanged, non-constant is changed
        # note: uniform field is changed slightly by dss
        @test sim.integrator.u.state_field1 ≈ integrator_copy.u.state_field1
        @test sim.integrator.u.state_field2 != integrator_copy.u.state_field2
        @test sim.integrator.p.cache_field == integrator_copy.p.cache_field
    end
end
