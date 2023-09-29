using Test
import ClimaCoupler
using ClimaCoupler.TestHelper: create_space
using ClimaCore: Fields, Spaces

include(pkgdir(ClimaCoupler, "experiments/AMIP/modular/components/land/bucket_init.jl"))
include(pkgdir(ClimaCoupler, "experiments/AMIP/modular/components/land/bucket_utils.jl"))

# struct DummySimulationBucket{I} <: BucketSimulation
#     integrator::I
# end

# TODO bucket doesn't currently work with Float32, but we want to eventually test with both FTs
for FT in (Float64,)
    @testset "dss_state! BucketSimulation for FT=$FT" begin
        # use TestHelper to create space, extract surface space
        subsurface_space = create_space(FT, nz = 2)
        surface_space = subsurface_space.horizontal_space

        # set up objects for test
        dss_buffer_3d = Spaces.create_dss_buffer(Fields.zeros(subsurface_space))
        dss_buffer_2d = Spaces.create_dss_buffer(Fields.zeros(surface_space))

        integrator = (;
            u = Fields.FieldVector(
                state_field1 = Fields.ones(surface_space),
                state_field_2d = Fields.zeros(surface_space),
                state_field_3d = Fields.zeros(subsurface_space),
            ),
            p = (;
                cache_field = Fields.zeros(surface_space),
                dss_buffer_2d = dss_buffer_2d,
                dss_buffer_3d = dss_buffer_3d,
            ),
            t = FT(0),
        )
        integrator_copy = deepcopy(integrator)
        # sim = DummySimulationBucket(integrator)
        sim = BucketSimulation(nothing, nothing, nothing, integrator, nothing)

        # make fields non-constant to check the impact of the dss step
        for i in eachindex(parent(sim.integrator.u.state_field_2d))
            parent(sim.integrator.u.state_field_2d)[i] = sin(i)
        end
        for i in eachindex(parent(sim.integrator.u.state_field_3d))
            parent(sim.integrator.u.state_field_3d)[i] = sin(i)
        end

        # apply DSS
        dss_state!(sim)

        # test that uniform field and cache are unchanged, non-constant is changed
        # note: uniform field is changed slightly by dss
        @test sim.integrator.u.state_field1 ≈ integrator_copy.u.state_field1
        @test sim.integrator.u.state_field_2d != integrator_copy.u.state_field_2d
        @test sim.integrator.u.state_field_3d != integrator_copy.u.state_field_3d
        @test sim.integrator.p.cache_field == integrator_copy.p.cache_field
    end
end
