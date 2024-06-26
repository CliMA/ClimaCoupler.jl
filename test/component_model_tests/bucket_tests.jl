import Test: @test, @testset
import ClimaCore as CC
import ClimaCoupler
import ClimaCoupler: TestHelper

include(pkgdir(ClimaCoupler, "experiments/ClimaEarth/components/land/climaland_bucket.jl"))

for FT in (Float32, Float64)
    @testset "dss_state! BucketSimulation for FT=$FT" begin
        # use TestHelper to create space, extract surface space
        subsurface_space = TestHelper.create_space(FT, nz = 2)
        surface_space = CC.Spaces.horizontal_space(subsurface_space)

        # set up objects for test
        dss_buffer_3d = CC.Spaces.create_dss_buffer(CC.Fields.zeros(subsurface_space))
        dss_buffer_2d = CC.Spaces.create_dss_buffer(CC.Fields.zeros(surface_space))

        integrator = (;
            u = CC.Fields.FieldVector(
                state_field1 = CC.Fields.ones(surface_space),
                state_field_2d = CC.Fields.zeros(surface_space),
                state_field_3d = CC.Fields.zeros(subsurface_space),
            ),
            p = (;
                cache_field = CC.Fields.zeros(surface_space),
                dss_buffer_2d = dss_buffer_2d,
                dss_buffer_3d = dss_buffer_3d,
            ),
            t = FT(0),
        )
        sim = BucketSimulation(nothing, nothing, nothing, integrator, nothing)

        # make fields non-constant to check the impact of the dss step
        coords_lat = CC.Fields.coordinate_field(sim.integrator.u.state_field_2d).lat
        @. sim.integrator.u.state_field_2d = sin(coords_lat)

        coords_lat = CC.Fields.coordinate_field(sim.integrator.u.state_field_3d).lat
        @. sim.integrator.u.state_field_3d = sin(coords_lat)

        # apply DSS
        integrator_copy = deepcopy(integrator)
        dss_state!(sim)

        # test that uniform field and cache are unchanged, non-constant is changed
        # note: uniform field is changed slightly by dss
        @test sim.integrator.u.state_field1 ≈ integrator_copy.u.state_field1
        @test sim.integrator.u.state_field_2d != integrator_copy.u.state_field_2d
        @test sim.integrator.u.state_field_3d != integrator_copy.u.state_field_3d
        @test sim.integrator.p.cache_field == integrator_copy.p.cache_field
    end
end
