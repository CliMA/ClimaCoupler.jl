import Test: @test, @testset
import ClimaCore as CC
import ClimaCoupler

include(joinpath("..", "..", "components", "ocean", "slab_ocean.jl"))

for FT in (Float32, Float64)
    @testset "dss_state! SlabOceanSimulation for FT=$FT" begin
        boundary_space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4)

        # construct dss buffer to put in cache
        dss_buffer = CC.Spaces.create_dss_buffer(CC.Fields.zeros(boundary_space))

        # set up objects for test
        u = CC.Fields.FieldVector(;
            state_field1 = CC.Fields.ones(boundary_space),
            state_field2 = CC.Fields.zeros(boundary_space),
        )
        p = (; cache_field = CC.Fields.zeros(boundary_space), dss_buffer = CC.Spaces.create_dss_buffer(u))
        integrator = (; u, p)
        sim = SlabOceanSimulation(nothing, integrator)

        # make field non-constant to check the impact of the dss step
        coords_lat = CC.Fields.coordinate_field(sim.integrator.u.state_field2).lat
        @. sim.integrator.u.state_field2 = sin(coords_lat)

        # apply DSS
        integrator_copy = deepcopy(integrator)
        dss_state!(sim)

        # test that uniform field and cache are unchanged, non-constant is changed
        # note: uniform field is changed slightly by dss
        @test sim.integrator.u.state_field1 ≈ integrator_copy.u.state_field1
        @test sim.integrator.u.state_field2 != integrator_copy.u.state_field2
        @test sim.integrator.p.cache_field == integrator_copy.p.cache_field
    end
end
