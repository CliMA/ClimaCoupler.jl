import Test: @test, @testset
import ClimaCore as CC
import ClimaCoupler

include(joinpath("..", "..", "components", "atmosphere", "climaatmos.jl"))

for FT in (Float32, Float64)
    @testset "dss_state! ClimaAtmosSimulation for FT=$FT" begin
        coupled_param_dict = CP.create_toml_dict(FT)

        boundary_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = coupled_param_dict["planet_radius"], # in meters
            n_quad_points = 4,
            h_elem = 4,
        )

        # set up objects for test
        integrator = (;
            u = (;
                state_field1 = CC.Fields.ones(boundary_space),
                state_field2 = CC.Fields.zeros(boundary_space),
            ),
            p = (; cache_field = CC.Fields.zeros(boundary_space)),
        )
        sim = ClimaAtmosSimulation(nothing, nothing, integrator, nothing)

        # make field non-constant to check the impact of the dss step
        coords_lat = CC.Fields.coordinate_field(sim.integrator.u.state_field2).lat
        @. sim.integrator.u.state_field2 = sin(coords_lat)

        integrator_copy = deepcopy(integrator)
        # apply DSS
        dss_state!(sim)

        # test that uniform field and cache are unchanged, non-constant is changed
        # note: uniform field is changed slightly by dss
        @test sim.integrator.u.state_field1 ≈ integrator_copy.u.state_field1
        @test sim.integrator.u.state_field2 != integrator_copy.u.state_field2
        @test sim.integrator.p.cache_field == integrator_copy.p.cache_field
    end
end
