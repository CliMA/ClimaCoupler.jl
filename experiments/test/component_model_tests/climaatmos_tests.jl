import Test: @test, @testset
import ClimaCore as CC
import ClimaCoupler
import ClimaCoupler.Utilities

import ClimaAtmos
import ClimaParams as CP
import ClimaUtilities.TimeManager: ITime

# Needed to construct ClimaAtmosSimulation
ClimaAtmosExt = Base.get_extension(ClimaCoupler, :ClimaCouplerClimaAtmosExt)
function ClimaCoupler.Utilities.apply_dss!(field::NamedTuple, buffer)
    if isnothing(buffer)
        return nothing
    else
        return map(
            key -> ClimaCoupler.Utilities.apply_dss!(getproperty(field, key), buffer[key]),
            propertynames(field),
        )
    end
end

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
        sim = ClimaAtmosExt.ClimaAtmosSimulation(
            nothing,
            nothing,
            integrator,
            nothing,
            nothing,
        )

        # make field non-constant to check the impact of the dss step
        coords_lat = CC.Fields.coordinate_field(sim.integrator.u.state_field2).lat
        @. sim.integrator.u.state_field2 = sin(coords_lat)

        integrator_copy = deepcopy(integrator)
        # apply DSS
        ClimaAtmosExt.dss_state!(sim)

        # test that uniform field and cache are unchanged, non-constant is changed
        # note: uniform field is changed slightly by dss
        @test sim.integrator.u.state_field1 ≈ integrator_copy.u.state_field1
        @test sim.integrator.u.state_field2 != integrator_copy.u.state_field2
        @test sim.integrator.p.cache_field == integrator_copy.p.cache_field
    end
end

mutable struct FakeAtmosIntegrator{T, P}
    t::T
    dt::T
    p::P
end
# Each step sets the surface rain flux to the new time, and the snow flux to its negative
function ClimaCoupler.Interfacer.step!(integrator::FakeAtmosIntegrator)
    integrator.t += integrator.dt
    integrator.p.precomputed.surface_rain_flux .= float(integrator.t)
    integrator.p.precomputed.surface_snow_flux .= -float(integrator.t)
    return nothing
end

for T in (Float64, ITime)
    @testset "step! ClimaAtmosSimulation averages precipitation over its steps ($T)" begin
        space = CC.CommonSpaces.CubedSphereSpace(
            Float64;
            radius = 1.0,
            n_quad_points = 2,
            h_elem = 2,
        )
        precomputed = (;
            surface_rain_flux = CC.Fields.zeros(space),
            surface_snow_flux = CC.Fields.zeros(space),
        )
        atmos = (; microphysics_model = ClimaAtmos.EquilibriumMicrophysics0M())
        integrator = FakeAtmosIntegrator(T(0), T(1), (; atmos, precomputed))
        precipitation = (; rain = CC.Fields.zeros(space), snow = CC.Fields.zeros(space))
        sim = ClimaAtmosExt.ClimaAtmosSimulation(
            nothing,
            nothing,
            integrator,
            nothing,
            precipitation,
        )
        rain(sim) =
            parent(ClimaCoupler.Interfacer.get_field(sim, Val(:liquid_precipitation)))
        snow(sim) = parent(ClimaCoupler.Interfacer.get_field(sim, Val(:snow_precipitation)))

        ClimaCoupler.Interfacer.step!(sim, T(3)) # fluxes at t = 1, 2, 3
        @test all(≈(2), rain(sim))
        @test all(≈(-2), snow(sim))

        ClimaCoupler.Interfacer.step!(sim, T(5)) # fluxes at t = 4, 5
        @test all(≈(4.5), rain(sim))
        @test all(≈(-4.5), snow(sim))

        ClimaCoupler.Interfacer.step!(sim, T(5)) # no step
        @test all(≈(4.5), rain(sim))
    end
end
