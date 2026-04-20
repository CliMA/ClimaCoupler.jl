using Test
using ClimaCoupler
import ClimaAtmos
ClimaAtmosExt = Base.get_extension(ClimaCoupler, :ClimaCouplerClimaAtmosExt)

@testset "Float64 time" begin
    FT = Float64
    t = FT(10^17)
    smallest_dt = FT(60)
    # use a fake sim and integrator because their functionality is not needed in this test
    struct FakeIntegrator{FT}
        dt::FT
        t::Ref{FT}
    end
    Base.getproperty(integrator::FakeIntegrator, s::Symbol) =
        s == :t ? getfield(integrator, s)[] : getfield(integrator, s)
    fake_sim = ClimaAtmosExt.ClimaAtmosSimulation(
        nothing,
        nothing,
        FakeIntegrator(6 * smallest_dt, Ref(t)),
        nothing,
    )
    function ClimaCoupler.Interfacer.step!(sim::FakeIntegrator, _, _)
        getfield(sim, :t)[] = t
        return
    end
    for i in 1:(10^5)
        t += smallest_dt
        ClimaCoupler.Interfacer.step!(fake_sim, t)
        if i % 6 == 0
            @test fake_sim.integrator.t == t
        else
            @test fake_sim.integrator.t < t
        end
    end
end
