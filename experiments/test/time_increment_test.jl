using Test
using ClimaCoupler
import ClimaAtmos
ClimaAtmosExt = Base.get_extension(ClimaCoupler, :ClimaCouplerClimaAtmosExt)

@testset "Floating Point time-stepping" begin
    TT = Float64
    t = TT(10^7)
    smallest_dt = TT(60)
    # use a fake sim and integrator because their functionality is not needed in this test
    struct FakeIntegrator{TT}
        dt::TT
        t::Ref{TT}
    end
    Base.getproperty(integrator::FakeIntegrator, s::Symbol) =
        s == :t ? getfield(integrator, s)[] : getfield(integrator, s)
    fake_sim = ClimaAtmosExt.ClimaAtmosSimulation(
        nothing,
        nothing,
        FakeIntegrator(smallest_dt, Ref(t)),
        nothing,
    )
    function ClimaCoupler.Interfacer.step!(sim::FakeIntegrator, Δt, _)
        getfield(sim, :t)[] += Δt
        return
    end
    # Take a sufficiently-large number of steps to test robustness against floating point drift
    n_steps = 9999999
    for i in 1:n_steps
        t += smallest_dt
        ClimaCoupler.Interfacer.step!(fake_sim, t)
        @test fake_sim.integrator.t == t
    end
end
