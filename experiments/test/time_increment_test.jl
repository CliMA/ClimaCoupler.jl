using Test
using ClimaCoupler
import ClimaAtmos
ClimaAtmosExt = Base.get_extension(ClimaCoupler, :ClimaCouplerClimaAtmosExt)

@testset "Floating Point time-stepping" begin
    FT = Float64
    t = FT(10^7)
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
        # component models can use float32 time
        FakeIntegrator(Float32(smallest_dt), Ref(Float32(t))),
        nothing,
    )
    function ClimaCoupler.Interfacer.step!(sim::FakeIntegrator, Δt, _)
        getfield(sim, :t)[] += Δt
        return
    end
    # 951815 is the number of steps after which the time incrementing test fails
    for i in 1:951815
        t += smallest_dt
        ClimaCoupler.Interfacer.step!(fake_sim, t)
        @test fake_sim.integrator.t == t broken = i >= 951815
    end
end
