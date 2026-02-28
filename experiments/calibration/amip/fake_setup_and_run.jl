import ClimaCoupler
import ClimaDiagnostics
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))

"""
    only_diagnostics_run!(sim, t_end)

Run a `sim`ulation that only compute diagnostics.
"""
function only_diagnostics_run!(sim, t_end)
    (; integrator) = sim
    diagnostics_handler =
        integrator.callback.discrete_callbacks[end].affect!.diagnostics_handler
    @info integrator.dt
    while sim.integrator.t < t_end
        sim.integrator.t += integrator.dt
        ClimaDiagnostics.orchestrate_diagnostics(integrator, diagnostics_handler)
    end
    return nothing
end

function only_diagnostics_run!(cs::ClimaCoupler.Interfacer.CoupledSimulation)
    t_end = cs.tspan[end]
    only_diagnostics_run!(cs.model_sims.atmos_sim, t_end)
end

function fake_setup_and_run(config_dict)
    cs = ClimaCoupler.Interfacer.CoupledSimulation(config_dict)
    only_diagnostics_run!(cs)
end
