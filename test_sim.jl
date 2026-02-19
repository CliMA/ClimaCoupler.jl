import ClimaCoupler
import ClimaDiagnostics
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))
config_file = "/home/kphan/Desktop/work/ClimaCoupler.jl/test-calibration/config/ci_configs/amip_albedo_function.yml"
config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)

# setup_and_run(config_dict)
cs = ClimaCoupler.Interfacer.CoupledSimulation(config_dict)

"""
    only_diagnostics_run!(sim, t_end)

Run a `sim`ulation that only compute diagnostics.
"""
function only_diagnostics_run!(sim, t_end)
    (; integrator) = sim
    diagnostics_handler = integrator.callback.discrete_callbacks[end].affect!.diagnostics_handler
    while float(sim.integrator.t) < float(t_end)
        sim.integrator.t += integrator.dt
        ClimaDiagnostics.orchestrate_diagnostics(integrator, diagnostics_handler)
    end
    return nothing
end

function only_diagnostics_run!(cs::ClimaCoupler.Interfacer.CoupledSimulation)
    t_end = cs.tspan[end]
    only_diagnostics_run!(cs.model_sims.atmos_sim, t_end)
end

@time only_diagnostics_run!(cs)

outdir_path = "/home/kphan/Desktop/work/ClimaCoupler.jl/test-calibration/output/amip_albedo_function/output_0018"

Plotting.compute_leaderboard("/home/kphan/Desktop/work/ClimaCoupler.jl/test-calibration/yooo", outdir_path, 3)
