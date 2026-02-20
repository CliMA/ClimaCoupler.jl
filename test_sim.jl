import ClimaCoupler
import ClimaDiagnostics
import ClimaCoupler: Plotting
include(joinpath(pkgdir(ClimaCoupler), "experiments", "ClimaEarth", "setup_run.jl"))
config_file = "config/ci_configs/amip_albedo_function.yml"
config_dict = ClimaCoupler.Input.get_coupler_config_dict(config_file)

cs = ClimaCoupler.Interfacer.CoupledSimulation(config_dict)

"""
    only_diagnostics_run!(sim, t_end)

Run a `sim`ulation that only compute diagnostics.
"""
function only_diagnostics_run!(sim, t_end)
    (; integrator) = sim
    diagnostics_handler = integrator.callback.discrete_callbacks[end].affect!.diagnostics_handler
    @info integrator.dt
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

@info "Beginning diagnostics only run"
@time only_diagnostics_run!(cs)

outdir_path = "output/amip_albedo_function/output_active"
leaderboard_path = "leaderboard_test/"

Plotting.compute_leaderboard(leaderboard_path, outdir_path, 3)
Plotting.compute_pfull_leaderboard(
    leaderboard_path,
    outdir_path,
    3,
)
