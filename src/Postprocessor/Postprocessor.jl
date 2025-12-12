module Postprocessor

import ClimaComms
import ClimaAnalysis as CAN
import ClimaDiagnostics as CD

# Include submodules
include("coupler_diagnostics.jl")

export simulated_years_per_day, walltime_per_coupling_step, save_sypd_walltime_to_disk

"""
    simulated_years_per_day(cs, walltime)

Compute the simulated years per walltime day for the given coupled simulation
`cs`, assuming that the simulation took `walltime` seconds.
"""
function simulated_years_per_day(cs, walltime)
    simulated_seconds_per_second = float(cs.tspan[end] - cs.tspan[begin]) / walltime
    return simulated_seconds_per_second / 365.25
end

"""
    walltime_per_coupling_step(cs, walltime)

Compute the average walltime needed to take one step for the given coupled simulation `cs`,
assuming that the simulation took `walltime` seconds. The result is in seconds.
"""
function walltime_per_coupling_step(cs, walltime)
    n_coupling_steps = (cs.tspan[end] - cs.tspan[begin]) / cs.Δt_cpl
    return walltime / n_coupling_steps
end

"""
    save_sypd_walltime_to_disk(cs, walltime)

Save the computed `sypd`, `walltime_per_coupling_step`,
and memory usage to text files in the `artifacts` directory.
"""
function save_sypd_walltime_to_disk(cs, walltime)
    if ClimaComms.iamroot(ClimaComms.context(cs))
        sypd = simulated_years_per_day(cs, walltime)
        walltime_per_step = walltime_per_coupling_step(cs, walltime)

        open(joinpath(cs.dir_paths.artifacts_dir, "sypd.txt"), "w") do sypd_filename
            println(sypd_filename, "$sypd")
        end

        open(
            joinpath(cs.dir_paths.artifacts_dir, "walltime_per_step.txt"),
            "w",
        ) do walltime_per_step_filename
            println(walltime_per_step_filename, "$(walltime_per_step)")
        end
    end
    return nothing
end

end # module Postprocessor
