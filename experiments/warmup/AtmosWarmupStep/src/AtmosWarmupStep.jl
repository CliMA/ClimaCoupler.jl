"""
    AtmosWarmupStep

Precompilation workload caching the rest of the first atmosphere step: implicit solve,
Jacobian and timestepper.
"""
module AtmosWarmupStep

import ClimaAtmos as CA
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers
import YAML
# Layered: reuse the code cached by the warmup packages beneath this one.
import AtmosWarmupSetup
import AtmosWarmupExplicit
using PrecompileTools

include(joinpath(@__DIR__, "..", "..", "warmup_utils.jl"))

# Rebuild the image when a config it reads changes, not only when the source does.
include_dependency(config_path("atmos_warmup.yml"))

@setup_workload begin
    outdir = mktempdir(; cleanup = true)
    @compile_workload begin
        status = "not attempted"
        try
            sim = build_sim(outdir)
            ClimaTimeSteppers.step!(sim.integrator)
            status = "ok"
        catch err
            status = "FAILED: " * first(sprint(showerror, err), 300)
        end
        write_report(joinpath(@__DIR__, "..", "workload_report.txt"), status)
    end
end

end # module
