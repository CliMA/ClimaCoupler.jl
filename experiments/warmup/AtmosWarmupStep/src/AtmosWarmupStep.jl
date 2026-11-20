"""
    AtmosWarmupStep

Precompilation workload that caches the rest of the first atmosphere step: the implicit
solve, Jacobian and timestepper.

Depends on `AtmosWarmupSetup` and `AtmosWarmupExplicit`, so it reuses their images rather
than caching that code again.

Precompile with `CLIMACOMMS_DEVICE=CUDA`, or the image holds CPU code. See
`experiments/warmup/README.md`.
"""
module AtmosWarmupStep

import ClimaAtmos as CA
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers
import YAML
import AtmosWarmupSetup  # layered: its cached code is reused, not re-serialised
import AtmosWarmupExplicit  # layered: its cached code is reused, not re-serialised
using PrecompileTools

include(joinpath(@__DIR__, "..", "..", "warmup_utils.jl"))

# Julia invalidates a package image when its *source* changes, not when data it read during
# precompilation changes. Register the config so editing it forces a rebuild.
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
            # A broken workload must not break precompilation; it only means startup
            # stays slow.
            status = "FAILED: " * first(sprint(showerror, err), 300)
        end
        write_report(joinpath(@__DIR__, "..", "workload_report.txt"), status)
    end
end

end # module
