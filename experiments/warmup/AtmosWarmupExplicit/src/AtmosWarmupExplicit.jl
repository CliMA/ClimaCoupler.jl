"""
    AtmosWarmupExplicit

Caches the explicit tendency path (`remaining_tendency!`), the largest single component of a cold first step.

Part of the layered atmosphere precompilation workloads in `experiments/warmup`. The
dependency edges between these packages mirror ClimaAtmos's call graph and are *required*:
without them each package re-caches the code it calls, and the resulting package image
exceeds the ~2 GiB limit imposed by 32-bit ELF relocations (see `experiments/warmup/README.md`).

`CLIMACOMMS_DEVICE=CUDA` must be set when precompiling, or the cached specializations are
for CPU-backed fields and a GPU run gains nothing.
"""
module AtmosWarmupExplicit

import ClimaAtmos as CA
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers
import YAML
import AtmosWarmupSetup  # layered: its cached code is reused, not re-serialised
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
            integ = sim.integrator
            Yₜ = similar(integ.u)
            CA.remaining_tendency!(Yₜ, Yₜ, integ.u, integ.p, integ.t)
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
