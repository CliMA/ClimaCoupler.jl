"""
    AtmosWarmupExplicit

Precompilation workload that caches the explicit tendency path,
`ClimaAtmos.remaining_tendency!`.

Depends on `AtmosWarmupSetup`, so it reuses that image rather than caching setup again.

Precompile with `CLIMACOMMS_DEVICE=CUDA`, or the image holds CPU code. See
`experiments/warmup/README.md`.
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
