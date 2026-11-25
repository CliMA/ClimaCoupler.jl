"""
    AtmosWarmupExplicit

Precompilation workload caching the explicit tendencies, `ClimaAtmos.remaining_tendency!`.
"""
module AtmosWarmupExplicit

import ClimaAtmos as CA
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers
import YAML
# Lower layers of the warmup stack, reused rather than re-cached.
import AtmosWarmupSetup
using PrecompileTools

include(joinpath(@__DIR__, "..", "..", "warmup_utils.jl"))

# Rebuild when the config changes, not only the source.
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
            status = "FAILED: " * first(sprint(showerror, err), 300)
        end
        write_report(joinpath(@__DIR__, "..", "workload_report.txt"), status)
    end
end

end # module
