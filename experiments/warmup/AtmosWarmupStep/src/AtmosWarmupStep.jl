"""
    AtmosWarmupStep

Caches the remaining first-step path: the implicit solve, Jacobian, and timestepper.

Part of the layered atmosphere precompilation workloads in `experiments/warmup`. The
dependency edges between these packages mirror ClimaAtmos's call graph and are *required*:
without them each package re-caches the code it calls, and the resulting package image
exceeds the ~2 GiB limit imposed by 32-bit ELF relocations (see `experiments/warmup/README.md`).

`CLIMACOMMS_DEVICE=CUDA` must be set when precompiling, or the cached specializations are
for CPU-backed fields and a GPU run gains nothing.
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

# experiments/warmup/<pkg>/src -> experiments/warmup
const WARMUP_DIR = normpath(joinpath(@__DIR__, "..", ".."))
# experiments/warmup -> repository root, for parameter files referenced by the config
const REPO_DIR = normpath(joinpath(WARMUP_DIR, "..", ".."))
const CONFIG = joinpath(WARMUP_DIR, "configs", "atmos_warmup.yml")

# Julia invalidates a package image when its *source* changes, not when data it read during
# precompilation changes. Register the config so editing it forces a rebuild.
include_dependency(CONFIG)

"Written so a caught workload failure is visible rather than silently yielding no speedup."
const REPORT = joinpath(@__DIR__, "..", "workload_report.txt")

"""
    build_sim(outdir)

Build the atmosphere simulation a coupled run would build.

Resolution is deliberately small: `h_elem`/`z_elem` change array sizes, not types. The model
*selections* must match the production config or the cached specializations are unusable --
see `ClimaCoupler.SimCoordinator` and the overrides replicated in the config.
"""
function build_sim(outdir)
    cfg = YAML.load_file(CONFIG)
    cfg["output_dir"] = outdir
    # Parameter files are given relative to the repository root.
    haskey(cfg, "toml") &&
        (cfg["toml"] = [isabspath(f) ? f : joinpath(REPO_DIR, f) for f in cfg["toml"]])
    return CA.get_simulation(CA.AtmosConfig(cfg; job_id = "atmos_warmup"))
end

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
        try
            open(REPORT, "w") do io
                println(io, "status = ", status)
                println(io, "device = ", get(ENV, "CLIMACOMMS_DEVICE", "(unset)"))
            end
        catch
        end
    end
end

end # module
