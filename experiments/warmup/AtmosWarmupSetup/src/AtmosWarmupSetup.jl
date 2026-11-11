"""
    AtmosWarmupSetup

Caches the atmosphere setup path (`ClimaAtmos.get_simulation`).

Part of the layered atmosphere precompilation workloads in `experiments/warmup`. The
dependency edges between these packages mirror ClimaAtmos's call graph and are *required*:
without them each package re-caches the code it calls, and the resulting package image
exceeds the ~2 GiB limit imposed by 32-bit ELF relocations (see `experiments/warmup/README.md`).

`CLIMACOMMS_DEVICE=CUDA` must be set when precompiling, or the cached specializations are
for CPU-backed fields and a GPU run gains nothing.
"""
module AtmosWarmupSetup

import ClimaAtmos as CA
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers
import YAML
using PrecompileTools

# Production loads the Makie stack before it builds anything -- `experiments/CMIP/code_loading.jl`
# does so to trigger `ClimaCouplerMakieExt`. Loading it *after* this workload's code has been
# compiled invalidates that code, so the cache has to be built in the same world the run uses.
# Measured: with the stack absent from the workload, `ClimaAtmos.get_simulation` compiled 0.2 s
# in a session without Makie and 227.1 s in one with it, and coupled setup went 4.4 s -> 1153.1 s.
# These imports are load-bearing even though nothing here calls them.
import CairoMakie, GeoMakie, Makie, Poppler_jll, Printf

# The rest of what `experiments/CMIP/code_loading.jl` loads. Like the Makie stack above,
# these are load-bearing even though nothing here calls them: a package loaded *after* this
# workload's code is compiled invalidates it, and GPU kernels are attached to the
# invalidated CodeInstances, so they are thrown away with it.
#
# Measured, the same atmosphere workload against this image: step 1 linked its 260 stepping
# kernels in 0.9 s with only the atmosphere chain loaded, and recompiled all 260 at a cost
# of 128.9 s once the coupled packages were loaded too. The cache has to be built in the
# world a production run executes in.
import ClimaCoupler
import Oceananigans, ClimaSeaIce, KernelAbstractions, ConservativeRegridding, Adapt
import ClimaLand


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
            @assert sim.integrator !== nothing
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
