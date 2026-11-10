"""
    CMIPWarmup

Thin top of the warmup DAG: the coupled setup and one coupling step.

Every component's specializations are already cached beneath it, so its incremental
content is the coupler itself -- field exchange, flux calculation, the exchange grid.

Part of the layered precompilation workloads in `experiments/warmup`. The dependency edges
mirror the call graph and are required, not an optimisation: without them a package
re-caches the code it calls and its image exceeds the ~2 GiB limit imposed by 32-bit ELF
relocations. See `experiments/warmup/README.md`.

`CLIMACOMMS_DEVICE=CUDA` must be set when precompiling.
"""
module CMIPWarmup

import ClimaCoupler
import ClimaCoupler: Interfacer, Input, Utilities
import ClimaUtilities.TimeManager: ITime
import Oceananigans, ClimaSeaIce, KernelAbstractions, ConservativeRegridding, Adapt
import ClimaLand, ClimaAtmos
import ClimaComms
ClimaComms.@import_required_backends
import YAML
using CouplerWarmupCommon
import AtmosWarmupSetup  # layered: its cached code is reused, not re-serialised
import AtmosWarmupExplicit  # layered: its cached code is reused, not re-serialised
import AtmosWarmupStep  # layered: its cached code is reused, not re-serialised
import OceanWarmup  # layered: its cached code is reused, not re-serialised
import LandWarmup  # layered: its cached code is reused, not re-serialised
import SeaIceWarmup  # layered: its cached code is reused, not re-serialised
using PrecompileTools

const WARMUP_DIR = normpath(joinpath(@__DIR__, "..", ".."))
const REPO_DIR = normpath(joinpath(WARMUP_DIR, "..", ".."))
const CONFIG = joinpath(WARMUP_DIR, "configs", "cmip_warmup.yml")
include_dependency(CONFIG)
include_dependency(joinpath(WARMUP_DIR, "configs", "atmos_warmup.yml"))

const REPORT = joinpath(@__DIR__, "..", "workload_report.txt")

"Load the warmup config with repository-relative paths resolved and output redirected."
"""
    small_gpu_atmos_config(path, outdir)

Path to the atmosphere config, with `h_elem` reduced when `WARMUP_H_ELEM` is set.

The coupled workloads hold the atmosphere, ocean and sea ice on the device at once, which
at production resolution exceeds a 15 GB card. Lowering `h_elem` cuts that cost without
costing cache validity, because horizontal resolution is *not* part of the type: a field's
data layout is `ClimaCore.DataLayouts.VIJFH{S, Nv, Nij, Nij, ...}`, which carries the
vertical level count and the quadrature size but not the number of elements. Measured
against an image built at `h_elem = 16`, `ClimaAtmos.get_simulation` at `h_elem = 8`
compiled 0.2 s -- the same as at 16.

`z_elem` is the opposite and must never be lowered for this purpose: it is the `Nv` type
parameter above, and against the same image `get_simulation` compiled 0.2 s at 79 but
250.6 s at 78 and 243.8 s at 63. One level of difference costs the whole cache.

The value is recorded in `workload_report.txt`.
"""
function small_gpu_atmos_config(path, outdir)
    h_elem = get(ENV, "WARMUP_H_ELEM", "")
    isempty(h_elem) && return path
    cfg = YAML.load_file(path)
    cfg["h_elem"] = parse(Int, h_elem)
    out = joinpath(outdir, "atmos_warmup_h$(h_elem).yml")
    YAML.write_file(out, cfg)
    return out
end

function warmup_config(outdir)
    cfg = YAML.load_file(CONFIG)
    cfg["coupler_output_dir"] = outdir
    cfg["job_id"] = "warmup"
    for key in ("atmos_config_file",)
        v = get(cfg, key, nothing)
        v isa AbstractString && !isabspath(v) && (cfg[key] = joinpath(REPO_DIR, v))
    end
    cfg["atmos_config_file"] = small_gpu_atmos_config(cfg["atmos_config_file"], outdir)
    tmp = joinpath(outdir, "cmip_warmup_resolved.yml")
    YAML.write_file(tmp, cfg)
    return Input.get_coupler_config_dict(tmp)
end

@setup_workload begin
    @compile_workload begin
        status = "not attempted"
        try
            outdir = mktempdir(; cleanup = true)
            config_dict = warmup_config(outdir)
            cs = ClimaCoupler.CoupledSimulation(config_dict)
            # `run!` is the entry point a production run uses, so caching it also covers
            # its own code: the timing path, `wait_slow_sims!` teardown, SYPD recording
            # and closing the output writers -- none of which a bare `step!` reaches.
            #
            # It evolves the whole span, which matters: a path on a longer cadence than
            # the workload is never compiled. `t_end` is set so that land diagnostics
            # (30mins) and gravity-wave drag (dt_ogw/dt_nogw, 1800secs) both fire; the
            # config shortens the remaining cadences so they fire within the first step.
            ClimaCoupler.run!(cs)
            status = "ok"
        catch err
            status = "FAILED: " * first(sprint(showerror, err), 300)
        end
        try
            open(REPORT, "w") do io
                println(io, "status = ", status)
                ho = get(ENV, "WARMUP_H_ELEM", "")
                isempty(ho) || println(io, "h_elem override = ", ho, " (type-neutral; see README)")
                println(io, "device = ", get(ENV, "CLIMACOMMS_DEVICE", "(unset)"))
            end
        catch
        end
    end
end

end # module
