"""
    OceanWarmup

Caches the Oceananigans ocean component as a coupled run constructs it.

Part of the layered precompilation workloads in `experiments/warmup`. The dependency edges
mirror the call graph and are required, not an optimisation: without them a package
re-caches the code it calls and its image exceeds the ~2 GiB limit imposed by 32-bit ELF
relocations. See `experiments/warmup/README.md`.

`CLIMACOMMS_DEVICE=CUDA` must be set when precompiling.
"""
module OceanWarmup

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
using PrecompileTools

const WARMUP_DIR = normpath(joinpath(@__DIR__, "..", ".."))
const REPO_DIR = normpath(joinpath(WARMUP_DIR, "..", ".."))
const CONFIG = joinpath(WARMUP_DIR, "configs", "cmip_warmup.yml")
include_dependency(CONFIG)
include_dependency(joinpath(WARMUP_DIR, "configs", "atmos_warmup.yml"))

const REPORT = joinpath(@__DIR__, "..", "workload_report.txt")

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

"Load the warmup config with repository-relative paths resolved and output redirected."
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
            P = prologue(config_dict)
            a = P.args
            ocean_sim = Interfacer.OceanSimulation(a.FT, a.ocean_model;
                dt = a.component_dt_dict["dt_ocean"], a.start_date, tspan = (a.t_start, a.t_end),
                P.coupled_param_dict, P.thermo_params, P.comms_ctx, P.boundary_space,
                output_dir = P.dir_paths.ocean_output_dir, a.simple_ocean, a.ocean_grid,
                a.use_intersection_grid, P.sst_path, a.sst_adjustment, a.saveat,
                evolving = a.evolving_ocean, a.ocean_diagnostic_interval, a.ocean_diagnostic_mode)
            Interfacer.step!(ocean_sim, a.component_dt_dict["dt_ocean"])
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
