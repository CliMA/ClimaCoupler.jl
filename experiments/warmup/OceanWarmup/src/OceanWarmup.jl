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

"Load the warmup config with repository-relative paths resolved and output redirected."
function warmup_config(outdir)
    cfg = YAML.load_file(CONFIG)
    cfg["coupler_output_dir"] = outdir
    cfg["job_id"] = "warmup"
    for key in ("atmos_config_file",)
        v = get(cfg, key, nothing)
        v isa AbstractString && !isabspath(v) && (cfg[key] = joinpath(REPO_DIR, v))
    end
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
                println(io, "device = ", get(ENV, "CLIMACOMMS_DEVICE", "(unset)"))
            end
        catch
        end
    end
end

end # module
