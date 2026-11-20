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

include(joinpath(@__DIR__, "..", "..", "warmup_utils.jl"))

# Julia invalidates a package image when its *source* changes, not when data it read during
# precompilation changes. Register the config so editing it forces a rebuild.
include_dependency(config_path("cmip_warmup.yml"))
include_dependency(config_path("atmos_warmup.yml"))

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
        write_report(
            joinpath(@__DIR__, "..", "workload_report.txt"),
            status;
            h_elem_override = get(ENV, "WARMUP_H_ELEM", ""),
        )
    end
end

end # module
