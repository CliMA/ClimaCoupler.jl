"""
    CMIPWarmup

Precompilation workload that caches a coupled `ClimaCoupler.run!`: field exchange, flux
calculation and the exchange grid.

Sits on top of the other warmup packages, so loading it loads all their images.

Precompile with `CLIMACOMMS_DEVICE=CUDA`, or the image holds CPU code. See
`experiments/warmup/README.md`.
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
            # `run!` rather than `step!`, so the timing, teardown and output-closing code
            # that production reaches is cached too. Code behind a schedule that never
            # fires is not compiled, so `t_end` covers the longest cadence (land
            # diagnostics and gravity-wave drag, 30 min) and the config shortens the rest.
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
