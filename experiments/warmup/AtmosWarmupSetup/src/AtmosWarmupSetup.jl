"""
    AtmosWarmupSetup

Precompilation workload that caches the atmosphere setup path, `ClimaAtmos.get_simulation`.

Precompile with `CLIMACOMMS_DEVICE=CUDA`, or the image holds CPU code. See
`experiments/warmup/README.md`.
"""
module AtmosWarmupSetup

import ClimaAtmos as CA
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers
import YAML
using PrecompileTools

# Everything `experiments/CMIP/code_loading.jl` loads, although nothing here calls it: a
# package loaded after this workload is compiled invalidates the cached code and its GPU
# kernels. See `experiments/warmup/precompile_info.md`.
import CairoMakie, GeoMakie, Makie, Poppler_jll, Printf
import ClimaCoupler
import Oceananigans, ClimaSeaIce, KernelAbstractions, ConservativeRegridding, Adapt
import ClimaLand


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
            @assert sim.integrator !== nothing
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
