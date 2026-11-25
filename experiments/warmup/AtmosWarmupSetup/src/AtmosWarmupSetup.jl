"""
    AtmosWarmupSetup

Precompilation workload caching atmosphere setup, `ClimaAtmos.get_simulation`.
"""
module AtmosWarmupSetup

import ClimaAtmos as CA
import ClimaComms
ClimaComms.@import_required_backends
import ClimaTimeSteppers
import YAML
using PrecompileTools

# Everything `experiments/CMIP/code_loading.jl` loads, though unused here: a package loaded
# later would invalidate the cached code. See `precompile_info.md`.
import CairoMakie, GeoMakie, Makie, Poppler_jll, Printf
import ClimaCoupler
import Oceananigans, ClimaSeaIce, KernelAbstractions, ConservativeRegridding, Adapt
import ClimaLand

include(joinpath(@__DIR__, "..", "..", "warmup_utils.jl"))

# Rebuild when the config changes, not only the source.
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
            status = "FAILED: " * first(sprint(showerror, err), 300)
        end
        write_report(joinpath(@__DIR__, "..", "workload_report.txt"), status)
    end
end

end # module
