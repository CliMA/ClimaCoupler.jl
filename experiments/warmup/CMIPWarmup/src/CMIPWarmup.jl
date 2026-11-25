"""
    CMIPWarmup

Precompilation workload caching a coupled `ClimaCoupler.run!`: field exchange, fluxes and
the exchange grid. Loading it loads every other warmup image.
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
# Lower layers of the warmup stack, reused rather than re-cached.
import AtmosWarmupSetup
import AtmosWarmupExplicit
import AtmosWarmupStep
import OceanWarmup
import LandWarmup
import SeaIceWarmup
using PrecompileTools

include(joinpath(@__DIR__, "..", "..", "warmup_utils.jl"))

# Rebuild when the config changes, not only the source.
include_dependency(config_path("cmip_warmup.yml"))
include_dependency(config_path("atmos_warmup.yml"))

@setup_workload begin
    @compile_workload begin
        status = "not attempted"
        try
            outdir = mktempdir(; cleanup = true)
            config_dict = warmup_config(outdir)
            cs = ClimaCoupler.CoupledSimulation(config_dict)
            # `run!`, not `step!`, to also cache teardown and output closing. `t_end` covers
            # the longest schedule (30 min); the config shortens the rest.
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
