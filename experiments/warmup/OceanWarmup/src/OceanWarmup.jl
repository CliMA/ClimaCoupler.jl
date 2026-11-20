"""
    OceanWarmup

Precompilation workload caching the Oceananigans ocean as a coupled run builds and steps it.
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
# Layered: reuse the code cached by the warmup packages beneath this one.
import AtmosWarmupSetup
using PrecompileTools

include(joinpath(@__DIR__, "..", "..", "warmup_utils.jl"))

# Rebuild the image when a config it reads changes, not only when the source does.
include_dependency(config_path("cmip_warmup.yml"))
include_dependency(config_path("atmos_warmup.yml"))

@setup_workload begin
    @compile_workload begin
        status = "not attempted"
        try
            outdir = mktempdir(; cleanup = true)
            config_dict = warmup_config(outdir)
            P = prologue(config_dict)
            a = P.args
            ocean_sim = Interfacer.OceanSimulation(
                a.FT,
                a.ocean_model;
                dt = a.component_dt_dict["dt_ocean"],
                a.start_date,
                tspan = (a.t_start, a.t_end),
                P.coupled_param_dict,
                P.thermo_params,
                P.comms_ctx,
                P.boundary_space,
                output_dir = P.dir_paths.ocean_output_dir,
                a.simple_ocean,
                a.ocean_grid,
                a.use_intersection_grid,
                P.sst_path,
                a.sst_adjustment,
                a.saveat,
                evolving = a.evolving_ocean,
                a.ocean_diagnostic_interval,
                a.ocean_diagnostic_mode,
            )
            Interfacer.step!(ocean_sim, a.component_dt_dict["dt_ocean"])
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
