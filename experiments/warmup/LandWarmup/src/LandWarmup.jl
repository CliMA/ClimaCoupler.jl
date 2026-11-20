"""
    LandWarmup

Caches the ClimaLand component as a coupled run constructs it.

Part of the layered precompilation workloads in `experiments/warmup`. The dependency edges
mirror the call graph and are required, not an optimisation: without them a package
re-caches the code it calls and its image exceeds the ~2 GiB limit imposed by 32-bit ELF
relocations. See `experiments/warmup/README.md`.

`CLIMACOMMS_DEVICE=CUDA` must be set when precompiling.
"""
module LandWarmup

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
            P = prologue(config_dict)
            a = P.args
            land_sim = Interfacer.LandSimulation(
                a.FT,
                a.land_model;
                dt = a.component_dt_dict["dt_land"],
                tspan = (a.t_start, a.t_end),
                a.start_date,
                output_dir = P.dir_paths.land_output_dir,
                area_fraction = P.land_fraction,
                surface_space = P.boundary_space,
                P.atmos_h,
                P.initial_T,
                a.use_land_diagnostics,
                a.land_diagnostics_period,
                a.land_diagnostics_reduction,
                P.coupled_param_dict,
                albedo_type = a.bucket_albedo_type,
                bucket_initial_condition = P.bucket_initial_condition,
                era5_albedo_file_path = P.albedo_path,
                a.land_spun_up_ic,
                land_ic_path = P.land_ic_path,
                a.lai_source,
                dt_drivers = ITime(Utilities.time_to_seconds(config_dict["dt_rad"])),
            )
            Interfacer.step!(land_sim, a.component_dt_dict["dt_land"])
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
