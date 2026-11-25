"""
    LandWarmup

Precompilation workload caching ClimaLand as a coupled run builds and steps it.
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
# Lower layers of the warmup stack, reused rather than re-cached.
import AtmosWarmupSetup
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
