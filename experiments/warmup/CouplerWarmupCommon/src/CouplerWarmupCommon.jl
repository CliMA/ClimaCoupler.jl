"""
    CouplerWarmupCommon

Shared setup for the component workloads (`OceanWarmup`, `LandWarmup`, `SeaIceWarmup`).
"""
module CouplerWarmupCommon

import ClimaCoupler
import ClimaCoupler: Interfacer, Input, Utilities, SimCoordinator
import ClimaParams as CP
import Thermodynamics.Parameters as TDP
import ClimaComms

# Loaded because `experiments/CMIP/code_loading.jl` loads it; see `AtmosWarmupSetup`.
import CairoMakie, GeoMakie, Makie, Poppler_jll, Printf

export prologue

"""
    prologue(config_dict; with_atmos = true)

Return, as a `NamedTuple`, the inputs a coupled run builds before its component models,
mirroring the start of `SimCoordinator._build_coupled_simulation`. With `with_atmos = false`
the atmosphere and everything derived from its boundary space are skipped.
"""
function prologue(config_dict::AbstractDict; with_atmos::Bool = true)
    comms_ctx = Utilities.get_comms_context(config_dict)
    args = Input.get_coupler_args(config_dict)
    (; FT, domain_type, column_latlon, parameter_files) = args

    override_file = CP.merge_toml_files(parameter_files; override = true)
    coupled_param_dict = CP.create_toml_dict(FT; override_file)
    thermo_params = TDP.ThermodynamicsParameters(coupled_param_dict)
    dir_paths = Utilities.setup_output_dirs(
        output_dir_root = args.output_dir_root,
        comms_ctx = comms_ctx,
    )

    with_atmos || return (;
        comms_ctx,
        args,
        coupled_param_dict,
        thermo_params,
        dir_paths,
        atmos_sim = nothing,
        boundary_space = nothing,
    )

    atmos_sim = Interfacer.AtmosSimulation(
        Val(:climaatmos);
        config_dict,
        atmos_output_dir = dir_paths.atmos_output_dir,
        coupled_param_dict,
        comms_ctx,
    )
    (; boundary_space, atmos_h, initial_T) =
        SimCoordinator.boundary_fields(FT, domain_type, atmos_sim, comms_ctx, column_latlon)

    # `era5_filepaths` is not part of `get_coupler_args`: `_build_coupled_simulation`
    # derives it separately and only in subseasonal mode. Mirrored here.
    sst_path = sic_path = land_ic_path = albedo_path = nothing
    bucket_initial_condition = args.bucket_initial_condition
    if args.sim_mode <: Interfacer.SubseasonalMode
        (; sst_path, sic_path, land_ic_path, albedo_path, bucket_initial_condition) =
            Input.get_era5_filepaths(
                config_dict["era5_initial_condition_dir"],
                args.start_date,
                bucket_initial_condition,
                comms_ctx,
            )
    end

    land_fraction = Input.get_land_fraction(
        boundary_space,
        comms_ctx;
        args.land_fraction_source,
        args.binary_area_fraction,
        args.sim_mode,
        domain_type,
        args.scm_surface_type,
    )

    return (;
        comms_ctx,
        args,
        coupled_param_dict,
        thermo_params,
        dir_paths,
        atmos_sim,
        boundary_space,
        atmos_h,
        initial_T,
        land_fraction,
        sst_path,
        sic_path,
        land_ic_path,
        albedo_path,
        bucket_initial_condition,
    )
end


end # module
