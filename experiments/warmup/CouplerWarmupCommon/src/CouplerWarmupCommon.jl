"""
    CouplerWarmupCommon

Shared setup for the component precompilation workloads (`OceanWarmup`, `LandWarmup`,
`SeaIceWarmup`).

`prologue` builds what a coupled run builds before its component models, mirroring the
start of `ClimaCoupler.SimCoordinator._build_coupled_simulation`. `check_matches_coupled`
checks that a workload built the same type as a coupled run.
"""
module CouplerWarmupCommon

import ClimaCoupler
import ClimaCoupler: Interfacer, Input, Utilities, SimCoordinator
import ClimaParams as CP
import Thermodynamics.Parameters as TDP
import ClimaComms

# Loaded, although nothing here calls it, because `experiments/CMIP/code_loading.jl` loads
# it: a package loaded after this workload is compiled invalidates the cached code. See
# `experiments/warmup/precompile_info.md`.
import CairoMakie, GeoMakie, Makie, Poppler_jll, Printf

export prologue, check_matches_coupled

"""
    prologue(config_dict; with_atmos = true)

Build the inputs shared by every component model, as a coupled run builds them.

# Arguments
- `config_dict`: coupler config, as returned by `Input.get_coupler_config_dict`.

# Keyword Arguments
- `with_atmos = true`: build the atmosphere simulation and everything derived from its
  boundary space. With `false`, `atmos_sim` and `boundary_space` are `nothing`.

# Returns
A `NamedTuple` with `comms_ctx`, `args` (from `Input.get_coupler_args`),
`coupled_param_dict`, `thermo_params`, `dir_paths`, `atmos_sim` and `boundary_space`. With
`with_atmos = true` it also holds `atmos_h`, `initial_T`, `land_fraction`, `sst_path`,
`sic_path`, `land_ic_path`, `albedo_path` and `bucket_initial_condition`.
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

"""
    check_matches_coupled(label, built, expected_type)

Return `true` if `typeof(built) === expected_type`, and throw an error naming `label`
otherwise.

Cached specializations are keyed on concrete types, so a workload that builds a different
type than a coupled run produces an image the run cannot use.
"""
function check_matches_coupled(label, built, expected_type)
    actual = typeof(built)
    actual === expected_type && return true
    error("""
          $label built $(actual)
          but a coupled run builds $(expected_type).
          The warmup config has drifted from the production config; the cached
          specializations would not be usable.
          """)
end

end # module
