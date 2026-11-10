"""
    CouplerWarmupCommon

Shared setup for the per-component precompilation workloads.

Each component's warmup package needs the same prerequisites a coupled run builds before it
can construct its model: a comms context, the coupled parameter dictionaries, an atmosphere
simulation (for the boundary space), and the land fraction. `prologue` assembles exactly
those, mirroring the opening of `ClimaCoupler.SimCoordinator._build_coupled_simulation`.

This duplicates a small amount of that constructor. The duplication is deliberate and
temporary: factoring the component construction out of `CoupledSimulation` into reusable
functions belongs in its own change. Until then, `check_matches_coupled` guards against
drift by failing loudly if what we build stops matching what a coupled run builds.
"""
module CouplerWarmupCommon

import ClimaCoupler
import ClimaCoupler: Interfacer, Input, Utilities, SimCoordinator
import ClimaParams as CP
import Thermodynamics.Parameters as TDP
import ClimaComms

# Production loads the Makie stack before it builds anything -- `experiments/CMIP/code_loading.jl`
# does so to trigger `ClimaCouplerMakieExt`. Loading it *after* this workload's code has been
# compiled invalidates that code, so the cache has to be built in the same world the run uses.
# Measured: with the stack absent from the workload, `ClimaAtmos.get_simulation` compiled 0.2 s
# in a session without Makie and 227.1 s in one with it, and coupled setup went 4.4 s -> 1153.1 s.
# These imports are load-bearing even though nothing here calls them.
import CairoMakie, GeoMakie, Makie, Poppler_jll, Printf


export prologue, check_matches_coupled

"""
    prologue(config_dict; with_atmos = true)

Build the prerequisites shared by every component model and return them as a NamedTuple,
together with the parsed coupler arguments in `args`.

`with_atmos = false` skips the atmosphere simulation (and therefore the boundary space);
only useful for components that need neither.
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

Fail loudly when a warmup builds a different type than a coupled run would.

Cached specializations are keyed on concrete types, so a mismatch silently yields an image
the production run cannot use — which looks exactly like "precompilation didn't help".
This is the check that caught the atmosphere workload defaulting to `Float32` and
`Val{:NoWarp}` because `FLOAT_TYPE` and `topography` are set in the *coupler* config.
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
