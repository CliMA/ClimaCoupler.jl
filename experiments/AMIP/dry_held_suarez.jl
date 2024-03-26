# # Dry Held-Suarez

redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))

#=
## Configuration initialization
Here we import standard Julia packages, ClimaESM packages, parse in command-line arguments (if none are specified then the defaults in `cli_options.jl` apply).
We then specify the input data file names. If these are not already downloaded, include `artifacts/download_artifacts.jl`.
=#

#=
### Package Import
=#

## standard packages
using Dates
import YAML

# ## ClimaESM packages
import ClimaAtmos as CA
using ClimaCore

# ## Coupler specific imports
using ClimaCoupler
using ClimaCoupler.BCReader: bcfile_info_init, update_midmonth_data!, next_date_in_file, interpolate_midmonth_to_daily
using ClimaCoupler.ConservationChecker:
    EnergyConservationCheck, WaterConservationCheck, check_conservation!, plot_global_conservation
using ClimaCoupler.Checkpointer: restart_model_state!
using ClimaCoupler.Diagnostics: init_diagnostics, accumulate_diagnostics!, save_diagnostics, TimeMean
using ClimaCoupler.FieldExchanger:
    import_atmos_fields!, import_combined_surface_fields!, update_model_sims!, reinit_model_sims!, step_model_sims!
using ClimaCoupler.FluxCalculator:
    PartitionedStateFluxes,
    CombinedStateFluxes,
    combined_turbulent_fluxes!,
    MoninObukhovScheme,
    partitioned_turbulent_fluxes!,
    water_albedo_from_atmosphere!
using ClimaCoupler.Interfacer: CoupledSimulation, SurfaceStub, get_field, update_field!
using ClimaCoupler.Regridder
using ClimaCoupler.Regridder: update_surface_fractions!, combine_surfaces!, binary_mask
using ClimaCoupler.TimeManager:
    current_date, Monthly, EveryTimestep, HourlyCallback, MonthlyCallback, update_firstdayofmonth!, trigger_callback!
import ClimaCoupler.Utilities: get_comms_context

pkg_dir = pkgdir(ClimaCoupler)

#=
### Helper Functions
These will be eventually moved to their respective component model and diagnostics packages, and so they should not
contain any internals of the ClimaCoupler source code, except extensions to the Interfacer functions.
=#

## helpers for component models
include("components/atmosphere/climaatmos.jl")

## helpers for user-specified IO
include("user_io/user_diagnostics.jl")
include("user_io/user_logging.jl")

include("driver_utils.jl")
#=
### Configuration Dictionaries
Each simulation mode has its own configuration dictionary. The `config_dict` of each simulation is a merge of the default configuration
dictionary and the simulation-specific configuration dictionary, which allows the user to override the default settings.

We can additionally pass the configuration dictionary to the component model initializers, which will then override the default settings of the component models.
=#

## run names
run_name = "dry_held_suarez"
coupler_output_dir = "$run_name"
const FT = Float64
restart_dir = "unspecified"
restart_t = Int(0)

## coupler simulation specific configuration
Δt_cpl = Float64(450)
t_end = "200days"
tspan = (Float64(0.0), Float64(time_to_seconds(t_end)))
start_date = "19790301"
hourly_checkpoint = true

## atmos arguments to override
config_dict = Dict(
    # file paths
    "atmos_config_file" => nothing,
    "coupler_toml_file" => nothing,
    "coupler_output_dir" => coupler_output_dir,
    "mode_name" => "",
    "run_name" => run_name,
    # timestepping
    "dt" => "$(Δt_cpl)secs",
    "dt_save_to_sol" => "10days",
    "t_end" => t_end,
    "start_date" => "19790301",
    # domain
    # "h_elem" => 6,#16,
    "z_elem" => 25,#63,
    "z_max" => 45000.0, # semi-high top
    "dz_bottom" => 300.0,
    # "dz_top" => 3000.0,
    # output
    "dt_save_state_to_disk" => "1000days",
    # numerics
    "apply_limiter" => false,
    "rayleigh_sponge" => true,
    # run
    "job_id" => run_name,
    "surface_setup" => "PrescribedSurface",
    # diagnostic (nested wirh period and short_name)
    "diagnostics" => [Dict("short_name" => ["mse", "lr", "mass_streamfunction", "stab", "vT", "egr"], "period" => "1days"),],
    # held-suarez specific
    # "forcing" => "held_suarez",
    "viscous_sponge" => true,
    "rayleigh_sponge" => true,
    "vert_diff" => "true",
)


# Alternative
# dt_save_state_to_disk: "10days"
# dt: "150secs"
# t_end: "300days"
# h_elem: 16
# z_elem: 63
# dz_bottom: 30.0
# dz_top: 3000.0
# z_max: 55000.0
# vert_diff: "true"
# moist: "equil"
# precip_model: "0M"
# rayleigh_sponge: true
# forcing: "held_suarez"
# job_id: "longrun_hs_rhoe_equil_55km_nz63_0M"
# toml: [toml/longrun_hs_rhoe_equil_55km_nz63_0M.toml]

# cirrently
#precip_model: "0M"
#moist: "equil"
# toml: [toml/sphere_held_suarez_rhoe_equilmoist_hightop_sponge.toml]


## merge dictionaries of command line arguments, coupler dictionary and component model dictionaries
## (if there are common keys, the last dictorionary in the `merge` arguments takes precedence)
config_dict_atmos = get_atmos_config_dict(config_dict)
config_dict = merge(config_dict_atmos, config_dict)
atmos_config_object = CA.AtmosConfig(config_dict_atmos)
# By passing `parsed_args` to `AtmosConfig`, `parsed_args` overwrites the default atmos config

# overriding parameter values
atmos_config_object.toml_dict["zd_viscous"]["value"] = 30000.0
atmos_config_object.toml_dict["zd_rayleigh"]["value"] = 30000.0

#=
## Setup Communication Context
We set up communication context for CPU single thread/CPU with MPI/GPU. If no device is passed to `ClimaComms.context()`
then `ClimaComms` automatically selects the device from which this code is called.
=#

using ClimaComms
comms_ctx = get_comms_context(Dict("device" => "auto"))
ClimaComms.init(comms_ctx)


#=
### I/O Directory Setup
`coupler_output_dir` is the directory where the output of the simulation will be saved, and `COUPLER_ARTIFACTS_DIR` is the directory where
the plots (from postprocessing and the conservation checks) of the simulation will be saved. `REGRID_DIR` is the directory where the regridding
temporary files will be saved.
=#

dir_paths = setup_output_dirs(output_dir = coupler_output_dir, comms_ctx = comms_ctx)

ClimaComms.iamroot(comms_ctx) ? @info(config_dict) : nothing

#=
## Data File Paths
The data files are downloaded from the `ClimaCoupler` artifacts directory. If the data files are not present, they are downloaded from the
original sources.
=#

include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))
co2_data = joinpath(co2_dataset_path(), "mauna_loa_co2.nc")

#=
## Component Model Initialization
Here we set initial and boundary conditions for each component model. Each component model is required to have an `init` function that
returns a `ComponentModelSimulation` object (see `Interfacer` docs for more details).

### Atmosphere
This uses the `ClimaAtmos.jl` model, with parameterization options specified in the `config_dict_atmos` dictionary.
=#

## init atmos model component
atmos_sim = atmos_init(FT, atmos_config_object);
thermo_params = get_thermo_params(atmos_sim)

#=
### Boundary Space
We use a common `Space` for all global surfaces. This enables the MPI processes to operate on the same columns in both
the atmospheric and surface components, so exchanges are parallelized. Note this is only possible when the
atmosphere and surface are of the same horizontal resolution.
=#

## init a 2D boundary space at the surface
boundary_space = ClimaCore.Spaces.horizontal_space(atmos_sim.domain.face_space) # TODO: specify this in the coupler and pass it to all component models #665

#=
## Coupler Initialization
The coupler needs to contain exchange information, manage the calendar and be able to access all component models. It can also optionally
save online diagnostics. These are all initialized here and saved in a global `CoupledSimulation` struct, `cs`.
=#

## coupler exchange fields
coupler_field_names = (
    :T_S,
    :z0m_S,
    :z0b_S,
    :ρ_sfc,
    :q_sfc,
    :surface_direct_albedo,
    :surface_diffuse_albedo,
    :beta,
    :F_turb_energy,
    :F_turb_moisture,
    :F_turb_ρτxz,
    :F_turb_ρτyz,
    :F_radiative,
    :P_liq,
    :P_snow,
    :radiative_energy_flux_toa,
    :P_net,
    :temp1,
    :temp2,
)
coupler_fields =
    NamedTuple{coupler_field_names}(ntuple(i -> ClimaCore.Fields.zeros(boundary_space), length(coupler_field_names)))

## model simulations
model_sims = (atmos_sim = atmos_sim,);

## dates
date0 = date = DateTime(start_date, dateformat"yyyymmdd")
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)], new_month = [false])

#=
### Online Diagnostics
User can write custom diagnostics in the `user_diagnostics.jl`.
Note, this will be replaced by the diagnostics framework currently in ClimaAtmos, once it is abstracted
into a more general package, so we can use it to save fields from surface models.
=#

# monthly_3d_diags = init_diagnostics(
#     (:T, :u, :q_tot, :q_liq_ice),
#     atmos_sim.domain.center_space;
#     save = Monthly(),
#     operations = (; accumulate = TimeMean([Int(0)])),
#     output_dir = dir_paths.output,
#     name_tag = "monthly_mean_3d_",
# )

# monthly_2d_diags = init_diagnostics(
#     (:precipitation_rate,),
#     boundary_space;
#     save = Monthly(),
#     operations = (; accumulate = TimeMean([Int(0)])),
#     output_dir = dir_paths.output,
#     name_tag = "monthly_mean_2d_",
# )

diagnostics = ()

#=
## Initialize Callbacks
Callbacks are used to update at a specified interval. The callbacks are initialized here and
saved in a global `Callbacks` struct, `callbacks`. The `trigger_callback!` function is used to call the callback during the simulation below.

The frequency of the callbacks is specified in the `HourlyCallback` and `MonthlyCallback` structs. The `func` field specifies the function to be called,
the `ref_date` field specifies the reference (first) date for the callback, and the `active` field specifies whether the callback is active or not.

The currently implemented callbacks are:
- `checkpoint_cb`: generates a checkpoint of all model states at a specified interval. This is mainly used for restarting simulations.
- `update_firstdayofmonth!_cb`: generates a callback to update the first day of the month for monthly message print (and other monthly operations).
- `albedo_cb`: for the amip mode, the water albedo is time varying (since the reflectivity of water depends on insolation and wave characteristics, with the latter
  being approximated from wind speed). It is updated at the same frequency as the atmospheric radiation.
  NB: Eventually, we will call all of radiation from the coupler, in addition to the albedo calculation.
=#
checkpoint_cb =
    HourlyCallback(dt = FT(480), func = checkpoint_sims, ref_date = [dates.date[1]], active = hourly_checkpoint) # 20 days TODO: not GPU friendly
update_firstdayofmonth!_cb =
    MonthlyCallback(dt = FT(1), func = update_firstdayofmonth!, ref_date = [dates.date1[1]], active = true)
callbacks =
    (; checkpoint = checkpoint_cb, update_firstdayofmonth! = update_firstdayofmonth!_cb)

cs = CoupledSimulation{FT}(
    comms_ctx,
    dates,
    boundary_space,
    coupler_fields,
    config_dict,
    nothing, # conservation checks
    [tspan[1], tspan[2]],
    atmos_sim.integrator.t,
    Δt_cpl,
    (; land = land_fraction, ocean = zeros(boundary_space), ice = zeros(boundary_space)),
    model_sims,
    (;), # mode_specifics
    diagnostics,
    callbacks,
    dir_paths,
    nothing, # turbulent_fluxes
    thermo_params,
);

#=
## Restart component model states if specified
If a restart directory is specified and contains output files from the `checkpoint_cb` callback, the component model states are restarted from those files. The restart directory
is specified in the `config_dict` dictionary. The `restart_t` field specifies the time step at which the restart is performed.
=#

if restart_dir !== "unspecified"
    for sim in cs.model_sims
        if get_model_prog_state(sim) !== nothing
            restart_model_state!(sim, comms_ctx, restart_t; input_dir = restart_dir)
        end
    end
end

#=
## Coupling Loop

The coupling loop is the main part of the simulation. It runs the component models sequentially for one coupling timestep (`Δt_cpl`), and exchanges combined fields and calculates fluxes using combined states.
Note that we want to implement this in a dispatchable function to allow for other forms of timestepping (e.g. leapfrog). (TODO: #610)
=#

function solve_coupler!(cs)
    (; model_sims, Δt_cpl, tspan, comms_ctx) = cs
    (; atmos_sim) = model_sims

    ClimaComms.iamroot(comms_ctx) ? @info("Starting coupling loop") : nothing
    ## step in time
    walltime = @elapsed for t in ((tspan[begin] + Δt_cpl):Δt_cpl:tspan[end])

        cs.dates.date[1] = current_date(cs, t)

        ## print date on the first of month
        if cs.dates.date[1] >= cs.dates.date1[1]
            ClimaComms.iamroot(comms_ctx) ? @show(cs.dates.date[1]) : nothing
        end

        ## step sims
        step_model_sims!(cs.model_sims, t)

        import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes) # radiative and/or turbulent

        ## callback to update the fist day of month if needed (for BCReader)
        trigger_callback!(cs, cs.callbacks.update_firstdayofmonth!)

        ## callback to checkpoint model state
        trigger_callback!(cs, cs.callbacks.checkpoint)

    end
    ClimaComms.iamroot(comms_ctx) ? @show(walltime) : nothing

    return cs
end

## run the coupled simulation
solve_coupler!(cs);

# Postprocessing
anim = true
energy_check = false
# include("postprocessing.jl")


# climate diagnostics: T, u, (q), rho v, N, MSE

# storm track diagnostics: [vT], [v][T], EGR = f/N dTdy



