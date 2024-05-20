# # AMIP Driver

#=
## Overview

AMIP is a standard experimental protocol of the Program for Climate Model Diagnosis & Intercomparison (PCMDI).
It is used as a model benchmark for the atmospheric and land model components, while sea-surface temperatures (SST) and sea-ice concentration (SIC)
are prescribed using time-interpolations between monthly observed data. We use standard data files with original sources:
- SST and SIC: https://gdex.ucar.edu/dataset/158_asphilli.html
- land-sea mask: https://www.ncl.ucar.edu/Applications/Data/#cdf

For more information, see the PCMDI's specifications for [AMIP I](https://pcmdi.github.io/mips/amip/) and [AMIP II](https://pcmdi.github.io/mips/amip2/).

This driver contains two modes. The full `AMIP` mode and a `SlabPlanet` (all surfaces are thermal slabs) mode. Since `AMIP` is not a closed system, the
`SlabPlanet` mode is useful for checking conservation properties of the coupling.

=#

#=
## Logging
When Julia 1.10+ is used interactively, stacktraces contain reduced type information to make them shorter.
Given that ClimaCore objects are heavily parametrized, non-abbreviated stacktraces are hard to read,
so we force abbreviated stacktraces even in non-interactive runs.
(See also `Base.type_limited_string_from_context()`)
=#

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
import Dates
import YAML

# ## ClimaESM packages
import ClimaAtmos as CA
import ClimaComms
import ClimaCore as CC

# ## Coupler specific imports
import ClimaCoupler
import ClimaCoupler:
    BCReader,
    ConservationChecker,
    Checkpointer,
    Diagnostics,
    FieldExchanger,
    FluxCalculator,
    Interfacer,
    Regridder,
    TimeManager,
    Utilities

pkg_dir = pkgdir(ClimaCoupler)

#=
### Helper Functions
These will be eventually moved to their respective component model and diagnostics packages, and so they should not
contain any internals of the ClimaCoupler source code, except extensions to the Interfacer functions.
=#

## helpers for component models
include("components/atmosphere/climaatmos.jl")
include("components/land/climaland_bucket.jl")
include("components/ocean/slab_ocean.jl")
include("components/ocean/prescr_seaice.jl")
include("components/ocean/eisenman_seaice.jl")

## helpers for user-specified IO
include("user_io/user_diagnostics.jl")
include("user_io/user_logging.jl")
include("user_io/debug_plots.jl")
include("user_io/io_helpers.jl")

#=
### Configuration Dictionaries
Each simulation mode has its own configuration dictionary. The `config_dict` of each simulation is a merge of the default configuration
dictionary and the simulation-specific configuration dictionary, which allows the user to override the default settings.

We can additionally pass the configuration dictionary to the component model initializers, which will then override the default settings of the component models.
=#

## coupler simulation default configuration
include("cli_options.jl")
parsed_args = parse_commandline(argparse_settings())

## modify parsed args for fast testing from REPL #hide
if isinteractive()
    parsed_args["config_file"] =
        isnothing(parsed_args["config_file"]) ? joinpath(pkg_dir, "config/ci_configs/interactive_debug.yml") :
        parsed_args["config_file"]
end

## read in config dictionary from file, overriding the coupler defaults in `parsed_args`
config_dict = YAML.load_file(parsed_args["config_file"])
config_dict = merge(parsed_args, config_dict)

## get component model dictionaries (if applicable)
atmos_config_dict, config_dict = get_atmos_config_dict(config_dict)
atmos_config_object = CA.AtmosConfig(atmos_config_dict)

## read in some parsed command line arguments, required by this script
mode_name = config_dict["mode_name"]
run_name = config_dict["run_name"]
energy_check = config_dict["energy_check"]
const FT = config_dict["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
land_sim_name = "bucket"
t_end = Float64(time_to_seconds(config_dict["t_end"]))
t_start = 0.0
tspan = (t_start, t_end)
Δt_cpl = Float64(config_dict["dt_cpl"])
saveat = Float64(time_to_seconds(config_dict["dt_save_to_sol"]))
date0 = date = Dates.DateTime(config_dict["start_date"], Dates.dateformat"yyyymmdd")
mono_surface = config_dict["mono_surface"]
hourly_checkpoint = config_dict["hourly_checkpoint"]
hourly_checkpoint_dt = config_dict["hourly_checkpoint_dt"]
restart_dir = config_dict["restart_dir"]
restart_t = Int(config_dict["restart_t"])
evolving_ocean = config_dict["evolving_ocean"]
dt_rad = config_dict["dt_rad"]
use_coupler_diagnostics = config_dict["use_coupler_diagnostics"]

#=
## Setup Communication Context
We set up communication context for CPU single thread/CPU with MPI/GPU. If no device is passed to `ClimaComms.context()`
then `ClimaComms` automatically selects the device from which this code is called.
=#

comms_ctx = Utilities.get_comms_context(parsed_args)
ClimaComms.init(comms_ctx)

#=
### I/O Directory Setup
`setup_output_dirs` returns `dir_paths.output = COUPLER_OUTPUT_DIR`, which is the directory where the output of the simulation will be saved, and `dir_paths.artifacts` is the directory where
the plots (from postprocessing and the conservation checks) of the simulation will be saved. `dir_paths.regrid` is the directory where the regridding
temporary files will be saved.
=#

COUPLER_OUTPUT_DIR = joinpath(config_dict["coupler_output_dir"], joinpath(mode_name, run_name))
dir_paths = setup_output_dirs(output_dir = COUPLER_OUTPUT_DIR, comms_ctx = comms_ctx)


if ClimaComms.iamroot(comms_ctx)
    @info(dir_paths.output)
    config_dict["print_config_dict"] && @info(config_dict)
end

#=
## Data File Paths
The data files are downloaded from the `ClimaCoupler` artifacts directory. If the data files are not present, they are downloaded from the
original sources.
=#

include(joinpath(pkgdir(ClimaCoupler), "artifacts", "artifact_funcs.jl"))
sst_data = artifact_data(sst_dataset_path(), "sst", "SST", dir_paths.regrid, date0, t_start, t_end, comms_ctx)
sic_data = artifact_data(sic_dataset_path(), "sic", "SEAICE", dir_paths.regrid, date0, t_start, t_end, comms_ctx)
co2_data = artifact_data(co2_dataset_path(), "mauna_loa_co2", "co2", dir_paths.regrid, date0, t_start, t_end, comms_ctx)
land_mask_data = artifact_data(mask_dataset_path(), "seamask")

#=
## Component Model Initialization
Here we set initial and boundary conditions for each component model. Each component model is required to have an `init` function that
returns a `ComponentModelSimulation` object (see `Interfacer` docs for more details).
=#

#=
### Atmosphere
This uses the `ClimaAtmos.jl` model, with parameterization options specified in the `atmos_config_object` dictionary.
=#

Utilities.show_memory_usage(comms_ctx)

## init atmos model component
atmos_sim = atmos_init(atmos_config_object);
Utilities.show_memory_usage(comms_ctx)

thermo_params = get_thermo_params(atmos_sim) # TODO: this should be shared by all models #610

#=
### Boundary Space
We use a common `Space` for all global surfaces. This enables the MPI processes to operate on the same columns in both
the atmospheric and surface components, so exchanges are parallelized. Note this is only possible when the
atmosphere and surface are of the same horizontal resolution.
=#

## init a 2D boundary space at the surface
boundary_space = CC.Spaces.horizontal_space(atmos_sim.domain.face_space) # TODO: specify this in the coupler and pass it to all component models #665

#=
### Land-sea Fraction
This is a static field that contains the area fraction of land and sea, ranging from 0 to 1. If applicable, sea ice is included in the sea fraction. at this stage.
Note that land-sea area fraction is different to the land-sea mask, which is a binary field (masks are used internally by the coupler to indicate passive cells that are not populated by a given component model).
=#

land_area_fraction =
    FT.(
        Regridder.land_fraction(
            FT,
            dir_paths.regrid,
            comms_ctx,
            land_mask_data,
            "LSMASK",
            boundary_space,
            mono = mono_surface,
        )
    )
Utilities.show_memory_usage(comms_ctx)

#=
### Surface Models: AMIP and SlabPlanet Modes
Both modes evolve `ClimaLand.jl`'s bucket model.

In the `AMIP` mode, all ocean properties are prescribed from a file, while sea-ice temperatures are calculated using observed
SIC and assuming a 2m thickness of the ice.

In the `SlabPlanet` mode, all ocean and sea ice are dynamical models, namely thermal slabs, with different parameters. We have several `SlabPlanet` versions
- `slabplanet` = land + slab ocean
- `slabplanet_aqua` = slab ocean everywhere
- `slabplanet_terra` = land everywhere
- `slabplanet_eisenman` = land + slab ocean + slab sea ice with an evolving thickness
=#

ClimaComms.iamroot(comms_ctx) && @info(mode_name)
if mode_name == "amip"
    ClimaComms.iamroot(comms_ctx) && @info("AMIP boundary conditions - do not expect energy conservation")

    ## land model
    land_sim = bucket_init(
        FT,
        tspan,
        config_dict["land_domain_type"],
        config_dict["land_albedo_type"],
        config_dict["land_temperature_anomaly"],
        dir_paths.regrid;
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_area_fraction,
        date_ref = date0,
        t_start = t_start,
        energy_check = energy_check,
    )

    ## ocean stub
    SST_info = BCReader.bcfile_info_init(
        FT,
        dir_paths.regrid,
        sst_data,
        "SST",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        scaling_function = scale_sst, ## convert to Kelvin
        land_fraction = land_area_fraction,
        date0 = date0,
        mono = mono_surface,
    )

    BCReader.update_midmonth_data!(date0, SST_info)
    SST_init = BCReader.interpolate_midmonth_to_daily(date0, SST_info)
    ocean_sim = Interfacer.SurfaceStub((;
        T_sfc = SST_init,
        ρ_sfc = CC.Fields.zeros(boundary_space),
        z0m = FT(1e-3),
        z0b = FT(1e-3),
        beta = FT(1),
        α_direct = CC.Fields.ones(boundary_space) .* FT(0.06),
        α_diffuse = CC.Fields.ones(boundary_space) .* FT(0.06),
        area_fraction = (FT(1) .- land_area_fraction),
        phase = TD.Liquid(),
        thermo_params = thermo_params,
    ))

    ## sea ice model
    SIC_info = BCReader.bcfile_info_init(
        FT,
        dir_paths.regrid,
        sic_data,
        "SEAICE",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        scaling_function = scale_sic, ## convert to fraction
        land_fraction = land_area_fraction,
        date0 = date0,
        mono = mono_surface,
    )
    BCReader.update_midmonth_data!(date0, SIC_info)
    SIC_init = BCReader.interpolate_midmonth_to_daily(date0, SIC_info)
    ice_fraction = get_ice_fraction.(SIC_init, mono_surface)
    ice_sim = ice_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = ice_fraction,
        thermo_params = thermo_params,
    )

    ## CO2 concentration from temporally varying file
    CO2_info = BCReader.bcfile_info_init(
        FT,
        dir_paths.regrid,
        co2_data,
        "co2",
        boundary_space,
        comms_ctx,
        interpolate_daily = true,
        land_fraction = ones(boundary_space),
        date0 = date0,
        mono = mono_surface,
    )

    BCReader.update_midmonth_data!(date0, CO2_info)
    CO2_init = BCReader.interpolate_midmonth_to_daily(date0, CO2_info)
    Interfacer.update_field!(atmos_sim, Val(:co2), CO2_init)

    mode_specifics = (; name = mode_name, SST_info = SST_info, SIC_info = SIC_info, CO2_info = CO2_info)
    Utilities.show_memory_usage(comms_ctx)

elseif mode_name in ("slabplanet", "slabplanet_aqua", "slabplanet_terra")


    land_area_fraction = mode_name == "slabplanet_aqua" ? land_area_fraction .* 0 : land_area_fraction
    land_area_fraction = mode_name == "slabplanet_terra" ? land_area_fraction .* 0 .+ 1 : land_area_fraction

    ## land model
    land_sim = bucket_init(
        FT,
        tspan,
        config_dict["land_domain_type"],
        config_dict["land_albedo_type"],
        config_dict["land_temperature_anomaly"],
        dir_paths.regrid;
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_area_fraction,
        date_ref = date0,
        t_start = t_start,
        energy_check = energy_check,
    )

    ## ocean model
    ocean_sim = ocean_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = (FT(1) .- land_area_fraction), ## NB: this ocean fraction includes areas covered by sea ice (unlike the one contained in the cs)
        thermo_params = thermo_params,
        evolving = evolving_ocean,
    )

    ## sea ice stub (here set to zero area coverage)
    ice_sim = Interfacer.SurfaceStub((;
        T_sfc = CC.Fields.ones(boundary_space),
        ρ_sfc = CC.Fields.zeros(boundary_space),
        z0m = FT(0),
        z0b = FT(0),
        beta = FT(1),
        α_direct = CC.Fields.ones(boundary_space) .* FT(1),
        α_diffuse = CC.Fields.ones(boundary_space) .* FT(1),
        area_fraction = CC.Fields.zeros(boundary_space),
        phase = TD.Ice(),
        thermo_params = thermo_params,
    ))

    mode_specifics = (; name = mode_name, SST_info = nothing, SIC_info = nothing)
    Utilities.show_memory_usage(comms_ctx)

elseif mode_name == "slabplanet_eisenman"

    ## land model
    land_sim = bucket_init(
        FT,
        tspan,
        config_dict["land_domain_type"],
        config_dict["land_albedo_type"],
        config_dict["land_temperature_anomaly"],
        dir_paths.regrid;
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = land_area_fraction,
        date_ref = date0,
        t_start = t_start,
        energy_check = energy_check,
    )

    ## ocean stub (here set to zero area coverage)
    ocean_sim = ocean_init(
        FT;
        tspan = tspan,
        dt = Δt_cpl,
        space = boundary_space,
        saveat = saveat,
        area_fraction = CC.Fields.zeros(boundary_space), # zero, since ML is calculated below
        thermo_params = thermo_params,
    )

    ## sea ice + ocean model
    ice_sim = eisenman_seaice_init(
        FT,
        tspan,
        space = boundary_space,
        area_fraction = (FT(1) .- land_area_fraction),
        dt = Δt_cpl,
        saveat = saveat,
        thermo_params = thermo_params,
    )

    mode_specifics = (; name = mode_name, SST_info = nothing, SIC_info = nothing)
    Utilities.show_memory_usage(comms_ctx)
end

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
    NamedTuple{coupler_field_names}(ntuple(i -> CC.Fields.zeros(boundary_space), length(coupler_field_names)))
Utilities.show_memory_usage(comms_ctx)

## model simulations
model_sims = (atmos_sim = atmos_sim, ice_sim = ice_sim, land_sim = land_sim, ocean_sim = ocean_sim);

## dates
dates = (; date = [date], date0 = [date0], date1 = [Dates.firstdayofmonth(date0)], new_month = [false])

#=
### Online Diagnostics
User can write custom diagnostics in the `user_diagnostics.jl`.
Note, this will be replaced by the diagnostics framework currently in ClimaAtmos, once it is abstracted
into a more general package, so we can use it to save fields from surface models.
=#
if use_coupler_diagnostics
    monthly_3d_diags = Diagnostics.init_diagnostics(
        (:T, :u, :q_tot, :q_liq_ice),
        atmos_sim.domain.center_space;
        save = TimeManager.Monthly(),
        operations = (; accumulate = Diagnostics.TimeMean([Int(0)])),
        output_dir = COUPLER_OUTPUT_DIR,
        name_tag = "monthly_mean_3d_",
    )

    monthly_3d_diags = Diagnostics.init_diagnostics(
        (:T, :u, :q_tot, :q_liq_ice),
        atmos_sim.domain.center_space;
        save = TimeManager.Monthly(),
        operations = (; accumulate = Diagnostics.TimeMean([Int(0)])),
        output_dir = dir_paths.output,
        name_tag = "monthly_mean_3d_",
    )

    monthly_2d_diags = Diagnostics.init_diagnostics(
        (:precipitation_rate, :toa_fluxes, :T_sfc, :tubulent_energy_fluxes),
        boundary_space;
        save = TimeManager.Monthly(),
        operations = (; accumulate = Diagnostics.TimeMean([Int(0)])),
        output_dir = dir_paths.output,
        name_tag = "monthly_mean_2d_",
    )

    diagnostics = (monthly_3d_diags, monthly_2d_diags)
    Utilities.show_memory_usage(comms_ctx)
else
    diagnostics = ()
end

#=
## Initialize Conservation Checks

The conservation checks are used to monitor the global energy and water conservation of the coupled system. The checks are only
applicable to the `slabplanet` mode, as the `amip` mode is not a closed system. The conservation checks are initialized here and
saved in a global `ConservationChecks` struct, `conservation_checks`.
=#

## init conservation info collector
conservation_checks = nothing
if energy_check
    @assert(
        mode_name[1:10] == "slabplanet" && !CA.is_distributed(ClimaComms.context(boundary_space)),
        "Only non-distributed slabplanet allowable for energy_check"
    )
    conservation_checks = (;
        energy = ConservationChecker.EnergyConservationCheck(model_sims),
        water = ConservationChecker.WaterConservationCheck(model_sims),
    )
end

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

checkpoint_cb = TimeManager.HourlyCallback(
    dt = hourly_checkpoint_dt,
    func = checkpoint_sims,
    ref_date = [dates.date[1]],
    active = hourly_checkpoint,
) # 20 days
update_firstdayofmonth!_cb = TimeManager.MonthlyCallback(
    dt = FT(1),
    func = TimeManager.update_firstdayofmonth!,
    ref_date = [dates.date1[1]],
    active = true,
)
dt_water_albedo = parse(FT, filter(x -> !occursin(x, "hours"), dt_rad))
albedo_cb = TimeManager.HourlyCallback(
    dt = dt_water_albedo,
    func = FluxCalculator.water_albedo_from_atmosphere!,
    ref_date = [dates.date[1]],
    active = mode_name == "amip",
)
callbacks =
    (; checkpoint = checkpoint_cb, update_firstdayofmonth! = update_firstdayofmonth!_cb, water_albedo = albedo_cb)

#=
## Initialize turbulent fluxes

Decide on the type of turbulent flux partition (see `FluxCalculator` documentation for more details).
=#
turbulent_fluxes = nothing
if config_dict["turb_flux_partition"] == "PartitionedStateFluxes"
    turbulent_fluxes = FluxCalculator.PartitionedStateFluxes()
elseif config_dict["turb_flux_partition"] == "CombinedStateFluxesMOST"
    turbulent_fluxes = FluxCalculator.CombinedStateFluxesMOST()
else
    error("turb_flux_partition must be either PartitionedStateFluxes or CombinedStateFluxesMOST")
end

#=
## Initialize Coupled Simulation

The coupled simulation is initialized here and saved in a global `CoupledSimulation` struct, `cs`. It contains all the information
required to run the coupled simulation, including the communication context, the dates, the boundary space, the coupler fields, the
configuration dictionary, the conservation checks, the time span, the time step, the land fraction, the model simulations, the mode
specifics, the diagnostics, the callbacks, and the directory paths.
=#

cs = Interfacer.CoupledSimulation{FT}(
    comms_ctx,
    dates,
    boundary_space,
    coupler_fields,
    config_dict,
    conservation_checks,
    [tspan[1], tspan[2]],
    atmos_sim.integrator.t,
    Δt_cpl,
    (; land = land_area_fraction, ocean = zeros(boundary_space), ice = zeros(boundary_space)),
    model_sims,
    mode_specifics,
    diagnostics,
    callbacks,
    dir_paths,
    turbulent_fluxes,
    thermo_params,
);
Utilities.show_memory_usage(comms_ctx)

#=
## Restart component model states if specified
If a restart directory is specified and contains output files from the `checkpoint_cb` callback, the component model states are restarted from those files. The restart directory
is specified in the `config_dict` dictionary. The `restart_t` field specifies the time step at which the restart is performed.
=#

if restart_dir !== "unspecified"
    for sim in cs.model_sims
        if Checkpointer.get_model_prog_state(sim) !== nothing
            Checkpointer.restart_model_state!(sim, comms_ctx, restart_t; input_dir = restart_dir)
        end
    end
end

#=
## Initialize Component Model Exchange

We need to ensure all models' initial conditions are shared to enable the coupler to calculate the first instance of surface fluxes. Some auxiliary variables (namely surface humidity and radiation fluxes)
depend on initial conditions of other component models than those in which the variables are calculated, which is why we need to step these models in time and/or reinitialize them.
The concrete steps for proper initialization are:
=#

# 1.coupler updates surface model area fractions
Regridder.update_surface_fractions!(cs)

# 2.surface density (`ρ_sfc`): calculated by the coupler by adiabatically extrapolating atmospheric thermal state to the surface.
# For this, we need to import surface and atmospheric fields. The model sims are then updated with the new surface density.
FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes)
FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes)
FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

# 3.surface vapor specific humidity (`q_sfc`): step surface models with the new surface density to calculate their respective `q_sfc` internally
## TODO: the q_sfc calculation follows the design of the bucket q_sfc, but it would be neater to abstract this from step! (#331)
Interfacer.step!(land_sim, Δt_cpl)
Interfacer.step!(ocean_sim, Δt_cpl)
Interfacer.step!(ice_sim, Δt_cpl)

# 4.turbulent fluxes: now we have all information needed for calculating the initial turbulent surface fluxes using the combined state
# or the partitioned state method
if cs.turbulent_fluxes isa FluxCalculator.CombinedStateFluxesMOST
    ## import the new surface properties into the coupler (note the atmos state was also imported in step 3.)
    FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes) # i.e. T_sfc, albedo, z0, beta, q_sfc
    ## calculate turbulent fluxes inside the atmos cache based on the combined surface state in each grid box
    FluxCalculator.combined_turbulent_fluxes!(cs.model_sims, cs.fields, cs.turbulent_fluxes) # this updates the atmos thermo state, sfc_ts
elseif cs.turbulent_fluxes isa FluxCalculator.PartitionedStateFluxes
    ## calculate turbulent fluxes in surface models and save the weighted average in coupler fields
    FluxCalculator.partitioned_turbulent_fluxes!(
        cs.model_sims,
        cs.fields,
        cs.boundary_space,
        FluxCalculator.MoninObukhovScheme(),
        cs.thermo_params,
    )

    ## update atmos sfc_conditions for surface temperature
    ## TODO: this is hard coded and needs to be simplified (req. CA modification) (#479)
    new_p = get_new_cache(atmos_sim, cs.fields)
    CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t) ## sets T_sfc (but SF calculation not necessary - requires split functionality in CA)
    atmos_sim.integrator.p.precomputed.sfc_conditions .= new_p.precomputed.sfc_conditions
end

# 5.reinitialize models + radiative flux: prognostic states and time are set to their initial conditions. For atmos, this also triggers the callbacks and sets a nonzero radiation flux (given the new sfc_conditions)
FieldExchanger.reinit_model_sims!(cs.model_sims)

# 6.update all fluxes: coupler re-imports updated atmos fluxes (radiative fluxes for both `turbulent_fluxes` types
# and also turbulent fluxes if `turbulent_fluxes isa CombinedStateFluxesMOST`,
# and sends them to the surface component models. If `turbulent_fluxes isa PartitionedStateFluxes`
# atmos receives the turbulent fluxes from the coupler.
FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes)
FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

#=
## Coupling Loop

The coupling loop is the main part of the simulation. It runs the component models sequentially for one coupling timestep (`Δt_cpl`), and exchanges combined fields and calculates fluxes using combined states.
Note that we want to implement this in a dispatchable function to allow for other forms of timestepping (e.g. leapfrog). (TODO: #610)
=#

function solve_coupler!(cs)
    (; model_sims, Δt_cpl, tspan, comms_ctx) = cs
    (; atmos_sim, land_sim, ocean_sim, ice_sim) = model_sims

    ClimaComms.iamroot(comms_ctx) && @info("Starting coupling loop")
    ## step in time
    for t in ((tspan[begin] + Δt_cpl):Δt_cpl:tspan[end])

        cs.dates.date[1] = TimeManager.current_date(cs, t)

        ## print date on the first of month
        if cs.dates.date[1] >= cs.dates.date1[1]
            ClimaComms.iamroot(comms_ctx) && @show(cs.dates.date[1])
        end

        if cs.mode.name == "amip"

            ## monthly read of boundary condition data for SST and SIC and CO2
            if cs.dates.date[1] >= BCReader.next_date_in_file(cs.mode.SST_info)
                BCReader.update_midmonth_data!(cs.dates.date[1], cs.mode.SST_info)
            end
            SST_current = BCReader.interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SST_info)
            Interfacer.update_field!(ocean_sim, Val(:surface_temperature), SST_current)

            if cs.dates.date[1] >= BCReader.next_date_in_file(cs.mode.SIC_info)
                BCReader.update_midmonth_data!(cs.dates.date[1], cs.mode.SIC_info)
            end
            SIC_current =
                get_ice_fraction.(
                    BCReader.interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.SIC_info),
                    cs.mode.SIC_info.mono,
                )
            Interfacer.update_field!(ice_sim, Val(:area_fraction), SIC_current)

            if cs.dates.date[1] >= BCReader.next_date_in_file(cs.mode.CO2_info)
                BCReader.update_midmonth_data!(cs.dates.date[1], cs.mode.CO2_info)
            end
            CO2_current = BCReader.interpolate_midmonth_to_daily(cs.dates.date[1], cs.mode.CO2_info)
            Interfacer.update_field!(atmos_sim, Val(:co2), CO2_current)

            ## calculate and accumulate diagnostics at each timestep, if we're using diagnostics in this run
            if !isempty(cs.diagnostics)
                ClimaComms.barrier(comms_ctx)
                Diagnostics.accumulate_diagnostics!(cs)

                ## save and reset monthly averages
                Diagnostics.save_diagnostics(cs)
            end
        end

        ## compute global energy
        !isnothing(cs.conservation_checks) && ConservationChecker.check_conservation!(cs)
        ClimaComms.barrier(comms_ctx)

        ## update water albedo from wind at dt_water_albedo (this will be extended to a radiation callback from the coupler)
        TimeManager.trigger_callback!(cs, cs.callbacks.water_albedo)

        ## run component models sequentially for one coupling timestep (Δt_cpl)
        Regridder.update_surface_fractions!(cs)
        FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

        ## step sims
        FieldExchanger.step_model_sims!(cs.model_sims, t)

        ## exchange combined fields and (if specified) calculate fluxes using combined states
        FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes) # i.e. T_sfc, surface_albedo, z0, beta
        if cs.turbulent_fluxes isa FluxCalculator.CombinedStateFluxesMOST
            FluxCalculator.combined_turbulent_fluxes!(cs.model_sims, cs.fields, cs.turbulent_fluxes) # this updates the surface thermo state, sfc_ts, in ClimaAtmos (but also unnecessarily calculates fluxes)
        elseif cs.turbulent_fluxes isa FluxCalculator.PartitionedStateFluxes
            ## calculate turbulent fluxes in surfaces and save the weighted average in coupler fields
            FluxCalculator.partitioned_turbulent_fluxes!(
                cs.model_sims,
                cs.fields,
                cs.boundary_space,
                FluxCalculator.MoninObukhovScheme(),
                cs.thermo_params,
            )

            ## update atmos sfc_conditions for surface temperature - TODO: this needs to be simplified (need CA modification)
            new_p = get_new_cache(atmos_sim, cs.fields)
            CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t) # to set T_sfc (but SF calculation not necessary - CA modification)
            atmos_sim.integrator.p.precomputed.sfc_conditions .= new_p.precomputed.sfc_conditions
        end

        FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes) # radiative and/or turbulent

        ## callback to update the fist day of month if needed (for BCReader)
        TimeManager.trigger_callback!(cs, cs.callbacks.update_firstdayofmonth!)

        ## callback to checkpoint model state
        TimeManager.trigger_callback!(cs, cs.callbacks.checkpoint)

    end

    return nothing
end

## exit if running performance anaysis #hide
if haskey(ENV, "CI_PERF_SKIP_COUPLED_RUN") #hide
    throw(:exit_profile_init) #hide
end #hide

## run the coupled simulation for one timestep to precompile everything before timing
cs.tspan[2] = Δt_cpl * 2
solve_coupler!(cs)

## run the coupled simulation for the full timespan and time it
cs.tspan[1] = Δt_cpl * 2
cs.tspan[2] = tspan[2]

## Run garbage collection before solving for more accurate memory comparison to ClimaAtmos
GC.gc()

## Use ClimaComms.@elapsed to time the simulation on both CPU and GPU
walltime = ClimaComms.@elapsed comms_ctx.device begin
    s = CA.@timed_str begin
        solve_coupler!(cs)
    end
end
ClimaComms.iamroot(comms_ctx) && @show(walltime)

## Use ClimaAtmos calculation to show the simulated years per day of the simulation (SYPD)
es = CA.EfficiencyStats(tspan, walltime)
sypd = CA.simulated_years_per_day(es)
@info "SYPD: $sypd"

## Save the SYPD and allocation information
if ClimaComms.iamroot(comms_ctx)
    sypd_filename = joinpath(dir_paths.artifacts, "sypd.txt")
    write(sypd_filename, "$sypd")

    cpu_allocs_GB = Utilities.show_memory_usage(comms_ctx)
    cpu_allocs_filename = joinpath(dir_paths.artifacts, "allocations_cpu.txt")
    write(cpu_allocs_filename, cpu_allocs_GB)
end

#=
## Postprocessing
Currently all postprocessing is performed using the root process only.
=#

if ClimaComms.iamroot(comms_ctx)

    ## energy check plots
    if !isnothing(cs.conservation_checks) && cs.mode.name[1:10] == "slabplanet"
        @info "Conservation Check Plots"
        ConservationChecker.plot_global_conservation(
            cs.conservation_checks.energy,
            cs,
            config_dict["conservation_softfail"],
            figname1 = joinpath(dir_paths.artifacts, "total_energy_bucket.png"),
            figname2 = joinpath(dir_paths.artifacts, "total_energy_log_bucket.png"),
        )
        ConservationChecker.plot_global_conservation(
            cs.conservation_checks.water,
            cs,
            config_dict["conservation_softfail"],
            figname1 = joinpath(dir_paths.artifacts, "total_water_bucket.png"),
            figname2 = joinpath(dir_paths.artifacts, "total_water_log_bucket.png"),
        )
    end

    ## sample animations (not compatible with MPI)
    if !CA.is_distributed(comms_ctx) && config_dict["anim"]
        @info "Animations"
        include("user_io/viz_explorer.jl")
        plot_anim(cs, dir_paths.artifacts)
    end

    ## plotting AMIP results
    if cs.mode.name == "amip" && !isempty(cs.diagnostics)
        ## plot data that correspond to the model's last save_hdf5 call (i.e., last month)
        @info "AMIP plots"

        ## ClimaESM
        include("user_io/amip_visualizer.jl")
        post_spec = (;
            T = (:regrid, :zonal_mean),
            u = (:regrid, :zonal_mean),
            q_tot = (:regrid, :zonal_mean),
            toa_fluxes = (:regrid, :horizontal_slice),
            precipitation_rate = (:regrid, :horizontal_slice),
            T_sfc = (:regrid, :horizontal_slice),
            tubulent_energy_fluxes = (:regrid, :horizontal_slice),
            q_liq_ice = (:regrid, :zonal_mean),
        )

        plot_spec = (;
            T = (; clims = (190, 320), units = "K"),
            u = (; clims = (-50, 50), units = "m/s"),
            q_tot = (; clims = (0, 30), units = "g/kg"),
            toa_fluxes = (; clims = (-250, 250), units = "W/m^2"),
            precipitation_rate = (clims = (0, 1e-4), units = "kg/m^2/s"),
            T_sfc = (clims = (225, 310), units = "K"),
            tubulent_energy_fluxes = (; clims = (-250, 250), units = "W/m^2"),
            q_liq_ice = (; clims = (0, 10), units = "g/kg"),
        )
        amip_data, fig_amip = amip_paperplots(
            post_spec,
            plot_spec,
            dir_paths.output,
            files_root = ".monthly",
            output_dir = dir_paths.artifacts,
        )

        ## NCEP reanalysis
        @info "NCEP plots"
        include("user_io/ncep_visualizer.jl")
        ncep_post_spec = (;
            T = (:zonal_mean,),
            u = (:zonal_mean,),
            q_tot = (:zonal_mean,),
            toa_fluxes = (:horizontal_slice,),
            precipitation_rate = (:horizontal_slice,),
            T_sfc = (:horizontal_slice,),
            tubulent_energy_fluxes = (:horizontal_slice,),
        )
        ncep_plot_spec = plot_spec
        ncep_data, fig_ncep = ncep_paperplots(
            ncep_post_spec,
            ncep_plot_spec,
            dir_paths.output,
            output_dir = dir_paths.artifacts,
            month_date = cs.dates.date[1],
        )

        ## combine AMIP and NCEP plots
        plot_combined = Plots.plot(fig_amip, fig_ncep, layout = (2, 1), size = (1400, 1800))
        Plots.png(joinpath(dir_paths.artifacts, "amip_ncep.png"))

        ## Compare against observations
        if t_end > 84600 && config_dict["output_default_diagnostics"]
            @info "Error against observations"
            output_dates = cs.dates.date0[] .+ Dates.Second.(atmos_sim.integrator.sol.t)

            include("user_io/leaderboard.jl")
            compare_vars = ["pr"]
            function compute_biases(dates)
                if isempty(dates)
                    return map(x -> 0.0, compare_vars)
                else
                    return Leaderboard.compute_biases(atmos_sim.integrator.p.output_dir, compare_vars, dates)
                end
            end

            function plot_biases(dates, biases, output_name)
                isempty(dates) && return nothing

                output_path = joinpath(dir_paths.artifacts, "bias_$(output_name).png")
                Leaderboard.plot_biases(biases; output_path)
            end

            ann_biases = compute_biases(output_dates)
            plot_biases(output_dates, ann_biases, "total")

            ## collect all days between cs.dates.date0 and cs.dates.date
            MAM, JJA, SON, DJF = Leaderboard.split_by_season(output_dates)

            MAM_biases = compute_biases(MAM)
            plot_biases(MAM, MAM_biases, "MAM")
            JJA_biases = compute_biases(JJA)
            plot_biases(JJA, JJA_biases, "JJA")
            SON_biases = compute_biases(SON)
            plot_biases(SON, SON_biases, "SON")
            DJF_biases = compute_biases(DJF)
            plot_biases(DJF, DJF_biases, "DJF")

            rmses = map(
                (index) -> Leaderboard.RMSEs(;
                    model_name = "CliMA",
                    ANN = ann_biases[index],
                    DJF = DJF_biases[index],
                    JJA = JJA_biases[index],
                    MAM = MAM_biases[index],
                    SON = SON_biases[index],
                ),
                1:length(compare_vars),
            )

            Leaderboard.plot_leaderboard(rmses; output_path = "bias_leaderboard.png")
        end
    end

    ## ci plots
    if config_dict["ci_plots"]
        @info "Generating CI plots"
        include("user_io/ci_plots.jl")
        make_plots(Val(:general_ci_plots), [atmos_sim.integrator.p.output_dir], dir_paths.artifacts)
    end

    ## plot all model states and coupler fields (useful for debugging)
    !(comms_ctx isa ClimaComms.MPICommsContext) && debug(cs, dir_paths.artifacts)

    if isinteractive()
        ## clean up for interactive runs, retain all output otherwise
        rm(dir_paths.output; recursive = true, force = true)
    end

end
