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
import DelimitedFiles

# ## ClimaESM packages
import ClimaAtmos as CA
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC
import ClimaParams as CP
import Thermodynamics.Parameters as TDP

# ## Coupler specific imports
using ClimaCoupler
import ClimaCoupler.Interfacer:
    AbstractSlabplanetSimulationMode,
    AMIPMode,
    CMIPMode,
    CoupledSimulation,
    SlabplanetAquaMode,
    SlabplanetMode,
    SlabplanetTerraMode,
    SubseasonalMode

import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.Utils: period_to_seconds_float
import ClimaUtilities.TimeManager: ITime, date
import Interpolations # triggers InterpolationsExt in ClimaUtilities
# Random is used by RRMTGP for some cloud properties
import Random

import ClimaDiagnostics as CD
import ClimaDiagnostics.Schedules: EveryCalendarDtSchedule, EveryStepSchedule

# Trigger ClimaCouplerMakieExt extension
using Makie, GeoMakie, CairoMakie, ClimaCoreMakie, NCDatasets, Poppler_jll

# Trigger ClimaCouplerCMIPExt extension
# Note we only need these if running CMIP, but for now we share one environment for all experiments
import Oceananigans, ClimaOcean, ClimaSeaIce, KernelAbstractions

# Trigger ClimaCouplerClimaLandExt extension
import ClimaLand

# Trigger ClimaCouplerClimaAtmosExt
import ClimaAtmos

#=
### Configuration Dictionaries
Each simulation mode has its own configuration dictionary. The `config_dict` of each simulation is a merge of the default configuration
dictionary and the simulation-specific configuration dictionary, which allows the user to override the default settings.

We can additionally pass the configuration dictionary to the component model initializers, which will then override the default settings of the component models.
=#


"""
    CoupledSimulation(config_file)
    CoupledSimulation(config_dict)

Set up a `CoupledSimulation` as prescribed by the given input.

This struct is defined in the Interfacer module and contains all information
about component models, diagnostics, timestepping, output directories, etc
needed to run a coupled simulation.
"""
function CoupledSimulation(
    config_file = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml"),
)
    config_dict = Input.get_coupler_config_dict(config_file)
    return CoupledSimulation(config_dict)
end

function CoupledSimulation(config_dict::AbstractDict)
    # Initialize communication context (do this first so all printing is only on root)
    comms_ctx = Utilities.get_comms_context(config_dict)

    (;
        job_id,
        sim_mode,
        random_seed,
        FT,
        t_end,
        t_start,
        start_date,
        Δt_cpl,
        component_dt_dict,
        share_surface_space,
        nh_poly,
        h_elem,
        saveat,
        checkpoint_dt,
        detect_restart_files,
        restart_dir,
        restart_t,
        restart_cache,
        save_cache,
        use_land_diagnostics,
        diagnostics_dt,
        evolving_ocean,
        land_model,
        land_temperature_anomaly,
        land_spun_up_ic,
        lai_source,
        bucket_albedo_type,
        energy_check,
        use_coupler_diagnostics,
        output_dir_root,
        parameter_files,
        era5_filepaths,
        ocean_model,
        simple_ocean,
        sst_adjustment,
        ice_model,
        land_fraction_source,
        binary_area_fraction,
    ) = Input.get_coupler_args(config_dict)

    # Get default shared parameters from ClimaParams.jl, overriding with any provided parameter files
    override_file = CP.merge_toml_files(parameter_files; override = true)
    coupled_param_dict = CP.create_toml_dict(FT; override_file)
    thermo_params = TDP.ThermodynamicsParameters(coupled_param_dict)

    #=
    ### I/O Directory Setup
    `setup_output_dirs` returns a NamedTuple with the paths to
    output directories for each component and the coupler,
    as well as paths to artifacts, regrid, and checkpoints directories.
    =#
    dir_paths = Utilities.setup_output_dirs(
        output_dir_root = output_dir_root,
        comms_ctx = comms_ctx,
    )

    ## set unique random seed if desired, otherwise use default
    Random.seed!(random_seed)
    @info "Random seed set to $(random_seed)"

    tspan = (t_start, t_end)
    @info "Starting from t_start $(t_start)"

    #=
    ## Component Model Initialization
    Here we set initial and boundary conditions for each component model. Each component model is required to have an `init` function that
    returns a `AbstractComponentSimulation` object (see `Interfacer` docs for more details).
    =#

    #=
    ### Atmosphere
    This uses the `ClimaAtmos.jl` model, with parameterization options specified in the `atmos_config_object` dictionary.
    =#

    ## init atmos model component
    ## Note this step must come after parsing the coupler config dictionary, since
    ## some parameters are passed from the coupler config to the component model configs
    atmos_sim = Interfacer.AtmosSimulation(
        Val(:climaatmos);
        config_dict,
        atmos_output_dir = dir_paths.atmos_output_dir,
        coupled_param_dict,
        comms_ctx,
    )

    #=
    ### Boundary Space
    We use a 2D boundary space at the surface for coupling operations.
    The boundary space is used to remap exchanged fields and compute fluxes between component models.

    We currently have two options for the boundary space:
    - `share_surface_space: true`: If `true`, we use the horizontal space of the atmosphere model as the boundary space.
      This is useful when the atmosphere and surface models are of the same resolution;
      this case makes remapping trivial, as the coupler can directly exchange fields and fluxes.
    - `share_surface_space: false`: If `false`, we create a boundary space independent of the component models,
      using the `ClimaCommonSpaces.CubedSphereSpace` constructor. This is useful when the
      atmosphere and surface models are of different resolutions or grid types. In this case,
      we need to remap exchanged fields and fluxes between component and boundary spaces.
    =#

    ## init a 2D boundary space at the surface
    if share_surface_space
        boundary_space = CC.Spaces.horizontal_space(atmos_sim.domain.face_space)
    else
        n_quad_points = nh_poly + 1
        radius = coupled_param_dict["planet_radius"] # in meters
        boundary_space = CC.CommonSpaces.CubedSphereSpace(FT; radius, n_quad_points, h_elem)
    end

    # Get surface elevation on the boundary space from `atmos` coordinate field
    surface_elevation = Interfacer.get_field(boundary_space, atmos_sim, Val(:height_sfc)) # on boundary space
    # Get atmospheric height relative to the surface directly from the atmosphere
    atmos_h = Interfacer.get_atmos_height_delta(
        Interfacer.get_field(atmos_sim, Val(:height_int)),
        surface_elevation,
    )

    land_fraction = Input.get_land_fraction(
        boundary_space,
        comms_ctx;
        land_fraction_source,
        binary_area_fraction,
        sim_mode,
    )

    #=
    ### Surface Models: AMIP and SlabPlanet Modes
    Both modes evolve `ClimaLand.jl`'s bucket model.

    In the `AMIP` mode, all ocean properties are prescribed from a file, while sea-ice temperatures are calculated using observed
    SIC and assuming a 2m thickness of the ice.

    In the `SlabPlanet` mode, all ocean and sea ice are dynamical models, namely thermal slabs, with different parameters. We have several `SlabPlanet` versions
    - `slabplanet` = land + slab ocean
    - `slabplanet_aqua` = slab ocean everywhere
    - `slabplanet_terra` = land everywhere

    In this section of the code, we initialize all component models and read in the prescribed data we'll be using.
    The specific models and data that are set up depend on which mode we're running.
    =#

    @info(sim_mode)
    land_sim = ice_sim = ocean_sim = nothing

    # Unpack ERA5-based file paths (populated for subseasonal mode, nothing otherwise)
    (; sst_path, sic_path, land_ic_path, albedo_path, bucket_initial_condition) =
        era5_filepaths

    ## Construct the land model component
    # Determine whether to use a shared surface space
    shared_surface_space = share_surface_space ? boundary_space : nothing
    land_sim = Interfacer.LandSimulation(
        FT,
        land_model;
        # Arguments used by multiple models
        dt = component_dt_dict["dt_land"],
        tspan,
        start_date,
        output_dir = dir_paths.land_output_dir,
        area_fraction = land_fraction,
        shared_surface_space,
        saveat,
        surface_elevation,
        atmos_h,
        land_temperature_anomaly,
        use_land_diagnostics,
        coupled_param_dict,
        # Arguments used by bucket model
        albedo_type = bucket_albedo_type,
        bucket_initial_condition,
        era5_albedo_file_path = albedo_path,
        # Arguments used by integrated model
        land_spun_up_ic,
        land_ic_path,
        lai_source,
    )

    ## Construct the ocean model component
    ocean_sim = Interfacer.OceanSimulation(
        FT,
        ocean_model;
        # Arguments used by multiple models
        dt = component_dt_dict["dt_ocean"],
        start_date,
        tspan,
        coupled_param_dict,
        thermo_params,
        comms_ctx,
        boundary_space,
        # Arguments used by Oceananigans
        output_dir = dir_paths.ocean_output_dir,
        simple_ocean,
        # Arguments used by prescribed ocean
        sst_path,
        sst_adjustment,
        # Arguments used by slab ocean
        saveat,
        evolving = evolving_ocean,
    )

    ## Construct the sea ice model component (note this must be constructed after the ocean model)
    ice_sim = Interfacer.SeaIceSimulation(
        FT,
        ice_model;
        # Arguments used by multiple models
        dt = component_dt_dict["dt_seaice"],
        start_date,
        coupled_param_dict,
        output_dir = dir_paths.ice_output_dir,
        # Arguments used by ClimaSeaIce model
        ocean = ocean_sim,
        # Arguments used by prescribed ice model
        tspan,
        saveat,
        boundary_space,
        thermo_params,
        comms_ctx,
        land_fraction,
        sic_path,
        binary_area_fraction,
    )

    #=
    ## Coupler Initialization
    The coupler needs to contain exchange information, access all component models, and manage the calendar,
    among other responsibilities.
    Objects containing information to enable these are initialized here and saved in the
    global `CoupledSimulation` struct, `cs`, below.
    =#

    ## collect component model simulations that have been initialized
    model_sims = (; atmos_sim, ice_sim, land_sim, ocean_sim)
    model_sims =
        NamedTuple{filter(key -> !isnothing(model_sims[key]), keys(model_sims))}(model_sims)
    @info "Component models initialized: $(keys(model_sims))"

    ## coupler exchange fields
    coupler_field_names = Interfacer.default_coupler_fields()
    foreach(sim -> Interfacer.add_coupler_fields!(coupler_field_names, sim), model_sims)

    # add coupler fields required to track conservation, if specified
    energy_check && push!(coupler_field_names, :P_net)

    # allocate space for the coupler fields
    coupler_fields = Interfacer.init_coupler_fields(FT, coupler_field_names, boundary_space)

    #=
    ## Initialize Conservation Checks

    The conservation checks are used to monitor the global energy and water conservation of the coupled system. The checks are only
    applicable to the `slabplanet` mode, as the `amip` mode is not a closed system. The conservation checks are initialized here and
    saved in a global `ConservationChecks` struct, `conservation_checks`, which is then stored as part of the larger `cs` struct.
    =#

    ## init conservation info collector
    conservation_checks = nothing
    if energy_check
        @assert(
            sim_mode <: AbstractSlabplanetSimulationMode &&
            comms_ctx isa ClimaComms.SingletonCommsContext,
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
    saved in a global `Callbacks` struct, `callbacks`. The `callbacks!` function is used to call the callback during the simulation below.

    The currently implemented callbacks are:
    - `checkpoint_cb`: generates a checkpoint of all model states at a specified interval. This is mainly used for restarting simulations.
    =#
    schedule_checkpoint =
        EveryCalendarDtSchedule(TimeManager.time_to_period(checkpoint_dt); start_date)
    checkpoint_cb =
        TimeManager.Callback(schedule_checkpoint, sim -> Checkpointer.checkpoint_sims(sim))

    # Don't use coupler walltime logging if atmos is using its own walltime logging is true
    if config_dict["atmos_log_progress"]
        callbacks = (checkpoint_cb,)
    else
        walltime_cb = TimeManager.capped_geometric_walltime_cb(t_start, t_end, Δt_cpl)
        callbacks = (checkpoint_cb, walltime_cb)
    end

    #= Set up default AMIP diagnostics
    Use ClimaDiagnostics for default AMIP diagnostics, which currently include turbulent energy fluxes.
    =#
    if use_coupler_diagnostics
        @info "Using default coupler diagnostics"
        diags_handler = SimOutput.diagnostics_setup(
            coupler_fields,
            dir_paths.coupler_output_dir,
            start_date,
            tspan[1],
            diagnostics_dt,
            Δt_cpl,
        )
    else
        diags_handler = nothing
    end

    #=
    ## Initialize Coupled Simulation

    The coupled simulation is initialized here and saved in a global `CoupledSimulation` struct, `cs`. It contains all the information
    required to run the coupled simulation, including the communication context, the dates, the boundary space, the coupler fields, the
    configuration dictionary, the conservation checks, the time span, the time step, the land fraction, the model simulations, the mode
    specifics, the callbacks, the directory paths, and diagnostics for AMIP simulations.
    =#

    prev_checkpoint_t = Ref(-1) # no checkpoint taken yet
    cs = CoupledSimulation{FT}(
        Ref(start_date),
        coupler_fields,
        conservation_checks,
        [tspan[1], tspan[2]],
        Δt_cpl,
        Ref(tspan[1]),
        prev_checkpoint_t,
        model_sims,
        callbacks,
        dir_paths,
        thermo_params,
        diags_handler,
        save_cache,
    )

    #=
    ## Restart component model states if specified
    If a restart directory is specified and contains output files from the `checkpoint_cb` callback,
    the component model states are restarted from those files. The restart directory is specified in
    the `config_dict` dictionary. The `restart_t` field specifies the time step at which the restart
    is performed.

    If `restart_cache` is true, the caches will be read from the restart file using `restore_cache!`.
    Otherwise, the caches will be initialized in each component model's constructor.
    When the caches are not read from the restart file, we have to perform the initial component
    model exchange so that `set_caches!` can be called to initialize the caches.
    =#
    if detect_restart_files
        isnothing(restart_t) &&
            (restart_t = Checkpointer.t_start_from_checkpoint(dir_paths.checkpoints_dir))
        isnothing(restart_dir) && (restart_dir = dir_paths.checkpoints_dir)
    end
    should_restart = !isnothing(restart_t) && !isnothing(restart_dir)
    should_restart && Checkpointer.restart!(cs, restart_dir, restart_t, restart_cache)

    # Make sure surface model area fractions sum to 1 everywhere.
    # Note that ocean and ice fractions are not accurate until after this call.
    # Area fractions are not saved/read in when restarting, so we need to update them here
    # whether or not we restart.
    FieldExchanger.update_surface_fractions!(cs)

    if !should_restart || !restart_cache
        #=
        ## Initialize Component Model Exchange

        The concrete steps for proper initialization are:
        =#

        # 1. Import static fields into the coupler fields
        FieldExchanger.import_static_fields!(cs.fields, cs.model_sims)

        # 2. Import atmospheric and surface fields into the coupler fields,
        #  then broadcast them back out to all components.
        FieldExchanger.exchange!(cs)

        # 3. Update any fields in the model caches that can only be filled after the initial exchange.
        FieldExchanger.set_caches!(cs)

        # 4. Calculate and update turbulent fluxes for each surface model,
        #  and save the weighted average in coupler fields
        FluxCalculator.turbulent_fluxes!(cs)

        # 4. Compute any ocean-sea ice fluxes
        FluxCalculator.ocean_seaice_fluxes!(cs)
    end
    Utilities.show_memory_usage()
    return cs
end

"""
    setup_and_run(config_dict)
    setup_and_run(config_file = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml"))

This function sets up and runs the coupled model simulation specified by the
input config file or dict. It initializes the component models, all coupler objects,
diagnostics, and conservation checks, and then runs the simulation.
"""
function setup_and_run(
    config_file = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml"),
)
    cs = CoupledSimulation(config_file)
    run!(cs)
    return cs
end

function setup_and_run(config_dict)
    cs = CoupledSimulation(config_dict)
    run!(cs)
    return cs
end
