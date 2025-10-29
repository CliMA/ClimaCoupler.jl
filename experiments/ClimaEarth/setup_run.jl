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

# ## Coupler specific imports
import ClimaCoupler
import ClimaCoupler:
    ConservationChecker,
    Checkpointer,
    FieldExchanger,
    FluxCalculator,
    Interfacer,
    TimeManager,
    Utilities
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
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaUtilities.TimeManager: ITime, date
import Interpolations # triggers InterpolationsExt in ClimaUtilities
# Random is used by RRMTGP for some cloud properties
import Random

# TODO: Move to ClimaUtilities once we move the Schedules to ClimaUtilities
import ClimaDiagnostics.Schedules: EveryCalendarDtSchedule, EveryStepSchedule

pkg_dir = pkgdir(ClimaCoupler)

#=
### Helper Functions
These will be eventually moved to their respective component model and utility packages, and so they should not
contain any internals of the ClimaCoupler source code, except extensions to the Interfacer functions.
=#

## helpers for component models
include("components/atmosphere/climaatmos.jl")
include("components/land/climaland_bucket.jl")
include("components/land/climaland_integrated.jl")
include("components/ocean/slab_ocean.jl")
include("components/ocean/prescr_ocean.jl")
include("components/ocean/prescr_seaice.jl")
include("components/ocean/oceananigans.jl")

#=
### Configuration Dictionaries
Each simulation mode has its own configuration dictionary. The `config_dict` of each simulation is a merge of the default configuration
dictionary and the simulation-specific configuration dictionary, which allows the user to override the default settings.

We can additionally pass the configuration dictionary to the component model initializers, which will then override the default settings of the component models.
=#

include("cli_options.jl")
include("user_io/arg_parsing.jl")
include("user_io/postprocessing.jl")
include("user_io/coupler_diagnostics.jl")

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
    config_dict = get_coupler_config_dict(config_file)
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
        saveat,
        checkpoint_dt,
        detect_restart_files,
        restart_dir,
        restart_t,
        use_land_diagnostics,
        diagnostics_dt,
        evolving_ocean,
        land_model,
        land_temperature_anomaly,
        land_spun_up_ic,
        bucket_albedo_type,
        bucket_initial_condition,
        energy_check,
        use_coupler_diagnostics,
        output_dir_root,
        parameter_files,
        era5_initial_condition_dir,
    ) = get_coupler_args(config_dict)

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

    ## get component model dictionaries (if applicable)
    ## Note this step must come after parsing the coupler config dictionary, since
    ##  some parameters are passed from the coupler config to the component model configs
    atmos_config_dict = get_atmos_config_dict(config_dict, dir_paths.atmos_output_dir)

    ## set unique random seed if desired, otherwise use default
    Random.seed!(random_seed)
    @info "Random seed set to $(random_seed)"

    if detect_restart_files
        isnothing(restart_t) &&
            (restart_t = Checkpointer.t_start_from_checkpoint(dir_paths.checkpoints_dir))
        isnothing(restart_dir) && (restart_dir = dir_paths.checkpoints_dir)
    end
    should_restart = !isnothing(restart_t) && !isnothing(restart_dir)
    if should_restart
        if t_start isa ITime
            t_start, _ = promote(ITime(restart_t), t_start)
        else
            t_start = restart_t
        end

        if pkgversion(CA) >= v"0.29.1"
            # We only support a round number of seconds
            isinteger(float(t_start)) ||
                error("Cannot restart from a non integer number of seconds")
            t_start_int = Int(float(t_start))
            atmos_config_dict.parsed_args["t_start"] = "$(t_start_int)secs"
        else
            # There was no `t_start`, so we have to use a workaround for this.
            # This does not support passing the command-line arguments (unless
            # restart_dir is exactly the same as output_dir_root)
            atmos_config_dict.parsed_args["restart_file"] =
                climaatmos_restart_path(output_dir_root, restart_t)
        end

        @info "Starting from t_start $(t_start)"
    end

    tspan = (t_start, t_end)

    #=
    ## Component Model Initialization
    Here we set initial and boundary conditions for each component model. Each component model is required to have an `init` function that
    returns a `ComponentModelSimulation` object (see `Interfacer` docs for more details).
    =#

    #=
    ### Atmosphere
    This uses the `ClimaAtmos.jl` model, with parameterization options specified in the `atmos_config_object` dictionary.
    =#

    ## init atmos model component
    atmos_sim = ClimaAtmosSimulation(atmos_config_dict)
    thermo_params = get_thermo_params(atmos_sim) # TODO: this should be shared by all models #342

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
        h_elem = config_dict["h_elem"]
        n_quad_points = 4
        boundary_space =
            CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points, h_elem)
    end

    # Get surface elevation on the boundary space from `atmos` coordinate field
    surface_elevation = Interfacer.get_field(atmos_sim, Val(:height_sfc), boundary_space) # on surface space

    # Get atmospheric height on the boundary space
    # Since the atmospheric height is defined on centers, we need to copy the values onto the boundary space
    # to be able to subtract the surface elevation
    # Note: This pattern is not reliable and should not be reused.
    atmos_h =
        ClimaCore.Fields.Field(
            ClimaCore.Fields.field_values(
                Interfacer.get_field(atmos_sim, Val(:height_int)),
            ),
            boundary_space,
        ) .- surface_elevation # atmos height relative to the surface, on the surface space

    #=
    ### Land-sea Fraction
    This is a static field that contains the area fraction of land and sea, ranging from 0 to 1.
    If applicable, sea ice is included in the sea fraction at this stage.
    Note that land-sea area fraction is different to the land-sea mask, which is a binary field
    (masks are used internally by the coupler to indicate passive cells that are not populated by a given component model).
    =#

    # Preprocess the file to be 1s and 0s before remapping into onto the grid
    land_mask_data =
        joinpath(@clima_artifact("landsea_mask_60arcseconds", comms_ctx), "landsea_mask.nc")
    land_fraction = SpaceVaryingInput(land_mask_data, "landsea", boundary_space)
    land_fraction = ifelse.(land_fraction .> eps(FT), FT(1), FT(0))

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
    if sim_mode <: AMIPMode || sim_mode <: CMIPMode || sim_mode <: SubseasonalMode
        @info("AMIP/CMIP boundary conditions - do not expect energy conservation")

        # Build ERA5-based file paths if subseasonal mode is selected
        subseasonal_sst = subseasonal_sic = subseasonal_land_ic = nothing
        if sim_mode <: SubseasonalMode
            isnothing(era5_initial_condition_dir) &&
                error("subseasonal mode requires --era5_initial_condition_dir")
            # Filenames inferred from start_date, which is YYYYMMDD
            datestr = Dates.format(start_date, Dates.dateformat"yyyymmdd")
            subseasonal_sst =
                joinpath(era5_initial_condition_dir, "sst_processed_$(datestr)_0000.nc")
            subseasonal_sic =
                joinpath(era5_initial_condition_dir, "sic_processed_$(datestr)_0000.nc")
            subseasonal_land_ic = joinpath(
                era5_initial_condition_dir,
                "era5_land_processed_$(datestr)_0000.nc",
            )
        end

        ## land model
        # Determine whether to use a shared surface space
        shared_surface_space = share_surface_space ? boundary_space : nothing
        if land_model == "bucket"
            land_sim = BucketSimulation(
                FT;
                dt = component_dt_dict["dt_land"],
                tspan,
                start_date,
                output_dir = dir_paths.land_output_dir,
                area_fraction = land_fraction,
                shared_surface_space,
                surface_elevation,
                land_temperature_anomaly,
                use_land_diagnostics,
                albedo_type = bucket_albedo_type,
                bucket_initial_condition,
                energy_check,
                parameter_files,
            )
        elseif land_model == "integrated"
            land_sim = ClimaLandSimulation(
                FT;
                dt = component_dt_dict["dt_land"],
                tspan,
                start_date,
                output_dir = dir_paths.land_output_dir,
                area_fraction = land_fraction,
                shared_surface_space,
                land_spun_up_ic,
                saveat,
                surface_elevation,
                atmos_h,
                land_temperature_anomaly,
                use_land_diagnostics,
                parameter_files,
                land_ic_path = subseasonal_land_ic,
            )
        else
            error("Invalid land model specified: $(land_model)")
        end

        ## sea ice model
        ice_sim = PrescribedIceSimulation(
            FT;
            tspan = tspan,
            dt = component_dt_dict["dt_seaice"],
            saveat = saveat,
            space = boundary_space,
            thermo_params = thermo_params,
            comms_ctx,
            start_date,
            land_fraction,
            sic_path = subseasonal_sic,
        )

        ## ocean model using prescribed data
        ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction))
        ocean_fraction = FT(1) .- ice_fraction .- land_fraction

        if sim_mode <: CMIPMode
            stop_date = date(tspan[end] - tspan[begin])
            ocean_sim = OceananigansSimulation(
                ocean_fraction,
                start_date,
                stop_date;
                output_dir = dir_paths.ocean_output_dir,
                comms_ctx,
            )
        else
            ocean_sim = PrescribedOceanSimulation(
                FT,
                boundary_space,
                start_date,
                t_start,
                ocean_fraction,
                thermo_params,
                comms_ctx;
                sst_path = subseasonal_sst,
            )
        end

    elseif (sim_mode <: AbstractSlabplanetSimulationMode)

        land_fraction = sim_mode <: SlabplanetAquaMode ? land_fraction .* 0 : land_fraction
        land_fraction =
            sim_mode <: SlabplanetTerraMode ? land_fraction .* 0 .+ 1 : land_fraction

        ## land model
        land_sim = BucketSimulation(
            FT;
            dt = component_dt_dict["dt_land"],
            tspan,
            start_date,
            output_dir = dir_paths.land_output_dir,
            area_fraction = land_fraction,
            surface_elevation,
            land_temperature_anomaly,
            use_land_diagnostics,
            albedo_type = bucket_albedo_type,
            bucket_initial_condition,
            energy_check,
            parameter_files,
        )

        ## ocean model
        ocean_sim = SlabOceanSimulation(
            FT;
            tspan = tspan,
            dt = component_dt_dict["dt_ocean"],
            space = boundary_space,
            saveat = saveat,
            area_fraction = (FT(1) .- land_fraction), ## NB: this ocean fraction includes areas covered by sea ice (unlike the one contained in the cs)
            thermo_params = thermo_params,
            evolving = evolving_ocean,
        )
    end

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
    checkpoint_cb = TimeManager.Callback(schedule_checkpoint, Checkpointer.checkpoint_sims)

    callbacks = (checkpoint_cb,)

    #= Set up default AMIP diagnostics
    Use ClimaDiagnostics for default AMIP diagnostics, which currently include turbulent energy fluxes.
    =#
    if use_coupler_diagnostics
        @info "Using default coupler diagnostics"
        diags_handler = coupler_diagnostics_setup(
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
    )

    #=
    ## Restart component model states if specified
    If a restart directory is specified and contains output files from the `checkpoint_cb` callback, the component model states are restarted from those files. The restart directory
    is specified in the `config_dict` dictionary. The `restart_t` field specifies the time step at which the restart is performed.
    =#
    should_restart && Checkpointer.restart!(cs, restart_dir, restart_t)

    if !should_restart
        #=
        ## Initialize Component Model Exchange

        The concrete steps for proper initialization are:
        =#

        # 1. Make sure surface model area fractions sum to 1 everywhere.
        FieldExchanger.update_surface_fractions!(cs)

        # 2. Import atmospheric and surface fields into the coupler fields,
        #  then broadcast them back out to all components.
        FieldExchanger.exchange!(cs)

        # 3. Update any fields in the model caches that can only be filled after the initial exchange.
        FieldExchanger.set_caches!(cs)

        # 4. Calculate and update turbulent fluxes for each surface model,
        #  and save the weighted average in coupler fields
        FluxCalculator.turbulent_fluxes!(cs)
    end
    Utilities.show_memory_usage()
    return cs
end

"""
    run!(cs::CoupledSimulation)

Evolve the given simulation, producing plots and other diagnostic information.

Keyword arguments
==================

`precompile`: If `true`, run the coupled simulations for two steps, so that most functions
              are precompiled and subsequent timing will be more accurate.
"""
function run!(
    cs::CoupledSimulation;
    precompile = (cs.tspan[end] > 2 * cs.Δt_cpl + cs.tspan[begin]),
)

    ## Precompilation of Coupling Loop
    # Here we run the entire coupled simulation for two timesteps to precompile several
    # functions for more accurate timing of the overall simulation.
    precompile && (step!(cs); step!(cs))

    ## Run garbage collection before solving for more accurate memory comparison to ClimaAtmos
    GC.gc()

    #=
    ## Solving and Timing the Full Simulation

    This is where the full coupling loop, `solve_coupler!` is called for the full timespan of the simulation.
    We use the `ClimaComms.@elapsed` macro to time the simulation on both CPU and GPU, and use this
    value to calculate the simulated years per day (SYPD) of the simulation.
    =#
    @info "Starting coupling loop"
    walltime = ClimaComms.@elapsed ClimaComms.device(cs) begin
        s = CA.@timed_str begin
            while cs.t[] <= cs.tspan[end]
                step!(cs)
            end
        end
    end
    @info "Simulation took $(walltime) seconds"

    sypd = simulated_years_per_day(cs, walltime)
    walltime_per_step = walltime_per_coupling_step(cs, walltime)
    @info "SYPD: $sypd"
    @info "Walltime per coupling step: $(walltime_per_step)"
    save_sypd_walltime_to_disk(cs, walltime)

    # Close all diagnostics file writers
    isnothing(cs.diags_handler) ||
        foreach(diag -> close(diag.output_writer), cs.diags_handler.scheduled_diagnostics)
    foreach(Interfacer.close_output_writers, cs.model_sims)

    return nothing
end

"""
    postprocess(cs; conservation_softfail = false, rmse_check = false)

Process the results after a simulation has completed, including generating
plots, checking conservation, and other diagnostics.
All postprocessing is performed using the root process only, if applicable.

When `conservation_softfail` is true, throw an error if conservation is not
respected.

When `rmse_check` is true, compute the RMSE against observations and test
that it is below a certain threshold.

The postprocessing includes:
- Energy and water conservation checks (if running SlabPlanet with checks enabled)
- Animations (if not running in MPI)
- AMIP plots of the final state of the model
- Error against observations
- Optional additional atmosphere diagnostics plots
- Plots of useful coupler and component model fields for debugging
"""
function postprocess(cs; conservation_softfail = false, rmse_check = false)
    if ClimaComms.iamroot(ClimaComms.context(cs)) && !isnothing(cs.diags_handler)
        postprocessing_vars = (; conservation_softfail, rmse_check)
        postprocess_sim(cs, postprocessing_vars)
    end
    return nothing
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

"""
    step!(cs::CoupledSimulation)

Take one coupling step forward in time.

This function runs the component models sequentially, and exchanges combined fields and
calculates fluxes using the selected turbulent fluxes option. Note, one coupling step might
require multiple steps in some of the component models.
"""
function step!(cs::CoupledSimulation)
    # Update the current time
    cs.t[] += cs.Δt_cpl

    ## compute global energy and water conservation checks
    ## (only for slabplanet if tracking conservation is enabled)
    ConservationChecker.check_conservation!(cs)

    ## step component model simulations sequentially for one coupling timestep (Δt_cpl)
    FieldExchanger.step_model_sims!(cs)

    ## update the surface fractions for surface models
    FieldExchanger.update_surface_fractions!(cs)

    ## exchange all non-turbulent flux fields between models
    FieldExchanger.exchange!(cs)

    ## calculate turbulent fluxes in the coupler and update the model simulations with them
    FluxCalculator.turbulent_fluxes!(cs)

    ## Maybe call the callbacks
    TimeManager.callbacks!(cs)

    # Compute coupler diagnostics
    CD.orchestrate_diagnostics(cs)
    return nothing
end
