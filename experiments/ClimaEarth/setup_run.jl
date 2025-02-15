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
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaCore as CC

# ## Coupler specific imports
import ClimaCoupler
import ClimaCoupler:
    ConservationChecker, Checkpointer, FieldExchanger, FluxCalculator, Interfacer, TimeManager, Utilities
import ClimaCoupler.Interfacer:
    AbstractSlabplanetSimulationMode,
    AMIPMode,
    SlabplanetAquaMode,
    SlabplanetEisenmanMode,
    SlabplanetMode,
    SlabplanetTerraMode

import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaUtilities.TimeVaryingInputs: TimeVaryingInput, evaluate!
import ClimaUtilities.Utils: period_to_seconds_float
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import Interpolations # triggers InterpolationsExt in ClimaUtilities
# Random is used by RRMTGP for some cloud properties
import Random

# TODO: Move to ClimaUtilities once we move the Schedules to ClimaUtilities
import ClimaDiagnostics.Schedules: EveryCalendarDtSchedule

pkg_dir = pkgdir(ClimaCoupler)

#=
### Helper Functions
These will be eventually moved to their respective component model and utility packages, and so they should not
contain any internals of the ClimaCoupler source code, except extensions to the Interfacer functions.
=#

## helpers for component models
include("components/atmosphere/climaatmos.jl")
include("components/land/climaland_bucket.jl")
include("components/ocean/slab_ocean.jl")
include("components/ocean/prescr_ocean.jl")
include("components/ocean/prescr_seaice.jl")
include("components/ocean/eisenman_seaice.jl")

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
    solve_coupler!(cs)

This function performs the coupling loop, which is the main part of the simulation.
It runs the component models sequentially for one coupling timestep (`Δt_cpl`) at a time,
and exchanges combined fields and calculates fluxes using the selected turbulent fluxes option.
Note that we want to implement this in a dispatchable function to allow for
other forms of timestepping (e.g. leapfrog).
"""

function solve_coupler!(cs)
    (; model_sims, Δt_cpl, tspan, comms_ctx) = cs
    (; atmos_sim, land_sim, ocean_sim, ice_sim) = model_sims

    @info("Starting coupling loop")
    ## step in time
    for t in ((tspan[begin] + Δt_cpl):Δt_cpl:tspan[end])
        # Update date
        cs.dates.date[] = TimeManager.current_date(cs, t)

        ## compute global energy and water conservation checks
        ## (only for slabplanet if tracking conservation is enabled)
        !isnothing(cs.conservation_checks) && ConservationChecker.check_conservation!(cs)
        ClimaComms.barrier(comms_ctx)

        ## update water albedo from wind at dt_water_albedo
        ## (this will be extended to a radiation callback from the coupler)
        TimeManager.maybe_trigger_callback(cs.callbacks.water_albedo, cs, t)

        ## update the surface fractions for surface models,
        ## and update all component model simulations with the current fluxes stored in the coupler
        FieldExchanger.update_surface_fractions!(cs)
        FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

        ## step component model simulations sequentially for one coupling timestep (Δt_cpl)
        FieldExchanger.step_model_sims!(cs.model_sims, t)

        ## update the coupler with the new surface properties and calculate the turbulent fluxes
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

        ## update the coupler with the new atmospheric properties
        FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes) # radiative and/or turbulent

        ## callback to checkpoint model state
        TimeManager.maybe_trigger_callback(cs.callbacks.checkpoint, cs, t)

        ## compute/output AMIP diagnostics if scheduled for this timestep
        ## wrap the current CoupledSimulation fields and time in a NamedTuple to match the ClimaDiagnostics interface
        cs_nt = (; u = cs.fields, p = nothing, t = t, step = round(t / Δt_cpl))
        !isnothing(cs.diags_handler) && CD.orchestrate_diagnostics(cs_nt, cs.diags_handler)
    end
    return nothing
end

"""
    setup_and_run(config_file = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml"))

This function sets up and runs the coupled model simulation specified by the
input config file. It initializes the component models, all coupler objects,
diagnostics, and conservation checks, and then runs the simulation.
"""
function setup_and_run(config_file = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml"))
    # Parse the configuration file
    config_dict = get_coupler_config_dict(config_file)
    # Initialize communication context (do this first so all printing is only on root)
    comms_ctx = Utilities.get_comms_context(config_dict)
    # Select the correct timestep for each component model based on which are available
    parse_component_dts!(config_dict)
    # Add extra diagnostics if specified
    add_extra_diagnostics!(config_dict)

    (;
        job_id,
        sim_mode,
        random_seed,
        FT,
        t_end,
        t_start,
        date0,
        date,
        Δt_cpl,
        component_dt_dict,
        saveat,
        checkpoint_dt,
        restart_dir,
        restart_t,
        use_coupler_diagnostics,
        use_land_diagnostics,
        calendar_dt,
        evolving_ocean,
        mono_surface,
        turb_flux_partition,
        land_domain_type,
        land_albedo_type,
        land_initial_condition,
        land_temperature_anomaly,
        energy_check,
        conservation_softfail,
        output_dir_root,
    ) = get_coupler_args(config_dict)

    #=
    ### I/O Directory Setup `setup_output_dirs` returns `dir_paths.output =
    COUPLER_OUTPUT_DIR`, which is the directory where the output of the simulation
    will be saved, `dir_paths.artifacts` is the directory where the plots (from
    postprocessing and the conservation checks) of the simulation will be saved,
    #and `dir_paths.checkpoints`, where restart files are saved.
    =#

    COUPLER_OUTPUT_DIR = joinpath(output_dir_root, job_id)
    dir_paths = Utilities.setup_output_dirs(output_dir = COUPLER_OUTPUT_DIR, comms_ctx = comms_ctx)
    @info "Coupler output directory $(dir_paths.output)"
    @info "Coupler artifacts directory $(dir_paths.artifacts)"
    @info "Coupler checkpoint directory $(dir_paths.checkpoints)"

    atmos_output_dir = joinpath(dir_paths.output, "clima_atmos")
    isdir(atmos_output_dir) || mkpath(atmos_output_dir)
    land_output_dir = joinpath(dir_paths.output, "clima_land")
    isdir(land_output_dir) || mkpath(land_output_dir)


    ## get component model dictionaries (if applicable)
    ## Note this step must come after parsing the coupler config dictionary, since
    ##  some parameters are passed from the coupler config to the component model configs
    atmos_config_dict = get_atmos_config_dict(config_dict, job_id, atmos_output_dir)
    (; dt_rad, output_default_diagnostics) = get_atmos_args(atmos_config_dict)

    ## set unique random seed if desired, otherwise use default
    Random.seed!(random_seed)
    @info "Random seed set to $(random_seed)"

    tspan = (t_start, t_end)

    #=
    ## Data File Paths
    =#
    land_mask_data = joinpath(@clima_artifact("landsea_mask_60arcseconds", comms_ctx), "landsea_mask.nc")

    #=
    ## Component Model Initialization
    Here we set initial and boundary conditions for each component model. Each component model is required to have an `init` function that
    returns a `ComponentModelSimulation` object (see `Interfacer` docs for more details).
    =#

    #=
    ### Atmosphere
    This uses the `ClimaAtmos.jl` model, with parameterization options specified in the `atmos_config_object` dictionary.
    =#

    Utilities.show_memory_usage()

    ## init atmos model component
    atmos_sim = ClimaAtmosSimulation(CA.AtmosConfig(atmos_config_dict))
    # Get surface elevation from `atmos` coordinate field
    surface_elevation = CC.Fields.level(CC.Fields.coordinate_field(atmos_sim.integrator.u.f).z, CC.Utilities.half)
    Utilities.show_memory_usage()

    thermo_params = get_thermo_params(atmos_sim) # TODO: this should be shared by all models #342

    #=
    ### Boundary Space
    We use a common `Space` for all global surfaces. This enables the MPI processes to operate on the same columns in both
    the atmospheric and surface components, so exchanges are parallelized. Note this is only possible when the
    atmosphere and surface are of the same horizontal resolution.

    Currently, we use the 2D surface space from the atmosphere model as our shared space,
    but ultimately we want this to specified within the coupler and passed to all component models. (see issue #665)
    =#

    ## init a 2D boundary space at the surface
    boundary_space = CC.Spaces.horizontal_space(atmos_sim.domain.face_space) # TODO: specify this in the coupler and pass it to all component models #665

    #=
    ### Land-sea Fraction
    This is a static field that contains the area fraction of land and sea, ranging from 0 to 1.
    If applicable, sea ice is included in the sea fraction at this stage.
    Note that land-sea area fraction is different to the land-sea mask, which is a binary field
    (masks are used internally by the coupler to indicate passive cells that are not populated by a given component model).
    =#

    # Preprocess the file to be 1s and 0s before remapping into onto the grid
    land_fraction = SpaceVaryingInput(land_mask_data, "landsea", boundary_space)
    if !mono_surface
        land_fraction = Utilities.binary_mask.(land_fraction)
    end
    Utilities.show_memory_usage()

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

    In this section of the code, we initialize all component models and read in the prescribed data we'll be using.
    The specific models and data that are set up depend on which mode we're running.
    =#

    @info(sim_mode)
    if sim_mode <: AMIPMode
        @info("AMIP boundary conditions - do not expect energy conservation")

        ## land model
        land_sim = BucketSimulation(
            FT,
            tspan,
            land_domain_type,
            land_albedo_type,
            land_initial_condition,
            land_temperature_anomaly,
            land_output_dir;
            dt = component_dt_dict["dt_land"],
            space = boundary_space,
            saveat = saveat,
            area_fraction = land_fraction,
            date_ref = date0,
            t_start = t_start,
            energy_check = energy_check,
            surface_elevation,
            use_land_diagnostics,
        )

        ## sea ice model
        ice_sim = PrescribedIceSimulation(
            FT;
            tspan = tspan,
            dt = component_dt_dict["dt_seaice"],
            saveat = saveat,
            space = boundary_space,
            thermo_params = thermo_params,
            comms_ctx,
            date0,
            mono_surface,
            land_fraction,
        )

        ## ocean model using prescribed data
        ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction))
        ocean_fraction = FT(1) .- ice_fraction .- land_fraction
        ocean_sim =
            PrescribedOceanSimulation(FT, boundary_space, date0, t_start, ocean_fraction, thermo_params, comms_ctx)

        Utilities.show_memory_usage()

    elseif (sim_mode <: AbstractSlabplanetSimulationMode) && !(sim_mode <: SlabplanetEisenmanMode)


        land_fraction = sim_mode <: SlabplanetAquaMode ? land_fraction .* 0 : land_fraction
        land_fraction = sim_mode <: SlabplanetTerraMode ? land_fraction .* 0 .+ 1 : land_fraction

        ## land model
        land_sim = BucketSimulation(
            FT,
            tspan,
            land_domain_type,
            land_albedo_type,
            land_initial_condition,
            land_temperature_anomaly,
            land_output_dir;
            dt = component_dt_dict["dt_land"],
            space = boundary_space,
            saveat = saveat,
            area_fraction = land_fraction,
            date_ref = date0,
            t_start = t_start,
            energy_check = energy_check,
            surface_elevation,
            use_land_diagnostics,
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

        ## sea ice stub (here set to zero area coverage)
        ice_sim = Interfacer.SurfaceStub((;
            T_sfc = ones(boundary_space),
            ρ_sfc = zeros(boundary_space),
            z0m = FT(0),
            z0b = FT(0),
            beta = FT(1),
            α_direct = ones(boundary_space),
            α_diffuse = ones(boundary_space),
            area_fraction = zeros(boundary_space),
            phase = TD.Ice(),
            thermo_params = thermo_params,
        ))

        Utilities.show_memory_usage()

    elseif sim_mode <: SlabplanetEisenmanMode

        ## land model
        land_sim = BucketSimulation(
            FT,
            tspan,
            land_domain_type,
            land_albedo_type,
            land_initial_condition,
            land_temperature_anomaly,
            land_output_dir;
            dt = component_dt_dict["dt_land"],
            space = boundary_space,
            saveat = saveat,
            area_fraction = land_fraction,
            date_ref = date0,
            t_start = t_start,
            energy_check = energy_check,
            surface_elevation,
            use_land_diagnostics,
        )

        ## ocean stub (here set to zero area coverage)
        ocean_sim = SlabOceanSimulation(
            FT;
            tspan = tspan,
            dt = component_dt_dict["dt_ocean"],
            space = boundary_space,
            saveat = saveat,
            area_fraction = zeros(boundary_space), # zero, since ML is calculated below
            thermo_params = thermo_params,
        )

        ## sea ice + ocean model
        ice_sim = EisenmanIceSimulation(
            FT,
            tspan,
            space = boundary_space,
            area_fraction = (FT(1) .- land_fraction),
            dt = component_dt_dict["dt_seaice"],
            saveat = saveat,
            thermo_params = thermo_params,
        )

        Utilities.show_memory_usage()
    end

    #=
    ## Coupler Initialization
    The coupler needs to contain exchange information, access all component models, and manage the calendar,
    among other responsibilities.
    Objects containing information to enable these are initialized here and saved in the
    global `CoupledSimulation` struct, `cs`, below.
    =#

    ## coupler exchange fields
    coupler_field_names = (
        :T_sfc,
        :z0m_sfc,
        :z0b_sfc,
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
    coupler_fields = NamedTuple{coupler_field_names}(ntuple(i -> zeros(boundary_space), length(coupler_field_names)))
    Utilities.show_memory_usage()

    ## model simulations
    model_sims = (atmos_sim = atmos_sim, ice_sim = ice_sim, land_sim = land_sim, ocean_sim = ocean_sim)

    ## dates
    dates = (; date = [date], date0 = [date0])

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
            sim_mode <: AbstractSlabplanetSimulationMode && !CA.is_distributed(ClimaComms.context(boundary_space)),
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

    The currently implemented callbacks are:
    - `checkpoint_cb`: generates a checkpoint of all model states at a specified interval. This is mainly used for restarting simulations.
    - `albedo_cb`: for the amip mode, the water albedo is time varying (since the reflectivity of water depends on insolation and wave characteristics, with the latter
    being approximated from wind speed). It is updated at the same frequency as the atmospheric radiation.
    NB: Eventually, we will call all of radiation from the coupler, in addition to the albedo calculation.
    =#
    schedule_checkpoint = EveryCalendarDtSchedule(TimeManager.time_to_period(checkpoint_dt); start_date = date0)
    checkpoint_cb = TimeManager.TimeManager.Callback(schedule_checkpoint, Checkpointer.checkpoint_sims)

    if sim_mode <: AMIPMode
        schedule_albedo = EveryCalendarDtSchedule(TimeManager.time_to_period(dt_rad); start_date = date0)
    else
        schedule_albedo = TimeManager.NeverSchedule()
    end
    albedo_cb = TimeManager.Callback(schedule_albedo, FluxCalculator.water_albedo_from_atmosphere!)

    callbacks = (; checkpoint = checkpoint_cb, water_albedo = albedo_cb)


    #=
    ## Initialize turbulent fluxes

    Decide on the type of turbulent flux partition, partitioned or combined (see `FluxCalculator` documentation for more details).
    =#
    turbulent_fluxes = nothing
    if turb_flux_partition == "PartitionedStateFluxes"
        turbulent_fluxes = FluxCalculator.PartitionedStateFluxes()
    elseif turb_flux_partition == "CombinedStateFluxesMOST"
        turbulent_fluxes = FluxCalculator.CombinedStateFluxesMOST()
    else
        error("turb_flux_partition must be either PartitionedStateFluxes or CombinedStateFluxesMOST")
    end

    #= Set up default AMIP diagnostics
    Use ClimaDiagnostics for default AMIP diagnostics, which currently include turbulent energy fluxes.
    =#
    if use_coupler_diagnostics
        @info "Using default coupler diagnostics"
        coupler_diags_path = joinpath(dir_paths.output, "coupler")
        isdir(coupler_diags_path) || mkpath(coupler_diags_path)
        diags_handler =
            coupler_diagnostics_setup(coupler_fields, coupler_diags_path, dates.date0[1], tspan[1], calendar_dt)
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

    cs = Interfacer.CoupledSimulation{FT}(
        comms_ctx,
        dates,
        boundary_space,
        coupler_fields,
        conservation_checks,
        [tspan[1], tspan[2]],
        Δt_cpl,
        model_sims,
        callbacks,
        dir_paths,
        turbulent_fluxes,
        thermo_params,
        diags_handler,
    )
    Utilities.show_memory_usage()

    #=
    ## Restart component model states if specified
    If a restart directory is specified and contains output files from the `checkpoint_cb` callback, the component model states are restarted from those files. The restart directory
    is specified in the `config_dict` dictionary. The `restart_t` field specifies the time step at which the restart is performed.
    =#

    if !isnothing(restart_dir)
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
    FieldExchanger.update_surface_fractions!(cs)

    # 2.surface density (`ρ_sfc`): calculated by the coupler by adiabatically extrapolating atmospheric thermal state to the surface.
    # For this, we need to import surface and atmospheric fields. The model sims are then updated with the new surface density.
    FieldExchanger.import_combined_surface_fields!(cs.fields, cs.model_sims, cs.turbulent_fluxes)
    FieldExchanger.import_atmos_fields!(cs.fields, cs.model_sims, cs.boundary_space, cs.turbulent_fluxes)
    FieldExchanger.update_model_sims!(cs.model_sims, cs.fields, cs.turbulent_fluxes)

    # 3.surface vapor specific humidity (`q_sfc`): step surface models with the new surface density to calculate their respective `q_sfc` internally
    ## TODO: the q_sfc calculation follows the design of the bucket q_sfc, but it would be neater to abstract this from step! (#331)
    Interfacer.step!(land_sim, tspan[1] + Δt_cpl)
    Interfacer.step!(ocean_sim, tspan[1] + Δt_cpl)
    Interfacer.step!(ice_sim, tspan[1] + Δt_cpl)

    # 4.turbulent fluxes: now we have all information needed for calculating the initial turbulent
    # surface fluxes using either the combined state or the partitioned state method
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
    ## Precompilation of Coupling Loop

    Here we run the entire coupled simulation for two timesteps to precompile everything
    for accurate timing of the overall simulation. After these two steps, we update the
    beginning and end of the simulation timespan to the correct values.
    =#

    ## run the coupled simulation for two timesteps to precompile
    cs.tspan[2] = tspan[1] + Δt_cpl * 2
    solve_coupler!(cs)

    ## update the timespan to the correct values
    cs.tspan[1] = tspan[1] + Δt_cpl * 2
    cs.tspan[2] = tspan[2]

    ## Run garbage collection before solving for more accurate memory comparison to ClimaAtmos
    GC.gc()

    #=
    ## Solving and Timing the Full Simulation

    This is where the full coupling loop, `solve_coupler!` is called for the full timespan of the simulation.
    We use the `ClimaComms.@elapsed` macro to time the simulation on both CPU and GPU, and use this
    value to calculate the simulated years per day (SYPD) of the simulation.
    =#
    walltime = ClimaComms.@elapsed comms_ctx.device begin
        s = CA.@timed_str begin
            solve_coupler!(cs)
        end
    end
    @info(walltime)

    ## Use ClimaAtmos calculation to show the simulated years per day of the simulation (SYPD)
    es = CA.EfficiencyStats((tspan[1], tspan[2]), walltime)
    sypd = CA.simulated_years_per_day(es)
    n_atmos_steps = atmos_sim.integrator.step
    walltime_per_atmos_step = es.walltime / n_atmos_steps
    @info "SYPD: $sypd"
    @info "Walltime per Atmos step: $(walltime_per_atmos_step)"

    ## Save the SYPD and allocation information
    if ClimaComms.iamroot(comms_ctx)
        open(joinpath(dir_paths.artifacts, "sypd.txt"), "w") do sypd_filename
            println(sypd_filename, "$sypd")
        end

        open(joinpath(dir_paths.artifacts, "walltime_per_atmos_step.txt"), "w") do walltime_per_atmos_step_filename
            println(walltime_per_atmos_step_filename, "$(walltime_per_atmos_step)")
        end

        open(joinpath(dir_paths.artifacts, "max_rss_cpu.txt"), "w") do cpu_max_rss_filename
            cpu_max_rss_GB = Utilities.show_memory_usage()
            println(cpu_max_rss_filename, cpu_max_rss_GB)
        end
    end

    #=
    ## Postprocessing
    All postprocessing is performed using the root process only, if applicable.
    Our postprocessing consists of outputting a number of plots to visualize the model output.

    The postprocessing includes:
    - Energy and water conservation checks (if running SlabPlanet with checks enabled)
    - Animations (if not running in MPI)
    - AMIP plots of the final state of the model
    - Error against observations
    - Optional additional atmosphere diagnostics plots
    - Plots of useful coupler and component model fields for debugging
    =#

    if ClimaComms.iamroot(comms_ctx) && use_coupler_diagnostics
        postprocessing_vars = (; output_default_diagnostics, t_end, conservation_softfail)
        postprocess_sim(cs, postprocessing_vars)
    end

    # Close all diagnostics file writers
    !isnothing(cs.diags_handler) && map(diag -> close(diag.output_writer), cs.diags_handler.scheduled_diagnostics)
end
