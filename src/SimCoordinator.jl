"""
    SimCoordinator

This module contains functions for coordinating the setup and execution of coupled
simulations. Key exports:

- `CoupledSimulation(config_file)` / `CoupledSimulation(config_dict)`: construct a
  fully initialized `CoupledSimulation` from a YAML config file path or dictionary.
  These are outer constructors for `Interfacer.CoupledSimulation`, defined here because
  they depend on modules loaded after `Interfacer` (`Input`, `Utilities`, etc.).
- `setup_and_run(config)`: convenience function that constructs a `CoupledSimulation`
  and immediately calls `run!` on it.
- `run!(cs)`: evolve a `CoupledSimulation` through its full time span.
- `step!(cs)`: advance a `CoupledSimulation` by one coupling time step.
"""
module SimCoordinator

import ClimaComms
import ClimaDiagnostics as CD
import ClimaUtilities.TimeManager: ITime
import ClimaCore as CC
import ClimaParams as CP
import Thermodynamics.Parameters as TDP
import Random

import ClimaCoupler
import ..Interfacer
import ..ConservationChecker
import ..FieldExchanger
import ..FluxCalculator
import ..TimeManager
import ..Input
import ..Utilities
import ..Checkpointer
import ..SimOutput

export run!, step!, setup_and_run

"""
    run!(cs::CoupledSimulation)

Evolve the given simulation, producing plots and other diagnostic information.

Keyword arguments
==================

`precompile`: If `true`, run the coupled simulations for two steps, so that most functions
              are precompiled and subsequent timing will be more accurate.
"""
function run!(
    cs::Interfacer.CoupledSimulation;
    precompile::Bool = cs.tspan[end] > 2 * cs.Δt_cpl + cs.tspan[begin],
)
    ## Precompilation of Coupling Loop
    # Here we run the entire coupled simulation for two timesteps to precompile several
    # functions for more accurate timing of the overall simulation.
    precompile && (step!(cs); step!(cs))

    ## Run garbage collection before solving for more accurate memory comparison to ClimaAtmos
    GC.gc()

    ## Solving and Timing the Full Simulation

    # This is where the full coupling loop is called for the full timespan of the simulation.
    # We use the `ClimaComms.@elapsed` macro to time the simulation on both CPU and GPU and use this
    # value to calculate the simulated years per day (SYPD) of the simulation.
    @info "Starting coupling loop"
    t_timed_start = cs.t[] # get t just before timing (equal to cs.tspan[begin] if `precompile` is false)
    walltime = ClimaComms.@elapsed ClimaComms.device(cs) begin
        while cs.t[] < cs.tspan[end]
            step!(cs)
        end
    end
    # Nothing may outlive the coupling loop with ice/ocean state in flight.
    FieldExchanger.wait_slow_sims!(cs)

    @info "Simulation took $(walltime) seconds"

    save_sypd_walltime_to_disk(cs, walltime, t_timed_start)

    # Close all diagnostics file writers
    isnothing(cs.diags_handler) ||
        foreach(diag -> close(diag.output_writer), cs.diags_handler.scheduled_diagnostics)
    foreach(Interfacer.close_output_writers, cs.model_sims)

    return nothing
end

"""
    step!(cs::CoupledSimulation)

Take one coupling step forward in time.

This function runs the component models sequentially, and exchanges combined fields and
calculates fluxes using the selected turbulent fluxes option. Note, one coupling step might
require multiple steps in some of the component models.
"""
function step!(cs::Interfacer.CoupledSimulation)
    # Update the current time and step number
    cs.t[] += cs.Δt_cpl
    cs.step[] += 1

    frozen = join_slow_if_due!(cs)

    # Compute global energy and water conservation checks
    # (only for slabplanet if tracking conservation is enabled)
    frozen || ConservationChecker.check_conservation!(cs)

    # Step component model simulations sequentially for one coupling timestep (Δt_cpl)
    FieldExchanger.step_model_sims!(cs; skip_slow = skip_slow_stepping(cs, frozen))

    # Update the surface fractions for surface models
    FieldExchanger.update_surface_fractions!(cs; slow_frozen = frozen)

    # Exchange all non-turbulent flux fields between models, including radiative and precipitation fluxes
    FieldExchanger.exchange!(cs; slow_frozen = frozen)

    # Calculate turbulent fluxes in the coupler and update the model simulations with them
    FluxCalculator.turbulent_fluxes!(
        cs;
        slow_frozen = frozen,
        force_slow_push = cs.prime_slow_surfaces && !frozen && at_slow_boundary(cs),
    )

    # Compute any ocean-sea ice fluxes
    frozen || FluxCalculator.ocean_seaice_fluxes!(cs)

    # The slow surfaces' forcing is now fully assembled, so their step can be
    # launched and left to run across the coupling steps that follow.
    launch_slow_if_due!(cs, frozen)

    # Maybe call the callbacks
    TimeManager.callbacks!(cs)

    # Compute and save coupler diagnostics
    CD.orchestrate_diagnostics(cs)
    return nothing
end

"""
    slow_step_due(cs)

Whether the overlapped group's step boundary falls on this coupling step.

Priming moves the component clocks away from coupler time, so the two schedules
ask different questions: without priming, whether the components are due to step;
with it, whether the coupler has reached the recorded window boundary.
"""
slow_step_due(cs::Interfacer.CoupledSimulation) =
    cs.prime_slow_surfaces ? at_slow_boundary(cs) : slow_surfaces_due(cs)

"""
    join_slow_if_due!(cs) -> frozen

Join an overlapped ice/ocean step if this coupling step needs its result, and
report whether one is still in flight afterwards.

Joining at the top of the step means the rest of it sees settled ice and ocean
state and behaves as it would without overlap. `frozen` is what every
slow-surface-touching call in `step!` keys off.
"""
function join_slow_if_due!(cs::Interfacer.CoupledSimulation)
    cs.overlap_slow_surfaces && slow_step_due(cs) && FieldExchanger.wait_slow_sims!(cs)
    return FieldExchanger.slow_step_in_flight(cs)
end

"""
    skip_slow_stepping(cs, frozen)

Whether `step_model_sims!` should leave the slow group alone.

Under overlap they are advanced only by the asynchronous task, never in the loop:
when the slow timestep equals the coupling timestep they would otherwise be
stepped both synchronously and again by the task.
"""
skip_slow_stepping(cs::Interfacer.CoupledSimulation, frozen::Bool) =
    frozen || cs.overlap_slow_surfaces

"""
    launch_slow_if_due!(cs, frozen) -> launched

Launch the next overlapped ice/ocean step if this coupling step is a boundary,
and advance the recorded boundary when priming.

Called after the slow group's forcing is fully assembled, so the step it starts
integrates with complete inputs.
"""
function launch_slow_if_due!(cs::Interfacer.CoupledSimulation, frozen::Bool)
    (cs.overlap_slow_surfaces && !frozen && slow_step_due(cs)) || return false
    FieldExchanger.launch_slow_sims!(cs)
    if cs.prime_slow_surfaces
        cs.slow_next_boundary[] =
            cs.slow_next_boundary[] + FieldExchanger.slow_window_steps(cs) * cs.Δt_cpl
    end
    return true
end

"""
    at_slow_boundary(cs)

Whether the coupler has just reached a slow-step boundary.

Used only when priming. The window is counted in coupling steps rather than
derived from a component clock, because priming deliberately moves those clocks
away from coupler time. Integer step counting also sidesteps comparing times as
floats, which is the reason `ITime` exists.
"""
function at_slow_boundary(cs::Interfacer.CoupledSimulation)
    nb = cs.slow_next_boundary[]
    isnothing(nb) && return false
    return Float64(float(cs.t[])) >= Float64(float(nb))
end

"""
    slow_surfaces_due(cs)

Whether the overlapped ice/ocean group would take a step at the next coupling
time. This is the same predicate `push_ready_accumulators!` uses to decide when
to deliver their window-averaged forcing, so launching on it guarantees the
forcing is complete before the step begins.
"""
function slow_surfaces_due(cs::Interfacer.CoupledSimulation)
    t_next = cs.t[] + cs.Δt_cpl
    target = FieldExchanger.slow_step_target(cs)
    if !isnothing(target)
        # A step is in flight, so the ocean clock is being written and must not
        # be read. Once the step lands the clock will read `target`, so the next
        # step falls due one slow timestep after that.
        return Float64(float(t_next)) >= target + FieldExchanger.slow_sim_dt(cs)
    end
    for sim in cs.model_sims
        Interfacer.is_overlapped(sim) || continue
        Interfacer.will_step(sim, t_next) && return true
    end
    return false
end

"""
    save_sypd_walltime_to_disk(cs, walltime, t_timed_start = cs.tspan[begin])

Save the computed `sypd`, `walltime_per_coupling_step`, and memory usage to text files 
in the `artifacts` directory. `t_timed_start` is the simulated time (in seconds) at which 
the timed portion of the run began (see the `precompile` flag in [run!](@ref)). 
"""
function save_sypd_walltime_to_disk(cs, walltime, t_timed_start = cs.tspan[begin])
    if ClimaComms.iamroot(ClimaComms.context(cs))
        sypd = TimeManager.simulated_years_per_day(t_timed_start, cs.tspan[end], walltime)
        walltime_per_step = walltime / ((cs.tspan[end] - t_timed_start) / cs.Δt_cpl)
        @info "SYPD: $sypd"
        @info "Walltime per coupling step: $(walltime_per_step)"

        open(joinpath(cs.dir_paths.artifacts_dir, "sypd.txt"), "w") do sypd_filename
            println(sypd_filename, "$sypd")
        end

        open(
            joinpath(cs.dir_paths.artifacts_dir, "walltime_per_step.txt"),
            "w",
        ) do walltime_per_step_filename
            println(walltime_per_step_filename, "$(walltime_per_step)")
        end
    end
    return nothing
end

"""
    CoupledSimulation(config_file = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml"))
    CoupledSimulation(config_dict)

Set up a `CoupledSimulation` as prescribed by the given input.

This struct is defined in the Interfacer module and contains all information
about component models, diagnostics, timestepping, output directories, etc
needed to run a coupled simulation.

If no arguments are provided, the default AMIP configuration is used,
which is defined in `config/ci_configs/amip_default.yml`.
"""
function Interfacer.CoupledSimulation(
    config_file::AbstractString = joinpath(
        pkgdir(ClimaCoupler),
        "config/ci_configs/amip_default.yml",
    ),
)
    config_dict = Input.get_coupler_config_dict(config_file)
    return Interfacer.CoupledSimulation(config_dict)
end

function Interfacer.CoupledSimulation(config_dict::AbstractDict)
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
        step_concurrently,
        overlap_slow_surfaces,
        prime_slow_surfaces,
        share_surface_space,
        nh_poly_coupler,
        h_elem_coupler,
        saveat,
        checkpoint_dt,
        walltime_dt,
        walltime_debug,
        atmos_progress_interval,
        detect_restart_files,
        restart_dir,
        restart_t,
        restart_cache,
        save_cache,
        use_land_diagnostics,
        land_diagnostics_period,
        land_diagnostics_reduction,
        land_progress_interval,
        evolving_ocean,
        land_model,
        land_spun_up_ic,
        lai_source,
        bucket_albedo_type,
        energy_check,
        use_coupler_diagnostics,
        coupler_diagnostics_period,
        coupler_diagnostics_reduction,
        output_dir_root,
        parameter_files,
        era5_filepaths,
        ocean_model,
        simple_ocean,
        ocean_grid,
        use_intersection_grid,
        sst_adjustment,
        ocean_progress_interval,
        ocean_diagnostic_interval,
        ocean_diagnostic_mode,
        ice_model,
        seaice_diagnostic_interval,
        seaice_diagnostic_mode,
        seaice_progress_interval,
        land_fraction_source,
        binary_area_fraction,
        domain_type,
        column_latlon,
        scm_surface_type,
    ) = Input.get_coupler_args(config_dict)

    override_file = CP.merge_toml_files(parameter_files; override = true)
    coupled_param_dict = CP.create_toml_dict(FT; override_file)
    thermo_params = TDP.ThermodynamicsParameters(coupled_param_dict)

    dir_paths = Utilities.setup_output_dirs(
        output_dir_root = output_dir_root,
        comms_ctx = comms_ctx,
    )

    Random.seed!(random_seed)
    @info "Random seed set to $(random_seed)"

    # Concurrent component stepping only pays off on a GPU, where each component
    # occupies a single Julia thread and submits to its own CUDA stream. On a CPU
    # device every component fans out over all threads through KernelAbstractions,
    # so the components oversubscribe each other and nothing is gained. Neither
    # case is incorrect, so warn rather than error.
    if step_concurrently
        if !(comms_ctx.device isa ClimaComms.CUDADevice)
            @warn "`step_concurrently` is set, but the device is \
                   $(nameof(typeof(comms_ctx.device))), not CUDADevice. Component models \
                   will be stepped in separate tasks that contend for the same threads, \
                   which adds overhead without speeding anything up."
        elseif Threads.nthreads() == 1
            @warn "`step_concurrently` is set, but Julia is running with a single thread, \
                   so the component tasks will time-share it and run effectively \
                   sequentially. Start Julia with `--threads=N` (N ≥ 2) to get concurrency."
        end
    end

    tspan = (t_start, t_end)
    @info "Starting from t_start $(t_start)"

    #=
    ## Component Model Initialization
    Each component model is required to have an `init` function that
    returns a `AbstractComponentSimulation` object (see `Interfacer` docs for more details).
    =#

    atmos_sim = Interfacer.AtmosSimulation(
        Val(:climaatmos);
        config_dict,
        atmos_output_dir = dir_paths.atmos_output_dir,
        coupled_param_dict,
        comms_ctx,
    )

    #=
    ### Boundary Space
    We use a boundary space at the surface for coupling operations (computing fluxes, regridding, etc).
    For column mode, this is a 1D PointSpace with lat/long coordinates.
    For global mode, this is a 2D CubedSphereSpace or the atmosphere's horizontal space
    (if `share_surface_space` is true).
    =#
    boundary_space = Utilities.create_boundary_space(
        FT,
        domain_type,
        atmos_sim,
        share_surface_space,
        comms_ctx;
        column_latlon,
        nh_poly_coupler,
        h_elem_coupler,
        coupled_param_dict,
    )

    surface_elevation = Interfacer.get_field(boundary_space, atmos_sim, Val(:height_sfc))
    atmos_bottom_center_height =
        Interfacer.get_field(boundary_space, atmos_sim, Val(:height_int))
    atmos_h =
        Interfacer.get_atmos_height_delta(atmos_bottom_center_height, surface_elevation)
    initial_T = CC.Fields.zeros(boundary_space)
    initial_T .= Interfacer.get_field(boundary_space, atmos_sim, Val(:air_temperature))

    land_fraction = Input.get_land_fraction(
        boundary_space,
        comms_ctx;
        land_fraction_source,
        binary_area_fraction,
        sim_mode,
        domain_type,
        scm_surface_type,
    )

    #=
    ### Surface Models
    Initialize land, ocean, and sea ice component models.
    =#
    @info(sim_mode)
    land_sim = ice_sim = ocean_sim = nothing

    (; sst_path, sic_path, land_ic_path, albedo_path, bucket_initial_condition) =
        era5_filepaths

    shared_surface_space =
        (share_surface_space || domain_type == "column") ? boundary_space : nothing
    land_sim = Interfacer.LandSimulation(
        FT,
        land_model;
        dt = component_dt_dict["dt_land"],
        tspan,
        start_date,
        output_dir = dir_paths.land_output_dir,
        area_fraction = land_fraction,
        shared_surface_space,
        atmos_h,
        initial_T,
        use_land_diagnostics,
        land_diagnostics_period,
        land_diagnostics_reduction,
        coupled_param_dict,
        albedo_type = bucket_albedo_type,
        bucket_initial_condition,
        era5_albedo_file_path = albedo_path,
        land_spun_up_ic,
        land_ic_path,
        lai_source,
        dt_drivers = ITime(Utilities.time_to_seconds(config_dict["dt_rad"])),
    )

    ocean_sim = Interfacer.OceanSimulation(
        FT,
        ocean_model;
        dt = component_dt_dict["dt_ocean"],
        start_date,
        tspan,
        coupled_param_dict,
        thermo_params,
        comms_ctx,
        boundary_space,
        output_dir = dir_paths.ocean_output_dir,
        simple_ocean,
        ocean_grid,
        use_intersection_grid,
        sst_path,
        sst_adjustment,
        saveat,
        evolving = evolving_ocean,
        ocean_diagnostic_interval,
        ocean_diagnostic_mode,
    )

    ice_sim = Interfacer.SeaIceSimulation(
        FT,
        ice_model;
        dt = component_dt_dict["dt_seaice"],
        start_date,
        coupled_param_dict,
        output_dir = dir_paths.ice_output_dir,
        ocean = ocean_sim,
        tspan,
        saveat,
        boundary_space,
        thermo_params,
        comms_ctx,
        land_fraction,
        sic_path,
        binary_area_fraction,
        domain_type,
        seaice_diagnostic_interval,
        seaice_diagnostic_mode,
    )

    #=
    ## Coupler Initialization
    =#
    model_sims = (; atmos_sim, ice_sim, land_sim, ocean_sim)
    model_sims =
        NamedTuple{filter(key -> !isnothing(model_sims[key]), keys(model_sims))}(model_sims)
    @info "Component models initialized: $(keys(model_sims))"
    @info "Component model types: $(nameof.(values(model_sims)))"

    coupler_field_names = Interfacer.default_coupler_fields()
    foreach(sim -> Interfacer.add_coupler_fields!(coupler_field_names, sim), model_sims)

    energy_check && push!(coupler_field_names, :P_net)
    overlap_slow_surfaces && append!(coupler_field_names, Interfacer.overlap_cache_fields())

    coupler_fields = Interfacer.init_coupler_fields(FT, coupler_field_names, boundary_space)

    # Allocate a FluxAccumulator per slow explicit surface (one whose timestep is
    # strictly greater than Δt_cpl). Other surfaces avoid this and receive turbulent
    # fluxes directly via `update_turbulent_fluxes!` each coupling step.
    Δt_cpl_secs = Float64(float(Δt_cpl))
    slow_surface_keys = Tuple(
        name for (name, sim) in pairs(model_sims) if
        sim isa Interfacer.AbstractSurfaceSimulation &&
            !(sim isa Interfacer.AbstractImplicitFluxSimulation) &&
            !(sim isa Interfacer.AbstractSurfaceStub) &&
            Interfacer.sim_dt(sim) > Δt_cpl_secs
    )
    flux_accumulators = NamedTuple{slow_surface_keys}(
        Tuple(FluxCalculator.FluxAccumulator(boundary_space) for _ in slow_surface_keys),
    )
    isempty(slow_surface_keys) ||
        @info "Allocated flux accumulators for slow surfaces: $(slow_surface_keys)"

    # set initial area fractions (remain set to 0 if model does not exist)
    if haskey(model_sims, :land_sim)
        coupler_fields.land_area_fraction .=
            Interfacer.get_field(model_sims.land_sim, Val(:area_fraction))
    end
    if haskey(model_sims, :ice_sim)
        coupler_fields.ice_area_fraction .=
            Interfacer.get_field(model_sims.ice_sim, Val(:area_fraction))
    end
    if haskey(model_sims, :ocean_sim)
        coupler_fields.ocean_area_fraction .=
            Interfacer.get_field(model_sims.ocean_sim, Val(:area_fraction))
    end

    ## Conservation checks (only applicable to global slabplanet mode)
    conservation_checks = nothing
    if energy_check && domain_type == "global"
        @assert(
            sim_mode <: Interfacer.AbstractSlabplanetSimulationMode &&
            comms_ctx isa ClimaComms.SingletonCommsContext,
            "Only non-distributed slabplanet allowable for energy_check"
        )
        conservation_checks = (;
            energy = ConservationChecker.EnergyConservationCheck(model_sims),
            water = ConservationChecker.WaterConservationCheck(model_sims),
        )
    elseif energy_check && domain_type == "column"
        @warn "Conservation checks are disabled for single-column mode."
    end

    ## Callbacks
    # TODO: Move callbacks code somewhere else (maybe in TimeManager?) so that it doesn't clutter up the constructor

    # checkpoint
    # Schedules are seeded with t_start so that a restarted simulation stays on
    # the same calendar boundaries as the original run (see calendar_dt_schedule).
    schedule_checkpoint =
        TimeManager.calendar_dt_schedule(checkpoint_dt, start_date, t_start)
    checkpoint_cb =
        TimeManager.Callback(schedule_checkpoint, sim -> Checkpointer.checkpoint_sims(sim))

    # walltime reporting
    schedule_walltime =
        TimeManager.walltime_schedule(walltime_dt, walltime_debug, start_date, t_start)
    if isnothing(schedule_walltime)
        callbacks = (checkpoint_cb,)
    else
        walltime_cb =
            TimeManager.Callback(schedule_walltime, TimeManager.WalltimeReporter())
        callbacks = (checkpoint_cb, walltime_cb)
    end

    # component model progress reporting
    progress_intervals = (;
        atmos_sim = atmos_progress_interval,
        ocean_sim = ocean_progress_interval,
        land_sim = land_progress_interval,
        ice_sim = seaice_progress_interval,
    )
    for (sim_name, interval) in pairs(progress_intervals)
        (haskey(model_sims, sim_name) && interval != "never") || continue
        schedule_progress = TimeManager.calendar_dt_schedule(interval, start_date, t_start)
        progress_cb = TimeManager.Callback(
            schedule_progress,
            let sim_name = sim_name
                cs -> begin
                    sim = cs.model_sims[sim_name]
                    # `progress` takes extrema and maxima over the model's own
                    # fields. For an overlapped sim those fields are being
                    # written by its in-flight step, so report from the snapshot
                    # the stepping task gathered when it last finished, rather
                    # than reducing over state that is moving underneath us.
                    if Interfacer.is_overlapped(sim) &&
                       FieldExchanger.slow_step_in_flight(cs)
                        snapshot = FieldExchanger.slow_progress_snapshot(cs, sim_name)
                        isnothing(snapshot) && return nothing
                        return Interfacer.progress(sim, cs, snapshot)
                    end
                    Interfacer.progress(sim, cs)
                end
            end,
        )
        callbacks = (callbacks..., progress_cb)
    end

    ## Coupler diagnostics
    if use_coupler_diagnostics
        @info "Using default coupler diagnostics"
        diags_handler = SimOutput.diagnostics_setup(
            coupler_fields,
            dir_paths.coupler_output_dir,
            start_date,
            tspan[1],
            coupler_diagnostics_period,
            Δt_cpl;
            reduction = coupler_diagnostics_reduction,
        )
    else
        diags_handler = nothing
    end

    ## Build the CoupledSimulation struct
    prev_checkpoint_t = Ref(-1)
    cs = Interfacer.CoupledSimulation{FT}(
        start_date,
        coupler_fields,
        conservation_checks,
        [tspan[1], tspan[2]],
        Δt_cpl,
        Ref(tspan[1]),
        Ref(0),
        prev_checkpoint_t,
        model_sims,
        callbacks,
        dir_paths,
        thermo_params,
        diags_handler,
        save_cache,
        step_concurrently,
        overlap_slow_surfaces,
        prime_slow_surfaces,
        Ref{Any}(nothing),
        Ref{Any}(nothing),
        Ref{Any}(nothing),
        flux_accumulators,
    )

    ## Restart component model states if specified
    if detect_restart_files
        isnothing(restart_t) &&
            (restart_t = Checkpointer.t_start_from_checkpoint(dir_paths.checkpoints_dir))
        isnothing(restart_dir) && (restart_dir = dir_paths.checkpoints_dir)
    end
    should_restart = !isnothing(restart_t) && !isnothing(restart_dir)
    should_restart && Checkpointer.restart!(cs, restart_dir, restart_t, restart_cache)

    FieldExchanger.update_surface_fractions!(cs)

    if !should_restart || !restart_cache
        ## Initialize Component Model Exchange
        FieldExchanger.import_static_fields!(cs.fields, cs.model_sims)
        FieldExchanger.exchange!(cs)
        FieldExchanger.set_caches!(cs)
        FluxCalculator.turbulent_fluxes!(cs)
        FluxCalculator.ocean_seaice_fluxes!(cs)

        # The initial `turbulent_fluxes!` above fills the flux accumulators.
        # Push the accumulated flux to each slow surface here.
        FluxCalculator.push_ready_accumulators!(
            cs.model_sims,
            cs.flux_accumulators,
            cs.t[];
            force = true,
        )

        # Priming: take the slow group's first step here, synchronously, so it
        # ends up one window ahead of the coupler. Thereafter each overlapped
        # step integrates the window that is about to happen rather than the one
        # that just did, and its result is ready when the atmosphere needs it --
        # removing the extra lag in the ocean state the atmosphere sees. The
        # forcing it uses is the previous window's, which the slow components
        # can absorb: their timestep is several coupling steps long precisely
        # because their physics is slow, so they already integrate under forcing
        # held constant across a whole window.
        if overlap_slow_surfaces && prime_slow_surfaces
            k = FieldExchanger.slow_window_steps(cs)
            FieldExchanger.step_slow_sims!(cs.model_sims, cs.t[] + k * cs.Δt_cpl)
            @info "Primed slow surfaces one step ($k coupling steps) ahead of the coupler"
        end
    end
    # Seed the overlapped-group schedule from where the slow components actually
    # are. This must come after both the restart restore and any priming step.
    #
    # It cannot be a count of coupling steps: `cs.step[]` starts again at zero on
    # a restart while the component clocks do not, and `checkpoint_sims` joins an
    # in-flight step before saving, so a checkpoint taken mid-window leaves the
    # slow group an arbitrary amount ahead rather than a whole window.
    if overlap_slow_surfaces && prime_slow_surfaces
        cs.slow_next_boundary[] = FieldExchanger.slow_step_boundary(cs)
        @info """Priming is on: the ocean and sea ice run ahead of coupler time.
                 Next overlapped launch at coupler time $(cs.slow_next_boundary[]).

                 Their own diagnostics are written by their own output writers on
                 their own clocks, so those files carry times that lead coupler
                 time by up to one slow step. That is the honest label: the state
                 in them really is the state at that model time, forced through
                 one window earlier. Coupler diagnostics are on coupler time and
                 hold the surface state the atmosphere actually saw, which priming
                 keeps aligned with a non-overlapped run.

                 Comparing a primed run against a non-primed one by output index
                 therefore compares different model times for the slow components;
                 compare by time, or expect an offset of one slow step."""
    end

    Utilities.show_memory_usage()
    return cs
end

"""
    setup_and_run(config_dict::AbstractDict)
    setup_and_run(config_file::AbstractString = joinpath(pkgdir(ClimaCoupler), "config/ci_configs/amip_default.yml"))

Set up and run the coupled model simulation specified by the input config
file or dict. Returns the `CoupledSimulation` after the run completes.

If no arguments are provided, the default AMIP configuration is used,
which is defined in `config/ci_configs/amip_default.yml`. This
is the same behavior as the `CoupledSimulation` constructor.
"""
function setup_and_run(
    config_file::AbstractString = joinpath(
        pkgdir(ClimaCoupler),
        "config/ci_configs/amip_default.yml",
    ),
)
    cs = Interfacer.CoupledSimulation(config_file)
    run!(cs)
    return cs
end

function setup_and_run(config_dict::AbstractDict)
    cs = Interfacer.CoupledSimulation(config_dict)
    run!(cs)
    return cs
end

end # module SimCoordinator
