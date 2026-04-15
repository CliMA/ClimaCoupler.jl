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
import ClimaDiagnostics.Schedules: EveryCalendarDtSchedule
import ClimaAtmos as CA
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

    #=
    ## Solving and Timing the Full Simulation

    This is where the full coupling loop, `solve_coupler!` is called for the full timespan of the simulation.
    We use the `ClimaComms.@elapsed` macro to time the simulation on both CPU and GPU, and use this
    value to calculate the simulated years per day (SYPD) of the simulation.
    =#
    @info "Starting coupling loop"
    walltime = ClimaComms.@elapsed ClimaComms.device(cs) begin
        s = CA.@timed_str begin
            while cs.t[] < cs.tspan[end]
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
    step!(cs::CoupledSimulation)

Take one coupling step forward in time.

This function runs the component models sequentially, and exchanges combined fields and
calculates fluxes using the selected turbulent fluxes option. Note, one coupling step might
require multiple steps in some of the component models.
"""
function step!(cs::Interfacer.CoupledSimulation)
    # Update the current time
    cs.t[] += cs.Δt_cpl

    # Compute global energy and water conservation checks
    # (only for slabplanet if tracking conservation is enabled)
    ConservationChecker.check_conservation!(cs)

    # Step component model simulations sequentially for one coupling timestep (Δt_cpl)
    FieldExchanger.step_model_sims!(cs)

    # Update the surface fractions for surface models
    FieldExchanger.update_surface_fractions!(cs)

    # Exchange all non-turbulent flux fields between models, including radiative and precipitation fluxes
    FieldExchanger.exchange!(cs)

    # Calculate turbulent fluxes in the coupler and update the model simulations with them
    FluxCalculator.turbulent_fluxes!(cs)

    # Compute any ocean-sea ice fluxes
    FluxCalculator.ocean_seaice_fluxes!(cs)

    # Maybe call the callbacks
    TimeManager.callbacks!(cs)

    # Compute and save coupler diagnostics
    CD.orchestrate_diagnostics(cs)
    return nothing
end

"""
    simulated_years_per_day(cs, walltime)

Compute the simulated years per walltime day for the given coupled simulation `cs`, assuming
that the simulation took `walltime`.
"""
function simulated_years_per_day(cs, walltime)
    simulated_seconds_per_second = float(cs.tspan[end] - cs.tspan[begin]) / walltime
    return simulated_seconds_per_second / 365.25
end

"""
    walltime_per_coupling_step(cs, walltime)

Compute the average walltime needed to take one step for the given coupled simulation `cs`,
assuming that the simulation took `walltime`. The result is in seconds.
"""
function walltime_per_coupling_step(cs, walltime)
    n_coupling_steps = (cs.tspan[end] - cs.tspan[begin]) / cs.Δt_cpl
    return walltime / n_coupling_steps
end

"""
    save_sypd_walltime_to_disk(cs, walltime)

Save the computed `sypd`, `walltime_per_coupling_step`,
and memory usage to text files in the `artifacts` directory.
"""
function save_sypd_walltime_to_disk(cs, walltime)
    if ClimaComms.iamroot(ClimaComms.context(cs))
        sypd = simulated_years_per_day(cs, walltime)
        walltime_per_step = walltime_per_coupling_step(cs, walltime)

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
        nh_poly,
        h_elem,
        coupled_param_dict,
    )

    surface_elevation = Interfacer.get_field(boundary_space, atmos_sim, Val(:height_sfc))
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
        surface_elevation,
        atmos_h,
        land_temperature_anomaly,
        use_land_diagnostics,
        coupled_param_dict,
        albedo_type = bucket_albedo_type,
        bucket_initial_condition,
        era5_albedo_file_path = albedo_path,
        land_spun_up_ic,
        land_ic_path,
        lai_source,
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
        sst_path,
        sst_adjustment,
        saveat,
        evolving = evolving_ocean,
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
    )

    #=
    ## Coupler Initialization
    =#
    model_sims = (; atmos_sim, ice_sim, land_sim, ocean_sim)
    model_sims =
        NamedTuple{filter(key -> !isnothing(model_sims[key]), keys(model_sims))}(model_sims)
    @info "Component models initialized: $(keys(model_sims))"

    coupler_field_names = Interfacer.default_coupler_fields()
    foreach(sim -> Interfacer.add_coupler_fields!(coupler_field_names, sim), model_sims)

    energy_check && push!(coupler_field_names, :P_net)

    coupler_fields = Interfacer.init_coupler_fields(FT, coupler_field_names, boundary_space)

    # set initial area fractions
    coupler_fields.land_area_fraction .= land_fraction
    coupler_fields.ice_area_fraction .= 0  # initialized as 0 since we start with no sea ice, but will evolve in time
    coupler_fields.ocean_area_fraction .= 1 .- land_fraction  # no sea ice
    @warn extrema(coupler_fields.land_area_fraction)
    @warn extrema(coupler_fields.ice_area_fraction)
    @warn extrema(coupler_fields.ocean_area_fraction)

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
    schedule_checkpoint =
        EveryCalendarDtSchedule(TimeManager.time_to_period(checkpoint_dt); start_date)
    checkpoint_cb =
        TimeManager.Callback(schedule_checkpoint, sim -> Checkpointer.checkpoint_sims(sim))

    if config_dict["atmos_log_progress"]
        callbacks = (checkpoint_cb,)
    else
        walltime_cb = TimeManager.capped_geometric_walltime_cb(t_start, t_end, Δt_cpl)
        callbacks = (checkpoint_cb, walltime_cb)
    end

    ## Coupler diagnostics
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

    ## Build the CoupledSimulation struct
    prev_checkpoint_t = Ref(-1)
    cs = Interfacer.CoupledSimulation{FT}(
        start_date,
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
