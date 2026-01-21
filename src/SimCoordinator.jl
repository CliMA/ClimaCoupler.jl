"""
    SimCoordinator

This module contains functions for coordinating the execution of coupled simulations,
including stepping through time and running full simulations.
"""
module SimCoordinator

import ClimaComms
import ClimaDiagnostics as CD
import ClimaAtmos as CA

import ..Interfacer
import ..ConservationChecker
import ..FieldExchanger
import ..FluxCalculator
import ..TimeManager

export run!, step!

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

end # module SimCoordinator
