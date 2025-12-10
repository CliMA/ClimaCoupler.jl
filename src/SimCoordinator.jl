"""
    SimCoordinator

This module handles timestepping and running of coupled simulations.
It coordinates the execution of component models and manages the simulation lifecycle.
"""
module SimCoordinator

import ClimaComms
import ..Interfacer
import ..ConservationChecker
import ..FieldExchanger
import ..FluxCalculator
import ..TimeManager
import ..Postprocessor

export run!, step!

"""
    run!(cs::Interfacer.CoupledSimulation; precompile = ...)

Run the coupled simulation for the full timespan.

Keyword arguments
==================

`precompile`: If `true`, run the coupled simulations for two steps, so that most functions
              are precompiled and subsequent timing will be more accurate.
"""
function run!(
    cs::Interfacer.CoupledSimulation;
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
    # Note: ClimaAtmos timing and helper functions are expected to be available in the calling context
    # This allows SimCoordinator to remain independent while still supporting ClimaEarth-specific functionality
    walltime = ClimaComms.@elapsed ClimaComms.device(cs) begin
        while cs.t[] < cs.tspan[end]
            step!(cs)
        end
    end
    @info "Simulation took $(walltime) seconds"

    # Compute and log performance metrics
    sypd = Postprocessor.simulated_years_per_day(cs, walltime)
    walltime_per_step = Postprocessor.walltime_per_coupling_step(cs, walltime)
    @info "SYPD: $sypd"
    @info "Walltime per coupling step: $(walltime_per_step)"
    Postprocessor.save_sypd_walltime_to_disk(cs, walltime)

    # Close all diagnostics file writers
    isnothing(cs.diags_handler) ||
        foreach(diag -> close(diag.output_writer), cs.diags_handler.scheduled_diagnostics)
    foreach(Interfacer.close_output_writers, cs.model_sims)

    return nothing
end

"""
    step!(cs::Interfacer.CoupledSimulation)

Take one coupling step forward in time.

This function runs the component models sequentially, and exchanges combined fields and
calculates fluxes using the selected turbulent fluxes option. Note, one coupling step might
require multiple steps in some of the component models.
"""
function step!(cs::Interfacer.CoupledSimulation)
    # Update the current time
    cs.t[] += cs.Δt_cpl

    ## compute global energy and water conservation checks
    ## (only for slabplanet if tracking conservation is enabled)
    ConservationChecker.check_conservation!(cs)

    ## step component model simulations sequentially for one coupling timestep (Δt_cpl)
    FieldExchanger.step_model_sims!(cs)

    ## update the surface fractions for surface models
    FieldExchanger.update_surface_fractions!(cs)

    ## exchange all non-turbulent flux fields between models, including radiative and precipitation fluxes
    FieldExchanger.exchange!(cs)

    ## calculate turbulent fluxes in the coupler and update the model simulations with them
    FluxCalculator.turbulent_fluxes!(cs)

    ## compute any ocean-sea ice fluxes
    FluxCalculator.ocean_seaice_fluxes!(cs)

    ## Maybe call the callbacks
    TimeManager.callbacks!(cs)

    # Compute coupler diagnostics
    Postprocessor.CD.orchestrate_diagnostics(cs)
    return nothing
end

end # module
