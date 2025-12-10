"""
    Postprocessor

This module handles postprocessing of coupled simulations, including generating
plots, checking conservation, and other diagnostics.
"""
module Postprocessor

import ClimaComms
import ClimaDiagnostics as CD
import ClimaAnalysis as CAN
import ..Interfacer
import ..ConservationChecker
import Dates

export postprocess, postprocess_sim, simulated_years_per_day, walltime_per_coupling_step, save_sypd_walltime_to_disk
export make_diagnostics_plots, make_ocean_diagnostics_plots, debug, plot_global_conservation
export compute_leaderboard, compute_pfull_leaderboard, test_rmse_thresholds
export coupler_diagnostics_setup

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
function postprocess(cs::Interfacer.CoupledSimulation; conservation_softfail = false, rmse_check = false)
    if ClimaComms.iamroot(ClimaComms.context(cs)) && !isnothing(cs.diags_handler)
        postprocessing_vars = (; conservation_softfail, rmse_check)
        postprocess_sim(cs, postprocessing_vars)
    end
    return nothing
end

"""
    postprocess_sim(cs, postprocessing_vars)

Perform all postprocessing operations. This includes plotting all available
diagnostics, plotting all model states and coupler fields for debugging,
producing the leaderboard if monthly data is available, performing
conservation checks if enabled, and closing all diagnostics file writers.

Note: This function depends on ClimaEarth-specific helper functions that are
expected to be available in the calling context (e.g., from ClimaEarth experiments).
"""
function postprocess_sim(cs::Interfacer.CoupledSimulation, postprocessing_vars)
    (; conservation_softfail, rmse_check) = postprocessing_vars
    (;
        coupler_output_dir,
        atmos_output_dir,
        land_output_dir,
        ocean_output_dir,
        artifacts_dir,
    ) = cs.dir_paths

    # Plot generic diagnostics (ClimaEarth-specific functions)
    @info "Plotting diagnostics for coupler, atmos, land, and ocean"
    make_diagnostics_plots(coupler_output_dir, artifacts_dir, output_prefix = "coupler_")
    make_diagnostics_plots(atmos_output_dir, artifacts_dir, output_prefix = "atmos_")
    make_diagnostics_plots(land_output_dir, artifacts_dir, output_prefix = "land_")
    make_ocean_diagnostics_plots(ocean_output_dir, artifacts_dir, output_prefix = "ocean_")

    # Plot all model states and coupler fields (useful for debugging)
    if ClimaComms.context(cs) isa ClimaComms.SingletonCommsContext
        debug(cs, artifacts_dir)
    end

    # If we have enough data (in time, but also enough variables), plot the leaderboard.
    # We need pressure to compute the leaderboard.
    simdir = CAN.SimDir(atmos_output_dir)
    if !isempty(simdir)
        pressure_in_output = "pfull" in CAN.available_vars(simdir)
        first_var = first(CAN.available_vars(simdir))
        times = CAN.times(CAN.get(simdir; short_name = first_var))
        t_end = times[end]
        if t_end > 84600 * 31 * 3 # 3 months for spin up
            leaderboard_base_path = artifacts_dir
            compute_leaderboard(leaderboard_base_path, atmos_output_dir, 3)
            rmse_check && test_rmse_thresholds(atmos_output_dir, 3)
            pressure_in_output &&
                compute_pfull_leaderboard(leaderboard_base_path, atmos_output_dir, 3)
        end
    end

    # Perform conservation checks if they exist
    if !isnothing(cs.conservation_checks)
        @info "Conservation Check Plots"
        plot_global_conservation(
            cs.conservation_checks.energy,
            cs,
            conservation_softfail,
            figname1 = joinpath(artifacts_dir, "total_energy_bucket.png"),
            figname2 = joinpath(artifacts_dir, "total_energy_log_bucket.png"),
        )
        plot_global_conservation(
            cs.conservation_checks.water,
            cs,
            conservation_softfail,
            figname1 = joinpath(artifacts_dir, "total_water_bucket.png"),
            figname2 = joinpath(artifacts_dir, "total_water_log_bucket.png"),
        )
    end

    # Close all diagnostics file writers
    !isnothing(cs.diags_handler) &&
        map(diag -> close(diag.output_writer), cs.diags_handler.scheduled_diagnostics)
    return nothing
end

"""
    simulated_years_per_day(cs, walltime)

Compute the simulated years per walltime day for the given coupled simulation `cs`, assuming
that the simulation took `walltime`.
"""
function simulated_years_per_day(cs::Interfacer.CoupledSimulation, walltime)
    simulated_seconds_per_second = float(cs.tspan[end] - cs.tspan[begin]) / walltime
    return simulated_seconds_per_second / 365.25
end

"""
    walltime_per_coupling_step(cs, walltime)

Compute the average walltime needed to take one step for the given coupled simulation `cs`,
assuming that the simulation took `walltime`. The result is in seconds.
"""
function walltime_per_coupling_step(cs::Interfacer.CoupledSimulation, walltime)
    n_coupling_steps = (cs.tspan[end] - cs.tspan[begin]) / cs.Δt_cpl
    return walltime / n_coupling_steps
end

"""
    save_sypd_walltime_to_disk(cs, walltime)

Save the computed `sypd`, `walltime_per_coupling_step`,
and memory usage to text files in the `artifacts` directory.
"""
function save_sypd_walltime_to_disk(cs::Interfacer.CoupledSimulation, walltime)
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
    CD.orchestrate_diagnostics(cs::CoupledSimulation)

Compute and output coupled diagnostics.
"""
function CD.orchestrate_diagnostics(cs::Interfacer.CoupledSimulation)
    ## wrap the current CoupledSimulation fields and time in a NamedTuple to match the ClimaDiagnostics interface
    cs_nt = (; u = cs.fields, p = nothing, t = cs.t[], step = round(cs.t[] / cs.Δt_cpl))
    !isnothing(cs.diags_handler) && CD.orchestrate_diagnostics(cs_nt, cs.diags_handler)
    return nothing
end

"""
    coupler_diagnostics_setup(fields, output_dir, start_date, t_start, diagnostics_dt, coupled_dt)

Set up the default diagnostics for an AMIP simulation, using ClimaDiagnostics.
The diagnostics are saved to NetCDF files. Currently, this just includes a
diagnostic for turbulent energy fluxes.

Return a DiagnosticsHandler object to coordinate the diagnostics.
"""
function coupler_diagnostics_setup(
    fields,
    output_dir,
    start_date,
    t_start,
    diagnostics_dt,
    coupled_dt,
)
    # Create schedules and writer
    schedule_everystep = CD.Schedules.EveryStepSchedule()
    schedule_calendar_dt = CD.Schedules.EveryCalendarDtSchedule(diagnostics_dt; start_date)
    netcdf_writer = CD.Writers.NetCDFWriter(axes(fields.F_sh), output_dir)

    # Create the diagnostic for turbulent energy fluxes
    F_turb_energy_diag = CD.DiagnosticVariable(;
        short_name = "F_turb_energy",
        long_name = "Turbulent energy fluxes",
        standard_name = "F_turb_energy",
        units = "W m^-2",
        comments = "Turbulent energy fluxes are calculated as the sum of sensible and latent heat fluxes,
                    weighted by surface simulation area.",
        compute! = (out, state, cache, time) -> begin
            if isnothing(out)
                return state.F_sh .+ state.F_lh
            else
                out .= state.F_sh .+ state.F_lh
            end
        end,
    )

    # Schedule the turbulent energy fluxes to save at every step, and output at the frequency calculated above
    F_turb_energy_diag_sched = CD.ScheduledDiagnostic(
        variable = F_turb_energy_diag,
        output_writer = netcdf_writer,
        reduction_time_func = (+),
        compute_schedule_func = schedule_everystep,
        output_schedule_func = schedule_calendar_dt,
        pre_output_hook! = CD.average_pre_output_hook!,
    )

    # Create the diagnostics handler containing the scheduled diagnostics
    scheduled_diags = [F_turb_energy_diag_sched]
    diags_handler =
        CD.DiagnosticsHandler(scheduled_diags, fields, nothing, t_start, dt = coupled_dt)
    return diags_handler
end

include("diagnostics_plots.jl")
include("debug_plots.jl")
include("benchmarks.jl")
include("leaderboard/leaderboard.jl")
include("leaderboard/test_rmses.jl")

end # module
