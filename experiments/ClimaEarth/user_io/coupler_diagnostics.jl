import ClimaDiagnostics as CD
import ClimaCoupler: Interfacer
import Dates

"""
    coupler_diagnostics_setup(fields, output_dir, start_date, t_start, diagnostics_dt)

Set up the default diagnostics for an AMIP simulation, using ClimaDiagnostics.
The diagnostics are saved to NetCDF files. Currently, this just includes a
diagnostic for turbulent energy fluxes.

Return a DiagnosticsHandler object to coordinate the diagnostics.
"""
function coupler_diagnostics_setup(fields, output_dir, start_date, t_start, diagnostics_dt, coupled_dt)
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
    diags_handler = CD.DiagnosticsHandler(scheduled_diags, fields, nothing, t_start, dt = coupled_dt)
    return diags_handler
end
