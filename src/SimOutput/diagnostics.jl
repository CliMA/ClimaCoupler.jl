import ClimaDiagnostics as CD
import ClimaCoupler: Interfacer
import Dates

export diagnostics_setup

#### Custom Schedule for diagnostics that only get output once at the beginning of the simulation

struct OnceSchedule <: CD.Schedules.AbstractSchedule end

function (::OnceSchedule)(integrator)
    return integrator.step == 1
end

function CD.Schedules.short_name(::OnceSchedule)
    return "once"
end

function CD.Schedules.long_name(::OnceSchedule)
    return "once at the first step of the simulation"
end

#### Diagnostics orchestration and setup functions

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
    diagnostics_setup(fields, output_dir, start_date, t_start, diagnostics_dt)

Set up the default diagnostics for an AMIP simulation, using ClimaDiagnostics.
The diagnostics are saved to NetCDF files. Currently, this just includes a
diagnostic for turbulent energy fluxes and diagnostics for each of the land, ocean,
and sea-ice area fractions.

Return a DiagnosticsHandler object to coordinate the diagnostics.
"""
function diagnostics_setup(
    fields,
    output_dir,
    start_date,
    t_start,
    diagnostics_dt,
    coupled_dt,
)
    # Create a list to hold the scheduled diagnostics
    scheduled_diags = []

    # Create output writer (shared across all diagnostics since they all live on the boundary space)
    boundary_space = axes(fields.F_lh)
    netcdf_writer = CD.Writers.NetCDFWriter(boundary_space, output_dir)

    # Create schedules for computing and outputting diagnostics
    schedule_once = OnceSchedule()
    schedule_everystep = CD.Schedules.EveryStepSchedule()
    schedule_calendar_dt = CD.Schedules.EveryCalendarDtSchedule(diagnostics_dt; start_date)

    #### Turbulent energy fluxes diagnostic

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

    # Schedule the turbulent energy fluxes to save at every step and output at the frequency calculated above
    F_turb_energy_diag_sched = CD.ScheduledDiagnostic(
        variable = F_turb_energy_diag,
        output_writer = netcdf_writer,
        reduction_time_func = (+),
        compute_schedule_func = schedule_everystep,
        output_schedule_func = schedule_calendar_dt,
        pre_output_hook! = CD.average_pre_output_hook!,
    )

    push!(scheduled_diags, F_turb_energy_diag_sched)

    #### Land area fraction diagnostic (only at beginning of the simulation since it's static)

    # Create the diagnostic for land fraction
    land_fraction_diag = CD.DiagnosticVariable(;
        short_name = "land_fraction",
        long_name = "Land area fraction",
        standard_name = "land_fraction",
        units = "1",
        comments = "Fraction of each grid cell that is land",
        compute! = (out, state, cache, time) -> begin
            if isnothing(out)
                return state.land_area_fraction
            else
                out .= state.land_area_fraction
            end
        end,
    )

    # Schedule the land fraction to save and output at the first step only
    land_fraction_diag_sched = CD.ScheduledDiagnostic(
        variable = land_fraction_diag,
        output_writer = netcdf_writer,
        compute_schedule_func = schedule_once, # only compute at the first step
        output_schedule_func = schedule_once,  # only output at the first step
    )

    push!(scheduled_diags, land_fraction_diag_sched)

    #### Ocean area fraction diagnostic

    # Create the diagnostic for ocean fraction
    ocean_fraction_diag = CD.DiagnosticVariable(;
        short_name = "ocean_fraction",
        long_name = "Ocean area fraction",
        standard_name = "ocean_fraction",
        units = "1",
        comments = "Fraction of each grid cell that is ocean",
        compute! = (out, state, cache, time) -> begin
            if isnothing(out)
                return state.ocean_area_fraction
            else
                out .= state.ocean_area_fraction
            end
        end,
    )

    # Schedule the ocean fraction to save and output at diagnostic frequency 
    # since it can change in time with evolving sea ice
    ocean_fraction_diag_sched = CD.ScheduledDiagnostic(
        variable = ocean_fraction_diag,
        output_writer = netcdf_writer,
        compute_schedule_func = schedule_calendar_dt,
        output_schedule_func = schedule_calendar_dt,
    )

    push!(scheduled_diags, ocean_fraction_diag_sched)

    #### Ice area fraction diagnostic

    # Create the diagnostic for ice fraction
    ice_fraction_diag = CD.DiagnosticVariable(;
        short_name = "ice_fraction",
        long_name = "Ice area fraction",
        standard_name = "ice_fraction",
        units = "1",
        comments = "Fraction of each grid cell that is ice",
        compute! = (out, state, cache, time) -> begin
            if isnothing(out)
                return state.ice_area_fraction
            else
                out .= state.ice_area_fraction
            end
        end,
    )

    # Schedule the ice fraction to save and output at diagnostic frequency 
    ice_fraction_diag_sched = CD.ScheduledDiagnostic(
        variable = ice_fraction_diag,
        output_writer = netcdf_writer,
        compute_schedule_func = schedule_calendar_dt,
        output_schedule_func = schedule_calendar_dt,
    )

    push!(scheduled_diags, ice_fraction_diag_sched)

    # Create the diagnostics handler containing the scheduled diagnostics
    diags_handler =
        CD.DiagnosticsHandler(scheduled_diags, fields, nothing, t_start, dt = coupled_dt)

    return diags_handler
end
