import Dates

"""
    diagnostic_schedule(mode, interval)

Return an Oceananigans output schedule corresponding to `mode`:
- `:average` → `OC.AveragedTimeInterval(interval)`
- `:instantaneous` → `OC.TimeInterval(interval)`

Throws for any other value.
"""
function diagnostic_schedule(mode::Symbol, interval)
    if mode === :average
        return OC.AveragedTimeInterval(interval)
    elseif mode === :instantaneous
        return OC.TimeInterval(interval)
    else
        error("Unknown diagnostic mode `$(mode)`. Expected `:average` or `:instantaneous`.")
    end
end

"""
    add_ocean_diagnostics!(ocean_sim::OceananigansSimulation;
                           output_dir = joinpath("output_active", "clima_ocean"),
                           interval = Dates.Day(1),
                           mode = :average,
                           filename_prefix = "ocean",
                           file_splitting_interval = Dates.Day(15))

Attach output writers to the underlying Oceananigans simulation inside an `OceananigansSimulation`.
Two writers are added to `ocean_sim.ocean.output_writers`, both with the same `interval` and `mode`:

1. **Surface diagnostics** (`<prefix>_surface.jld2`): 2-D fields —
   SST, SSS, SSH, surface velocities, squared surface fields for variance, surface momentum/heat/freshwater fluxes.
2. **3-D field diagnostics** (`<prefix>_fields.jld2`): full 3-D temperature, salinity, velocity (and TKE if a `:e` tracer is present).

`mode` selects the reduction:
- `:average` uses `Oceananigans.AveragedTimeInterval(interval)` (time-averaged fields).
- `:instantaneous` uses `Oceananigans.TimeInterval(interval)` (snapshots).

!!! note
    `tauuo`, `tauvo`, `hfds`, and `wfo` here are the *combined* surface fluxes on the ocean (atmosphere + sea-ice contributions),
    not the atmosphere-only contribution. Sensible/latent heat (`hfss`/`hfls`) are not stored as standalone fields in CMIP —
    they are folded into `oc_flux_T` inside `update_turbulent_fluxes!` — and are therefore omitted.
"""
function add_ocean_diagnostics!(
    ocean_sim::OceananigansSimulation;
    output_dir = joinpath("output_active", "clima_ocean"),
    interval = Dates.Day(1),
    mode = :average,
    filename_prefix = "ocean",
    file_splitting_interval = Dates.Day(15),
)
    ocean = ocean_sim.ocean
    Nz = size(ocean.model.grid, 3)
    file_splitting = OC.TimeInterval(file_splitting_interval)

    T, S = ocean.model.tracers.T, ocean.model.tracers.S
    u, v, w = ocean.model.velocities
    η = ocean.model.free_surface.displacement

    τx = surface_flux(u)
    τy = surface_flux(v)
    JT = surface_flux(T)
    Js = surface_flux(S)

    # Seed FieldStatus with the clock time so computed fields adopt the clock's time type
    computed_field_status = OC.Fields.FieldStatus(ocean.model.clock.time)

    surface_indices = (:, :, Nz)
    sea_surface_temperature = OC.Field(T; indices = surface_indices)
    sea_surface_salinity = OC.Field(S; indices = surface_indices)
    sea_surface_zonal_velocity = OC.Field(u; indices = surface_indices)
    sea_surface_meridional_velocity = OC.Field(v; indices = surface_indices)
    sea_surface_temperature_squared =
        OC.Field(T * T; indices = surface_indices, status = computed_field_status)
    sea_surface_salinity_squared =
        OC.Field(S * S; indices = surface_indices, status = computed_field_status)
    sea_surface_height_squared = OC.Field(η * η; status = computed_field_status)

    surface_outputs = Dict{Symbol, Any}(
        :sea_surface_temperature => sea_surface_temperature,
        :sea_surface_salinity => sea_surface_salinity,
        :sea_surface_height => η,
        :sea_surface_zonal_velocity => sea_surface_zonal_velocity,
        :sea_surface_meridional_velocity => sea_surface_meridional_velocity,
        :sea_surface_temperature_squared => sea_surface_temperature_squared,
        :sea_surface_salinity_squared => sea_surface_salinity_squared,
        :sea_surface_height_squared => sea_surface_height_squared,
        :zonal_wind_stress => τx,
        :meridional_wind_stress => τy,
        :surface_heat_flux => JT,
        :surface_salinity_flux => Js,
    )

    schedule = diagnostic_schedule(mode, interval)

    ocean.output_writers[:surface_diagnostics] = OC.JLD2Writer(
        ocean.model,
        surface_outputs;
        schedule,
        filename = joinpath(output_dir, filename_prefix * "_surface"),
        file_splitting,
        overwrite_existing = true,
    )

    field_outputs = Dict{Symbol, Any}(
        :temperature => T,
        :salinity => S,
        :zonal_velocity => u,
        :meridional_velocity => v,
        :vertical_velocity => w,
    )

    if haskey(ocean.model.tracers, :e)
        field_outputs[:turbulent_kinetic_energy] = ocean.model.tracers.e
    end

    ocean.output_writers[:field_diagnostics] = OC.JLD2Writer(
        ocean.model,
        field_outputs;
        schedule,
        filename = joinpath(output_dir, filename_prefix * "_fields"),
        file_splitting,
        overwrite_existing = true,
    )

    @info "Ocean diagnostics attached:" *
          " surface ($(length(surface_outputs)) fields, $mode, every $interval)," *
          " 3-D ($(length(field_outputs)) fields, $mode, every $interval)"

    return nothing
end
