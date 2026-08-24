import Dates

"""
    add_seaice_diagnostics!(ice_sim::ClimaSeaIceSimulation;
                            output_dir = joinpath("output_active", "clima_seaice"),
                            interval = Dates.Day(1),
                            mode = :average,
                            filename_prefix = "seaice",
                            file_splitting_interval = Dates.Day(15))

Attach an output writer to the underlying ClimaSeaIce simulation inside a `ClimaSeaIceSimulation`.
A single writer is added to `ice_sim.ice.output_writers`:

1. **Surface diagnostics** (`<prefix>_surface.jld2`): sea-ice concentration, thickness, velocities, and top
   surface temperature, sampled every `interval`.

`mode` selects the reduction:
- `:average` uses `Oceananigans.AveragedTimeInterval(interval)` (time-averaged fields).
- `:instantaneous` uses `Oceananigans.TimeInterval(interval)` (snapshots).

`interval` and `file_splitting_interval` accept either a `Dates.Period` or a `Real` number of seconds.
The value is automatically converted to match the model clock's time type via [`to_clock_interval`](@ref):
a `Dates.Period` is converted to seconds for Float64 clocks, or kept as-is for DateTime clocks.
"""
function add_seaice_diagnostics!(
    ice_sim::ClimaSeaIceSimulation;
    output_dir = joinpath("output_active", "clima_seaice"),
    interval = Dates.Day(1),
    mode = :average,
    filename_prefix = "seaice",
    file_splitting_interval = Dates.Day(15),
)
    ice = ice_sim.ice
    clock_time = ice.model.clock.time
    resolved_interval = to_clock_interval(interval, clock_time)
    file_splitting = OC.TimeInterval(to_clock_interval(file_splitting_interval, clock_time))

    hi = ice.model.ice_thickness
    ℵi = ice.model.ice_concentration
    ui, vi = ice.model.velocities

    sitemptop = ice.model.ice_thermodynamics.top_surface_temperature

    surface_outputs = Dict{Symbol, Any}(
        :ice_concentration => ℵi,
        :ice_thickness => hi,
        :ice_zonal_velocity => ui,
        :ice_meridional_velocity => vi,
        :ice_top_temperature => sitemptop,
    )

    ice.output_writers[:surface_diagnostics] = OC.JLD2Writer(
        ice.model,
        surface_outputs;
        schedule = diagnostic_schedule(mode, resolved_interval),
        filename = joinpath(output_dir, filename_prefix * "_surface"),
        file_splitting,
        overwrite_existing = true,
    )

    @info "Sea-ice diagnostics attached:" *
          " surface ($(length(surface_outputs)) fields, $mode, every $interval)"

    return nothing
end
