import Dates

"""
    add_ocean_diagnostics!(ocean_sim::OceananigansSimulation;
                           output_dir = ".",
                           surface_averaging_interval = Dates.Day(1),
                           field_averaging_interval = Dates.Day(15),
                           checkpoint_interval = Dates.Day(90),
                           filename_prefix = "ocean",
                           file_splitting_interval = Dates.Day(15))

Attach averaged-output writers to the underlying Oceananigans simulation inside an `OceananigansSimulation`. 
Three writers are added to `ocean_sim.ocean.output_writers`:

1. **Surface diagnostics** (`<prefix>_surface.jld2`): 2-D fields averaged over `surface_averaging_interval` —
   SST, SSS, SSH, surface velocities, squared surface fields for variance, surface momentum/heat/freshwater fluxes.
2. **3-D field diagnostics** (`<prefix>_fields.jld2`): full 3-D temperature, salinity, velocity (and TKE if a `:e` tracer is present), 
   averaged over `field_averaging_interval`.
3. **Checkpointer** (`<prefix>_checkpoint_iteration<N>.jld2`): JLD2 checkpoint of the ocean model at `checkpoint_interval`.

!!! note
    `tauuo`, `tauvo`, `hfds`, and `wfo` here are the *combined* surface fluxes on the ocean (atmosphere + sea-ice contributions), 
    not the atmosphere-only contribution. Sensible/latent heat (`hfss`/`hfls`) are not stored as standalone fields in CMIP — 
    they are folded into `oc_flux_T` inside `update_turbulent_fluxes!` — and are therefore omitted.
"""
function add_ocean_diagnostics!(
    ocean_sim::OceananigansSimulation;
    output_dir = ".",
    surface_averaging_interval = Dates.Day(1),
    field_averaging_interval = Dates.Day(15),
    checkpoint_interval = Dates.Day(90),
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

    ocean.output_writers[:surface_averages] = OC.JLD2Writer(
        ocean.model,
        surface_outputs;
        schedule = OC.AveragedTimeInterval(surface_averaging_interval),
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

    ocean.output_writers[:field_averages] = OC.JLD2Writer(
        ocean.model,
        field_outputs;
        schedule = OC.AveragedTimeInterval(field_averaging_interval),
        filename = joinpath(output_dir, filename_prefix * "_fields"),
        file_splitting,
        overwrite_existing = true,
    )

    ocean.output_writers[:checkpointer] = OC.Checkpointer(
        ocean.model;
        schedule = OC.TimeInterval(checkpoint_interval),
        prefix = joinpath(output_dir, filename_prefix * "_checkpoint"),
        cleanup = false,
        verbose = true,
    )

    @info "Ocean diagnostics attached:" *
          " surface ($(length(surface_outputs)) fields, every $surface_averaging_interval)," *
          " 3-D ($(length(field_outputs)) fields, every $field_averaging_interval)," *
          " checkpointer (every $checkpoint_interval)"

    return nothing
end
