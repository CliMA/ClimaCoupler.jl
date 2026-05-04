"""
    add_seaice_diagnostics!(ice_sim::ClimaSeaIceSimulation;
                            output_dir = ".",
                            surface_averaging_interval = 1days,
                            checkpoint_interval = 90days,
                            filename_prefix = "seaice",
                            file_splitting_interval = 15days)

Attach averaged-output writers to the underlying ClimaSeaIce simulation inside a `ClimaSeaIceSimulation`. 
Two writers are added to `ice_sim.ice.output_writers`:

1. **Surface diagnostics** (`<prefix>_surface.jld2`): sea-ice concentration, thickness, velocities, and (when present) top 
   surface temperature, averaged over `surface_averaging_interval`.
2. **Checkpointer** (`<prefix>_checkpoint`): JLD2 checkpoint of the sea-ice model at `checkpoint_interval`.
"""
function add_seaice_diagnostics!(
    ice_sim::ClimaSeaIceSimulation;
    output_dir = ".",
    surface_averaging_interval = 1days,
    checkpoint_interval = 90days,
    filename_prefix = "seaice",
    file_splitting_interval = 15days,
)
    ice = ice_sim.ice
    file_splitting = OC.TimeInterval(file_splitting_interval)

    hi = ice.model.ice_thickness
    ℵi = ice.model.ice_concentration
    ui, vi = ice.model.velocities

    sitemptop = try
        ice.model.ice_thermodynamics.top_surface_temperature
    catch
        nothing
    end

    surface_outputs = Dict{Symbol, Any}(
        :ice_concentration => ℵi,
        :ice_thickness => hi,
        :ice_zonal_velocity => ui,
        :ice_meridional_velocity => vi,
    )

    if !isnothing(sitemptop)
        surface_outputs[:ice_top_temperature] = sitemptop
    end

    ice.output_writers[:surface_averages] = OC.JLD2Writer(
        ice.model,
        surface_outputs;
        schedule = OC.AveragedTimeInterval(surface_averaging_interval),
        filename = joinpath(output_dir, filename_prefix * "_surface"),
        file_splitting,
        overwrite_existing = true,
    )

    ice.output_writers[:checkpointer] = OC.Checkpointer(
        ice.model;
        schedule = OC.TimeInterval(checkpoint_interval),
        prefix = joinpath(output_dir, filename_prefix * "_checkpoint"),
        cleanup = false,
        verbose = true,
    )

    @info "Sea-ice diagnostics attached:" *
          " surface ($(length(surface_outputs)) fields, every $(prettytime(surface_averaging_interval)))," *
          " checkpointer (every $(prettytime(checkpoint_interval)))"

    return nothing
end
