import Oceananigans as OC
import ClimaOcean as CO
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaComms
import ClimaCore as CC
import Thermodynamics as TD
import ClimaOcean.EN4: download_dataset
using KernelAbstractions: @kernel, @index, @inbounds

include("climaocean_helpers.jl")

"""
    ClimaSeaIceSimulation{SIM, A, OPROP, REMAP}

The ClimaCoupler simulation object used to run with ClimaSeaIce.
This type is used by the coupler to indicate that this simulation
is an surface/ocean simulation for dispatch.

It contains the following objects:
- `sea_ice::SIM`: The ClimaSeaIce simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
- `melting_speed::MS`: An constant characteristic speed for melting/freezing.
- `remapping::REMAP`: Objects needed to remap from the exchange (spectral) grid to Oceananigans spaces.
- `ocean_sea_ice_fluxes::NT`: A NamedTuple of fluxes between the ocean and sea ice, computed at each coupling step.
"""
struct ClimaSeaIceSimulation{SIM, A, MS, REMAP, NT} <: Interfacer.SeaIceModelSimulation
    sea_ice::SIM
    area_fraction::A
    melting_speed::MS
    remapping::REMAP
    ocean_sea_ice_fluxes::NT
end

"""
    ClimaSeaIceSimulation()

Creates an ClimaSeaIceSimulation object containing a model, an integrator, and
a surface area fraction field.
This type is used to indicate that this simulation is an ocean simulation for
dispatch in coupling.

Specific details about the default model configuration
can be found in the documentation for `ClimaOcean.ocean_simulation`.
"""
function ClimaSeaIceSimulation(area_fraction, ocean; output_dir)
    # Initialize the sea ice with the same grid as the ocean
    # Initially no sea ice is present
    grid = ocean.ocean.model.grid
    arch = grid.architecture
    advection = ocean.ocean.model.advection.T
    sea_ice = CO.sea_ice_simulation(grid, ocean.ocean; advection)

    melting_speed = 1e-4

    # Since ocean and sea ice share the same grid, we can also share the remapping objects
    remapping = ocean.remapping

    # Before version 0.96.22, the NetCDFWriter was broken on GPU
    # TODO this fails with OC Field MethodError
    # if arch isa OC.CPU || pkgversion(OC) >= v"0.96.22"
    #     # Save all tracers and velocities to a NetCDF file at daily frequency
    #     outputs = OC.prognostic_fields(sea_ice.model)
    #     netcdf_writer = OC.NetCDFWriter(
    #         sea_ice.model,
    #         outputs;
    #         schedule = OC.TimeInterval(86400), # Daily output
    #         filename = joinpath(output_dir, "seaice_diagnostics.nc"),
    #         indices = (:, :, grid.Nz),
    #         overwrite_existing = true,
    #         array_type = Array{Float32},
    #     )
    #     sea_ice.output_writers[:diagnostics] = netcdf_writer
    # end

    # Allocate space for the sea ice-ocean (io) fluxes
    io_bottom_heat_flux = OC.Field{Center, Center, Nothing}(grid)
    io_frazil_heat_flux = OC.Field{Center, Center, Nothing}(grid)
    io_salt_flux = OC.Field{Center, Center, Nothing}(grid)
    x_momentum = OC.Field{Face, Center, Nothing}(grid)
    y_momentum = OC.Field{Center, Face, Nothing}(grid)

    ocean_sea_ice_fluxes = (
        interface_heat = io_bottom_heat_flux,
        frazil_heat = io_frazil_heat_flux,
        salt = io_salt_flux,
        x_momentum = x_momentum,
        y_momentum = y_momentum,
    )

    sim = ClimaSeaIceSimulation(
        sea_ice,
        area_fraction,
        melting_speed,
        remapping,
        ocean_sea_ice_fluxes,
    )
    return sim
end

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

# Timestep the simulation forward to time `t`
Interfacer.step!(sim::ClimaSeaIceSimulation, t) =
    OC.time_step!(sim.sea_ice, float(t) - sim.sea_ice.model.clock.time)

Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:area_fraction}) = sim.area_fraction

# TODO: Better values for this
# At the moment, we return always Float32. This is because we always want to run
# Oceananingans with Float64, so we have no way to know the float type here. Sticking with
# Float32 ensures that nothing is accidentally promoted to Float64. We will need to change
# this anyway.
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:roughness_buoyancy}) =
    Float32(5.8e-5)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:roughness_momentum}) =
    Float32(5.8e-5)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:beta}) = Float32(1)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_direct_albedo}) =
    Float32(0.011)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_diffuse_albedo}) =
    Float32(0.069)

# Approximate the sea ice surface temperature from the bulk temperature and ocean surface temperature,
#  assuming a linear profile in the top layer
# TODO how to get ocean surface temp here? Is it in the sea ice or do we need to pass the ocean sim too?
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_temperature}) =
    273.15 + sim.ocean.model.tracers.T

"""
    FluxCalculator.update_turbulent_fluxes!(sim::ClimaSeaIceSimulation, fields)

Update the turbulent fluxes in the simulation using the values stored in the coupler fields.
These include latent heat flux, sensible heat flux, momentum fluxes, and moisture flux.

A note on sign conventions:
SurfaceFluxes and ClimaSeaIce both use the convention that a positive flux is an upward flux.
No sign change is needed during the exchange, except for moisture/salinity fluxes:
SurfaceFluxes provides moisture moving from atmosphere to ocean as a negative flux at the surface,
and ClimaSeaIce represents moisture moving from atmosphere to ocean as a positive salinity flux,
so a sign change is needed when we convert from moisture to salinity flux.
"""
function FluxCalculator.update_turbulent_fluxes!(sim::ClimaSeaIceSimulation, fields)
    # Only LatitudeLongitudeGrid are supported because otherwise we have to rotate the vectors

    (; F_lh, F_sh, F_turb_ρτxz, F_turb_ρτyz, F_turb_moisture) = fields
    grid = sim.sea_ice.model.grid

    # Remap momentum fluxes onto reduced 2D Center, Center fields using scratch arrays and fields
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        F_turb_ρτxz,
    )
    OC.set!(sim.remapping.scratch_cc1, sim.remapping.scratch_arr1) # zonal momentum flux
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr2,
        sim.remapping.remapper_cc,
        F_turb_ρτyz,
    )
    OC.set!(sim.remapping.scratch_cc2, sim.remapping.scratch_arr2) # meridional momentum flux

    # Rename for clarity; these are now Center, Center Oceananigans fields
    F_turb_ρτxz_cc = sim.remapping.scratch_cc1
    F_turb_ρτyz_cc = sim.remapping.scratch_cc2

    # Set the momentum flux BCs at the correct locations using the remapped scratch fields
    oc_flux_u = surface_flux(sim.sea_ice.model.velocities.u)
    oc_flux_v = surface_flux(sim.sea_ice.model.velocities.v)
    set_from_extrinsic_vectors!(
        (; u = oc_flux_u, v = oc_flux_v),
        grid,
        F_turb_ρτxz_cc,
        F_turb_ρτyz_cc,
    )

    # Remap the latent and sensible heat fluxes using scratch arrays
    CC.Remapping.interpolate!(sim.remapping.scratch_arr1, sim.remapping.remapper_cc, F_lh) # latent heat flux
    CC.Remapping.interpolate!(sim.remapping.scratch_arr2, sim.remapping.remapper_cc, F_sh) # sensible heat flux

    # Rename for clarity; recall F_turb_energy = F_lh + F_sh
    remapped_F_lh = sim.remapping.scratch_arr1
    remapped_F_sh = sim.remapping.scratch_arr2


    # TODO update this for sea ice
    # TODO what is sim.sea_ice.model.ice_thermodynamics.top_surface_temperature? where is it set?
    # TODO ocean_reference_density -> sea_ice.model.ice_density ?
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    OC.interior(oc_flux_T, :, :, 1) .=
        OC.interior(oc_flux_T, :, :, 1) .+
        (remapped_F_lh .+ remapped_F_sh) ./ (ocean_reference_density * ocean_heat_capacity)

    # Add the part of the salinity flux that comes from the moisture flux. We also need to
    # add the component due to precipitation (that was done with the radiative fluxes)
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        F_turb_moisture,
    )
    moisture_fresh_water_flux = sim.remapping.scratch_arr1 ./ ocean_fresh_water_density # TODO do we need this for sea ice too?
    oc_flux_S = surface_flux(sim.sea_ice.model.tracers.S)
    surface_salinity = OC.interior(sim.sea_ice.model.tracers.S, :, :, 1)
    OC.interior(oc_flux_S, :, :, 1) .=
        OC.interior(oc_flux_S, :, :, 1) .- surface_salinity .* moisture_fresh_water_flux
    return nothing
end

function Interfacer.update_field!(sim::ClimaSeaIceSimulation, ::Val{:area_fraction}, field)
    sim.area_fraction .= field
    return nothing
end

"""
    FieldExchanger.update_sim!(sim::ClimaSeaIceSimulation, csf, area_fraction)

Update the sea ice simulation with the provided fields, which have been filled in
by the coupler.

Update the portion of the surface_fluxes for T and S that is due to radiation and
precipitation. The rest will be updated in `update_turbulent_fluxes!`.

A note on sign conventions:
ClimaAtmos and ClimaSeaIce both use the convention that a positive flux is an upward flux.
No sign change is needed during the exchange, except for precipitation/salinity fluxes.
ClimaAtmos provides precipitation as a negative flux at the surface, and
ClimaSeaIce represents precipitation as a positive salinity flux,
so a sign change is needed when we convert from precipitation to salinity flux.
"""
function FieldExchanger.update_sim!(sim::ClimaSeaIceSimulation, csf, area_fraction)
    # TODO update this for sea ice
    (; ocean_reference_density, ocean_heat_capacity, ocean_fresh_water_density) =
        sim.ocean_properties

    # Remap radiative flux onto scratch array; rename for clarity
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        csf.F_radiative,
    )
    remapped_F_radiative = sim.remapping.scratch_arr1

    # Update only the part due to radiative fluxes. For the full update, the component due
    # to latent and sensible heat is missing and will be updated in update_turbulent_fluxes.
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    OC.interior(oc_flux_T, :, :, 1) .=
        remapped_F_radiative ./ (ocean_reference_density * ocean_heat_capacity)

    # Remap precipitation fields onto scratch arrays; rename for clarity
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        csf.P_liq,
    )
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr2,
        sim.remapping.remapper_cc,
        csf.P_snow,
    )
    remapped_P_liq = sim.remapping.scratch_arr1
    remapped_P_snow = sim.remapping.scratch_arr2

    # Virtual salt flux
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    precipitating_fresh_water_flux =
        (remapped_P_liq .+ remapped_P_snow) ./ ocean_fresh_water_density
    surface_salinity_flux =
        OC.interior(sim.ocean.model.tracers.S, :, :, 1) .* precipitating_fresh_water_flux
    OC.interior(oc_flux_S, :, :, 1) .= .-surface_salinity_flux
    return nothing
end

"""
    ocean_seaice_fluxes!(ocean_sim, ice_sim)

Compute the fluxes between the ocean and sea ice, storing them in the `ocean_sea_ice_fluxes`
fields of the ocean and sea ice simulations.
"""
function FluxCalculator.ocean_seaice_fluxes!(
    ocean_sim::OceananigansSimulation,
    ice_sim::ClimaSeaIceSimulation,
)
    # TODO unify sea_ice vs ice naming convention in this file
    melting_speed = ice_sim.melting_speed
    ocean_properties = ocean_sim.ocean_properties

    # Compute the fluxes and store them in the both simulations
    OC.compute_sea_ice_ocean_fluxes!(
        ice_sim.ocean_sea_ice_fluxes,
        ocean_sim.ocean,
        ice_sim.sea_ice,
        melting_speed,
        ocean_properties,
    )
    ocean_sim.ocean_sea_ice_fluxes = ice_sim.ocean_sea_ice_fluxes

    # TODO what do we do with these fluxes now? They need to be passed to the component sims somehow
    return nothing
end


"""
    get_model_prog_state(sim::ClimaSeaIceSimulation)

Returns the model state of a simulation as a `ClimaCore.FieldVector`.
It's okay to leave this unimplemented for now, but we won't be able to use the
restart system.

TODO extend this for non-ClimaCore states.
"""
function Checkpointer.get_model_prog_state(sim::ClimaSeaIceSimulation)
    @warn "get_model_prog_state not implemented for ClimaSeaIceSimulation"
end
