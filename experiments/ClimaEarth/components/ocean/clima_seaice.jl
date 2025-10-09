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

Initially, no sea ice is present.

Since this model does not solve for prognostic temperature, we use a
prescribed heat flux boundary condition at the top, which is used to solve for
temperature at the surface of the sea ice. The surface temperature from the
previous step is provided to the coupler to be used in computing fluxes.

Specific details about the default model configuration
can be found in the documentation for `ClimaOcean.ocean_simulation`.
"""
function ClimaSeaIceSimulation(land_fraction, ocean; output_dir)
    # Initialize the sea ice with the same grid as the ocean
    grid = ocean.ocean.model.grid
    arch = grid.architecture
    advection = ocean.ocean.model.advection.T
    top_heat_boundary_condition = CO.MeltingConstrainedFluxBalance()

    # TODO use branch ss-js/top-heat-bc for this constructor
    sea_ice =
        CO.sea_ice_simulation(grid, ocean.ocean; advection, top_heat_boundary_condition)

    melting_speed = 1e-4

    # Since ocean and sea ice share the same grid, we can also share the remapping objects
    remapping = ocean.remapping

    # Before version 0.96.22, the NetCDFWriter was broken on GPU
    if arch isa OC.CPU || pkgversion(OC) >= v"0.96.22"
        # Save all tracers and velocities to a NetCDF file at daily frequency
        outputs = OC.prognostic_fields(sea_ice.model)
        jld_writer = OC.JLDWriter(
            sea_ice.model,
            outputs;
            schedule = OC.TimeInterval(86400), # Daily output
            filename = joinpath(output_dir, "seaice_diagnostics.jld2"),
            overwrite_existing = true,
            array_type = Array{Float32},
        )
        sea_ice.output_writers[:diagnostics] = jld_writer
    end

    # Allocate space for the sea ice-ocean (io) fluxes
    ocean_sea_ice_fluxes = ocean.ocean_sea_ice_fluxes

    # TODO clean this up
    # Get the area fraction from the fractional sea ice concentration
    area_fraction = sea_ice.model.ice_concentration

    # Overwrite ice fraction with the static land area fraction anywhere we have nonzero land area
    #  max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
    boundary_space = axes(ocean.area_fraction)
    FT = CC.Spaces.undertype(boundary_space)

    area_fraction_boundary_space = Interfacer.remap(area_fraction, boundary_space)
    @. area_fraction_boundary_space =
        max(min(area_fraction_boundary_space, FT(1) - land_fraction), FT(0))

    sim = ClimaSeaIceSimulation(
        sea_ice,
        area_fraction_boundary_space,
        melting_speed,
        remapping,
        ocean_sea_ice_fluxes,
    )
    return sim
end

function update_sic!(area_fraction, sea_ice)
    # TODO
end

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

# Timestep the simulation forward to time `t`
Interfacer.step!(sim::ClimaSeaIceSimulation, t) =
    OC.time_step!(sim.sea_ice, float(t) - sim.sea_ice.model.clock.time)

Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:area_fraction}) = sim.area_fraction
Interfacer.get_field(sim::ClimaSeaIceSimulation, Val(:sea_ice_concentration)) =
    sim.sea_ice.model.ice_concentration

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
    Float32(0.7)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_diffuse_albedo}) =
    Float32(0.7)

# Approximate the sea ice surface temperature as the temperature computed from the
#  fluxes at the previous timestep.
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_temperature}) =
    273.15 .+
    OC.interior(sim.sea_ice.model.ice_thermodynamics.top_surface_temperature, :, :, 1)

"""
    FluxCalculator.update_turbulent_fluxes!(sim::ClimaSeaIceSimulation, fields)

Update the turbulent fluxes in the simulation using the values stored in the coupler fields.
These include latent heat flux, sensible heat flux, momentum fluxes, and moisture flux.

Note that currently the moisture flux has no effect on the sea ice model, which has
constant salinity.

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
    # Note that this requires the sea ice model to always be run with dynamics turned on
    si_flux_u = sim.sea_ice.model.dynamics.external_stresses.top.u
    si_flux_v = sim.sea_ice.model.dynamics.external_stresses.top.v
    set_from_extrinsic_vectors!(
        (; u = si_flux_u, v = si_flux_v),
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

    # Update the sea ice only where the concentration is greater than zero.
    si_flux_heat = ice_sim.sea_ice.model.external_heat_fluxes.top
    OC.interior(si_flux_heat, :, :, 1) .=
        OC.interior(si_flux_heat, :, :, 1) .+ (
            (OC.interior(sim.sea_ice.model.ice_concentration, :, :, 1) .> 0) .*
            (remapped_F_lh .+ remapped_F_sh)
        )

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

Note that currently precipitation has no effect on the sea ice model, which has
constant salinity.

A note on sign conventions:
ClimaAtmos and ClimaSeaIce both use the convention that a positive flux is an upward flux.
No sign change is needed during the exchange, except for precipitation/salinity fluxes.
ClimaAtmos provides precipitation as a negative flux at the surface, and
ClimaSeaIce represents precipitation as a positive salinity flux,
so a sign change is needed when we convert from precipitation to salinity flux.
"""
function FieldExchanger.update_sim!(sim::ClimaSeaIceSimulation, csf, area_fraction)
    # Remap radiative flux onto scratch array; rename for clarity
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        csf.F_radiative,
    )
    remapped_F_radiative = sim.remapping.scratch_arr1

    # Update only the part due to radiative fluxes. For the full update, the component due
    # to latent and sensible heat is missing and will be updated in update_turbulent_fluxes.
    # Update the sea ice only where the concentration is greater than zero.
    si_flux_heat = sim.sea_ice.model.external_heat_fluxes.top
    OC.interior(si_flux_heat, :, :, 1) .=
        (OC.interior(sim.sea_ice.model.ice_concentration, :, :, 1) .> 0) .*
        remapped_F_radiative

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
    CO.compute_sea_ice_ocean_fluxes!(
        ice_sim.ocean_sea_ice_fluxes,
        ocean_sim.ocean,
        ice_sim.sea_ice,
        melting_speed,
        ocean_properties,
    )
    ocean_sim.ocean_sea_ice_fluxes = ice_sim.ocean_sea_ice_fluxes

    ## Update the internals of the sea ice model
    # Set the bottom heat flux to the sum of the frazil and interface heat fluxes
    bottom_heat_flux = ice_sim.sea_ice.model.external_heat_fluxes.bottom

    Qf = sea_ice_ocean_fluxes.frazil_heat        # frazil heat flux
    Qi = sea_ice_ocean_fluxes.interface_heat     # interfacial heat flux
    bottom_heat_flux .= Qf .+ Qi

    ## Update the internals of the ocean model
    ρₒ⁻¹ = 1 / ocean_sim.ocean_properties.reference_density
    cₒ = ocean_sim.ocean_properties.heat_capacity

    Jᵀio = Qio * ρₒ⁻¹ / cₒ
    Jˢio = sea_ice_ocean_fluxes.salt[i, j, 1] * ℵᵢ

    # ℑxᶠᵃᵃ: interpolate faces to centers
    τxio = ρτxio[i, j, 1] * ρₒ⁻¹ * ℑxᶠᵃᵃ(i, j, 1, grid, ℵ)
    τyio = ρτyio[i, j, 1] * ρₒ⁻¹ * ℑyᵃᶠᵃ(i, j, 1, grid, ℵ)

    # TODO finish this

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
