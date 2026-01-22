import Oceananigans as OC
import ClimaOcean as CO
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaComms
import ClimaCore as CC
import Thermodynamics as TD
import ClimaParams as CP
import ClimaOcean.EN4: download_dataset
using KernelAbstractions: @kernel, @index, @inbounds
import ConservativeRegridding as CR

include("climaocean_helpers.jl")
include("remapping.jl")

"""
    OceananigansSimulation{SIM, A, OPROP, REMAP, SIC}

The ClimaCoupler simulation object used to run with Oceananigans.
This type is used by the coupler to indicate that this simulation
is a surface/ocean simulation for dispatch.

It contains the following objects:
- `ocean::SIM`: The Oceananigans simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
- `ocean_properties::OPROP`: A NamedTuple of ocean properties and parameters
- `remapping::REMAP`: Objects needed to remap from the exchange (spectral) grid to Oceananigans spaces.
- `ice_concentration::SIC`: An Oceananigans Field representing the sea ice concentration on the ocean/sea ice grid.
"""
struct OceananigansSimulation{SIM, A, OPROP, REMAP, SIC} <: Interfacer.OceanModelSimulation
    ocean::SIM
    area_fraction::A
    ocean_properties::OPROP
    remapping::REMAP
    ice_concentration::SIC
end

"""
    OceananigansSimulation()

Creates an OceananigansSimulation object containing a model, an integrator, and
a surface area fraction field.
This type is used to indicate that this simulation is an ocean simulation for
dispatch in coupling.

Specific details about the default model configuration
can be found in the documentation for `ClimaOcean.ocean_simulation`.
"""
function OceananigansSimulation(
    ::Type{FT},
    boundary_space,
    start_date,
    stop_date;
    output_dir,
    ice_model,
    Δt = nothing,
    comms_ctx = ClimaComms.context(),
    coupled_param_dict = CP.create_toml_dict(eltype(area_fraction)),
) where {FT}
    arch = comms_ctx.device isa ClimaComms.CUDADevice ? OC.GPU() : OC.CPU()

    # Use Float64 for the ocean to avoid precision issues
    FT_ocean = Float64
    OC.Oceananigans.defaults.FloatType = FT_ocean

    # Retrieve EN4 data (monthly)
    # (It requires username and password)
    dates = range(start_date, step = Dates.Month(1), stop = stop_date)
    en4_temperature = CO.Metadata(:temperature; dates, dataset = CO.EN4.EN4Monthly())
    en4_salinity = CO.Metadata(:salinity; dates, dataset = CO.EN4.EN4Monthly())
    download_dataset(en4_temperature)
    download_dataset(en4_salinity)

    # Set up ocean grid (1 degree)
    resolution_points = (360, 160, 32)
    Nz = last(resolution_points)
    depth = 4000 # meters
    z = OC.ExponentialDiscretization(Nz, -depth, 0; scale = 0.85 * depth)

    # Regular LatLong because we know how to do interpolation there
    underlying_grid = OC.LatitudeLongitudeGrid(
        arch;
        size = resolution_points,
        longitude = (-180, 180),
        latitude = (-80, 80),   # NOTE: Don't goo to high up when using LatLongGrid, or the cells will be too small
        z,
        halo = (7, 7, 7),
    )

    bottom_height = CO.regrid_bathymetry(
        underlying_grid;
        minimum_depth = 30,
        interpolation_passes = 20,
        major_basins = 1,
    )

    grid = OC.ImmersedBoundaryGrid(
        underlying_grid,
        OC.GridFittedBottom(bottom_height);
        active_cells_map = true,
    )

    # Restore the ocean to the EN4 state periodically if running for more than one month and not using ClimaSeaIce
    use_restoring = start_date + Dates.Month(1) < stop_date && ice_model != "clima_seaice"

    if use_restoring
        # When we use EN4 data, the forcing takes care of everything, including
        # the initial conditions
        restoring_rate = 1 / (3 * 86400)
        mask = CO.LinearlyTaperedPolarMask(
            southern = (-80, -70),
            northern = (70, 90),
            z = (z(1), 0),
        )

        forcing_T = CO.DatasetRestoring(en4_temperature, grid; mask, rate = restoring_rate)
        forcing_S = CO.DatasetRestoring(en4_salinity, grid; mask, rate = restoring_rate)
        forcing = (T = forcing_T, S = forcing_S)
    else
        forcing = (;)
    end

    # Create ocean simulation
    free_surface = OC.SplitExplicitFreeSurface(grid; substeps = 70)
    momentum_advection = OC.VectorInvariant()
    tracer_advection = OC.WENO(order = 5)
    horizontal_viscosity = OC.HorizontalScalarBiharmonicDiffusivity(ν = 1e11)

    # Use Float32 for the vertical mixing parameters to avoid parameter memory limits
    vertical_mixing = OC.CATKEVerticalDiffusivity(Float32)

    Δt = isnothing(Δt) ? CO.OceanSimulations.estimate_maximum_Δt(grid) : Δt
    ocean = CO.ocean_simulation(
        grid;
        Δt,
        forcing,
        momentum_advection,
        tracer_advection,
        free_surface,
        closure = (horizontal_viscosity, vertical_mixing),
    )

    # Set initial condition to EN4 state estimate at start_date
    OC.set!(ocean.model, T = en4_temperature[1], S = en4_salinity[1])

    # Construct the remapper object and allocate scratch space
    remapping = construct_remappers(grid, boundary_space)

    # Get some ocean properties and parameters
    ocean_properties = (;
        reference_density = 1020,
        heat_capacity = 3991,
        σ = coupled_param_dict["stefan_boltzmann_constant"],
        C_to_K = coupled_param_dict["temperature_water_freeze"],
    )

    # Before version 0.96.22, the NetCDFWriter was broken on GPU
    if arch isa OC.CPU
        # Save all tracers and velocities to a NetCDF file at daily frequency
        outputs = merge(ocean.model.tracers, ocean.model.velocities)
        surface_writer = OC.NetCDFWriter(
            ocean.model,
            outputs;
            schedule = OC.TimeInterval(86400), # Daily output
            filename = joinpath(output_dir, "ocean_diagnostics.nc"),
            indices = (:, :, grid.Nz),
            overwrite_existing = true,
            array_type = Array{Float32},
        )
        free_surface_writer = OC.NetCDFWriter(
            ocean.model,
            (; η = ocean.model.free_surface.η); # The free surface (.η) will change to .displacement after version 0.104.0
            schedule = OC.TimeInterval(3600), # hourly snapshots
            filename = joinpath(output_dir, "ocean_free_surface.nc"),
            overwrite_existing = true,
            array_type = Array{Float32},
        )
        Tflux = ocean.model.tracers.T.boundary_conditions.top.condition
        Sflux = ocean.model.tracers.S.boundary_conditions.top.condition
        uflux = ocean.model.velocities.u.boundary_conditions.top.condition
        vflux = ocean.model.velocities.v.boundary_conditions.top.condition
        fluxes_writer = OC.NetCDFWriter(
            ocean.model,
            (; Tflux, Sflux, uflux, vflux);
            schedule = OC.TimeInterval(3600), # hourly snapshots
            filename = joinpath(output_dir, "ocean_fluxes.nc"),
            overwrite_existing = true,
            array_type = Array{Float32},
        )

        ocean.output_writers[:surface] = surface_writer
        ocean.output_writers[:free_surface] = free_surface_writer
        ocean.output_writers[:fluxes] = fluxes_writer
    end

    # Initialize with 0 ice concentration; this will be updated in `resolve_area_fractions!`
    # if the ocean is coupled to a non-prescribed sea ice model.
    ice_concentration = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    # Create a dummy area fraction that will get overwritten in `update_surface_fractions!`
    area_fraction = ones(boundary_space)

    return OceananigansSimulation(
        ocean,
        area_fraction,
        ocean_properties,
        remapping,
        ice_concentration,
    )
end

"""
    construct_remappers(grid_oc, boundary_space)

Given an Oceananigans grid and a ClimaCore boundary space, construct the
remappers needed to remap between the two grids in both directions.

Returns a remapper from the Oceananigans grid to the ClimaCore boundary space.
To regrid from Oceananigans to ClimaCore, use `CR.regrid!(dest_vector, remapper_oc_to_cc, src_vector)`.
To regrid from ClimaCore to Oceananigans, use `CR.regrid!(dest_vector, transpose(remapper_oc_to_cc), src_vector)`.
"""
function construct_remappers(grid_oc, boundary_space)
    # Create the remapper from the Oceananigans grid to the ClimaCore boundary space
    remapper_oc_to_cc =
        CR.Regridder(boundary_space, grid_oc; normalize = false, threaded = false)

    # Create a field of ones on the boundary space so we can compute element areas
    field_ones_cc = CC.Fields.ones(boundary_space)

    # Allocate a vector with length equal to the number of elements in the target space
    # To be used as a temp field for remapping
    FT = CC.Spaces.undertype(boundary_space)
    value_per_element_cc = zeros(FT, CC.Meshes.nelements(boundary_space.grid.topology.mesh))

    # Construct two 2D Oceananigans Center/Center fields to use as scratch space while remapping
    scratch_field_oc1 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)
    scratch_field_oc2 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)

    # Allocate space for a Field of UVVectors, which we need for remapping momentum fluxes
    temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)
    return (;
        remapper_oc_to_cc,
        field_ones_cc,
        value_per_element_cc,
        scratch_field_oc1,
        scratch_field_oc2,
        temp_uv_vec,
    )
end

"""
    FieldExchanger.resolve_area_fractions!(ocean_sim, ice_sim, land_fraction)

Ensure the ocean and ice area fractions are consistent with each other.
This matters in the case of a LatitudeLongitudeGrid, which is only
defined between -80 and 80 degrees latitude. In this case, we set the ice
and ocean area fractions to 0 and the land fraction to 1 on [78°S, 90°S]
and on [78°N, 90°N]. Note the overlap between 78° and 80° to avoid any gaps
introduced by the cubed sphere not aligning with latitude lines.

This function also updates the ice concentration field in the ocean simulation
so that it can be used for weighting flux updates.
"""
function FieldExchanger.resolve_area_fractions!(
    ocean_sim::OceananigansSimulation,
    ice_sim,
    land_fraction,
)
    if ocean_sim.ocean.model.grid.underlying_grid isa OC.LatitudeLongitudeGrid
        ocean_fraction = Interfacer.get_field(ocean_sim, Val(:area_fraction))
        ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction))

        # Create a "polar" mask that's 1 at latitudes in [-90, -80] and [80, 90] degrees
        boundary_space = axes(ocean_fraction)
        FT = CC.Spaces.undertype(boundary_space)
        lat = CC.Fields.coordinate_field(boundary_space).lat
        polar_mask = CC.Fields.zeros(boundary_space)
        polar_mask .= abs.(lat) .>= FT(78)

        # Set land fraction to 1 and ice/ocean fraction to 0 where polar_mask is 1
        @. land_fraction = ifelse.(polar_mask == FT(1), FT(1), land_fraction)
        @. ice_fraction = ifelse.(polar_mask == FT(1), FT(0), ice_fraction)
        @. ocean_fraction = ifelse.(polar_mask == FT(1), FT(0), ocean_fraction)
    end

    # Update the ice concentration field in the ocean simulation
    ice_sim isa ClimaSeaIceSimulation && (
        ocean_sim.ice_concentration .=
            Interfacer.get_field(ice_sim, Val(:ice_concentration))
    )
    return nothing
end

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

# Timestep the simulation forward to time `t`
Interfacer.step!(sim::OceananigansSimulation, t) =
    OC.time_step!(sim.ocean, float(t) - sim.ocean.model.clock.time)

Interfacer.get_field(sim::OceananigansSimulation, ::Val{:area_fraction}) = sim.area_fraction

# TODO: Better values for this
# At the moment, we return always Float32. This is because we always want to run
# Oceananingans with Float64, so we have no way to know the float type here. Sticking with
# Float32 ensures that nothing is accidentally promoted to Float64. We will need to change
# this anyway.
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_buoyancy}) =
    Float32(5.8e-5)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_momentum}) =
    Float32(5.8e-5)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:emissivity}) = Float32(0.97)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_direct_albedo}) =
    Float32(0.011)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_diffuse_albedo}) =
    Float32(0.069)

# NOTE: This is 3D, but it will be remapped to 2D
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_temperature}) =
    sim.ocean.model.tracers.T + sim.ocean_properties.C_to_K # convert from Celsius to Kelvin

"""
    FluxCalculator.update_turbulent_fluxes!(sim::OceananigansSimulation, fields)

Update the turbulent fluxes in the simulation using the values computed at this time step.
These include latent heat flux, sensible heat flux, momentum fluxes, and moisture flux.

Rather than setting the surface fluxes and overwriting previous values, this function adds only
the contributions from the turbulent fluxes. `update_sim!` sets the surface fluxes due to
radiation and precipitation. Additional contributions may be made in `ocean_seaice_fluxes!`.
An exception is the momentum fluxes, which are set directly here since they are not updated
in `update_sim!`.

A note on sign conventions:
SurfaceFluxes and Oceananigans both use the convention that a positive flux is an upward flux.
No sign change is needed during the exchange, except for moisture/salinity fluxes:
SurfaceFluxes provides moisture moving from atmosphere to ocean as a negative flux at the surface,
and Oceananigans represents moisture moving from atmosphere to ocean as a positive salinity flux,
so a sign change is needed when we convert from moisture to salinity flux.
"""
function FluxCalculator.update_turbulent_fluxes!(sim::OceananigansSimulation, fields)
    (; F_lh, F_sh, F_turb_ρτxz, F_turb_ρτyz, F_turb_moisture) = fields
    (; reference_density, heat_capacity) = sim.ocean_properties
    grid = sim.ocean.model.grid
    ice_concentration = sim.ice_concentration

    # Convert the momentum fluxes from contravariant to Cartesian basis
    contravariant_to_cartesian!(sim.remapping.temp_uv_vec, F_turb_ρτxz, F_turb_ρτyz)
    F_turb_ρτxz_uv = sim.remapping.temp_uv_vec.components.data.:1
    F_turb_ρτyz_uv = sim.remapping.temp_uv_vec.components.data.:2

    # Remap momentum fluxes onto reduced 2D Center, Center fields
    Interfacer.remap!(sim.remapping.scratch_field_oc1, F_turb_ρτxz, sim.remapping) # zonal momentum flux
    Interfacer.remap!(sim.remapping.scratch_field_oc2, F_turb_ρτyz, sim.remapping) # meridional momentum flux

    # Rename for clarity; these are now Center, Center Oceananigans fields
    oc_F_turb_ρτxz = sim.remapping.scratch_field_oc1
    oc_F_turb_ρτyz = sim.remapping.scratch_field_oc2

    # Weight by (1 - sea ice concentration)
    OC.interior(oc_F_turb_ρτxz, :, :, 1) .=
        OC.interior(oc_F_turb_ρτxz, :, :, 1) .* (1.0 .- ice_concentration) ./
        reference_density
    OC.interior(oc_F_turb_ρτyz, :, :, 1) .=
        OC.interior(oc_F_turb_ρτyz, :, :, 1) .* (1.0 .- ice_concentration) ./
        reference_density

    # Set the momentum flux BCs at the correct locations using the remapped scratch fields
    oc_flux_u = surface_flux(sim.ocean.model.velocities.u)
    oc_flux_v = surface_flux(sim.ocean.model.velocities.v)
    set_from_extrinsic_vector!(
        (; u = oc_flux_u, v = oc_flux_v),
        grid,
        oc_F_turb_ρτxz,
        oc_F_turb_ρτyz,
    )

    # Remap the latent and sensible heat fluxes using scratch arrays
    Interfacer.remap!(sim.remapping.scratch_field_oc1, F_lh, sim.remapping) # latent heat flux
    Interfacer.remap!(sim.remapping.scratch_field_oc2, F_sh, sim.remapping) # sensible heat flux

    # Rename for clarity; recall F_turb_energy = F_lh + F_sh
    remapped_F_lh = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)
    remapped_F_sh = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)

    # TODO: Note, SW radiation penetrates the surface. Right now, we just put
    # everything on the surface, but later we will need to account for this.
    # One way we can do this is using directly ClimaOcean
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    OC.interior(oc_flux_T, :, :, 1) .=
        OC.interior(oc_flux_T, :, :, 1) .+
        (1.0 .- ice_concentration) .* (remapped_F_lh .+ remapped_F_sh) ./
        (reference_density * heat_capacity)

    # Add the part of the salinity flux that comes from the moisture flux, we also need to
    # add the component due to precipitation (that was done with the radiative fluxes)
    Interfacer.remap!(sim.remapping.scratch_field_oc1, F_turb_moisture, sim.remapping) # moisture flux
    moisture_fresh_water_flux =
        OC.interior(sim.remapping.scratch_field_oc1, :, :, 1) ./ reference_density
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    surface_salinity = OC.interior(sim.ocean.model.tracers.S, :, :, grid.Nz)
    OC.interior(oc_flux_S, :, :, 1) .=
        OC.interior(oc_flux_S, :, :, 1) .-
        (1.0 .- ice_concentration) .* surface_salinity .* moisture_fresh_water_flux
    return nothing
end

function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:area_fraction}, field)
    sim.area_fraction .= field
    return nothing
end

"""
    FieldExchanger.update_sim!(sim::OceananigansSimulation, csf)

Update the ocean simulation with the provided fields, which have been filled in
by the coupler.

Update the portion of the surface_fluxes for T and S that is due to radiation and
precipitation. The rest will be updated in `update_turbulent_fluxes!`.

This function sets the surface fluxes directly, overwriting any previous values.
Additional contributions will be made in `update_turbulent_fluxes!` and `ocean_seaice_fluxes!`.

A note on sign conventions:
ClimaAtmos and Oceananigans both use the convention that a positive flux is an upward flux.
No sign change is needed during the exchange, except for precipitation/salinity fluxes.
ClimaAtmos provides precipitation as a negative flux at the surface, and
Oceananigans represents precipitation as a positive salinity flux,
so a sign change is needed when we convert from precipitation to salinity flux.
"""
function FieldExchanger.update_sim!(sim::OceananigansSimulation, csf)
    (; reference_density, heat_capacity) = sim.ocean_properties
    ice_concentration = sim.ice_concentration
    Nz = sim.ocean.model.grid.Nz

    # Remap radiative flux onto scratch array; rename for clarity
    Interfacer.remap!(sim.remapping.scratch_field_oc1, csf.SW_d, sim.remapping) # shortwave radiation
    remapped_SW_d = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)

    Interfacer.remap!(sim.remapping.scratch_field_oc2, csf.LW_d, sim.remapping) # longwave radiation
    remapped_LW_d = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)

    # Update only the part due to radiative fluxes. For the full update, the component due
    # to latent and sensible heat is missing and will be updated in update_turbulent_fluxes.
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    (; σ, C_to_K) = sim.ocean_properties
    α = Interfacer.get_field(sim, Val(:surface_direct_albedo)) # scalar
    ϵ = Interfacer.get_field(sim, Val(:emissivity)) # scalar
    OC.interior(oc_flux_T, :, :, 1) .=
        (1.0 .- ice_concentration) .* (
            -(1 - α) .* remapped_SW_d .-
            ϵ * (
                remapped_LW_d .-
                σ .* (C_to_K .+ OC.interior(sim.ocean.model.tracers.T, :, :, Nz)) .^ 4
            )
        ) ./ (reference_density * heat_capacity)

    # Remap precipitation fields onto scratch fields; rename for clarity
    Interfacer.remap!(sim.remapping.scratch_field_oc1, csf.P_liq, sim.remapping) # liquid precipitation
    remapped_P_liq = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)

    Interfacer.remap!(sim.remapping.scratch_field_oc2, csf.P_snow, sim.remapping) # snow precipitation
    remapped_P_snow = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)

    # Virtual salt flux
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    OC.interior(oc_flux_S, :, :, 1) .=
        OC.interior(oc_flux_S, :, :, 1) .-
        OC.interior(sim.ocean.model.tracers.S, :, :, Nz) .* (1.0 .- ice_concentration) .*
        (remapped_P_liq .+ remapped_P_snow) ./ reference_density
    return nothing
end

"""
    get_model_prog_state(sim::OceananigansSimulation)

Returns the model state of a simulation as a `ClimaCore.FieldVector`.
It's okay to leave this unimplemented for now, but we won't be able to use the
restart system.

TODO extend this for non-ClimaCore states.
"""
function Checkpointer.get_model_prog_state(sim::OceananigansSimulation)
    @warn "get_model_prog_state not implemented for OceananigansSimulation"
end
