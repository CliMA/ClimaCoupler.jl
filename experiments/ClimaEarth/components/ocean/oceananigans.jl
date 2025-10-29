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
    OceananigansSimulation{SIM, A, OPROP, REMAP, SIC}

The ClimaCoupler simulation object used to run with Oceananigans.
This type is used by the coupler to indicate that this simulation
is an surface/ocean simulation for dispatch.

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
    area_fraction,
    start_date,
    stop_date;
    output_dir,
    comms_ctx = ClimaComms.context(),
)
    arch = comms_ctx.device isa ClimaComms.CUDADevice ? OC.GPU() : OC.CPU()

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
        interpolation_passes = 1, # TODO revert
        major_basins = 1,
    )

    grid = OC.ImmersedBoundaryGrid(
        underlying_grid,
        OC.GridFittedBottom(bottom_height);
        active_cells_map = true,
    )

    # TODO use_restoring only if not using ClimaSeaIce
    use_restoring = start_date + Dates.Month(1) < stop_date

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
    momentum_advection = OC.WENOVectorInvariant(order = 5)
    tracer_advection = OC.WENO(order = 5)
    eddy_closure = OC.TurbulenceClosures.IsopycnalSkewSymmetricDiffusivity(
        κ_skew = 1e3,
        κ_symmetric = 1e3,
    )
    vertical_mixing = CO.OceanSimulations.default_ocean_closure()
    horizontal_viscosity = OC.HorizontalScalarBiharmonicDiffusivity(ν = 1e11)

    ocean = CO.ocean_simulation(
        grid;
        forcing,
        momentum_advection,
        tracer_advection,
        free_surface,
        closure = (eddy_closure, horizontal_viscosity, vertical_mixing),
    )

    # Set initial condition to EN4 state estimate at start_date
    OC.set!(ocean.model, T = en4_temperature[1], S = en4_salinity[1])

    long_cc = OC.λnodes(grid, OC.Center(), OC.Center(), OC.Center())
    lat_cc = OC.φnodes(grid, OC.Center(), OC.Center(), OC.Center())

    # TODO: Go from 0 to Nx+1, Ny+1 (for halos) (for LatLongGrid)

    # Construct a remapper from the exchange grid to `Center, Center` fields
    long_cc = reshape(long_cc, length(long_cc), 1)
    lat_cc = reshape(lat_cc, 1, length(lat_cc))
    target_points_cc = @. CC.Geometry.LatLongPoint(lat_cc, long_cc)
    # TODO: We can remove the `nothing` after CC > 0.14.33
    remapper_cc = CC.Remapping.Remapper(axes(area_fraction), target_points_cc, nothing)

    # Construct two 2D Center/Center fields to use as scratch space while remapping
    scratch_cc1 = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    scratch_cc2 = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    # Construct two scratch arrays to use while remapping
    # We get the array type, float type, and dimensions from the remapper object to maintain consistency
    ArrayType = ClimaComms.array_type(remapper_cc.space)
    FT = CC.Spaces.undertype(remapper_cc.space)
    interpolated_values_dim..., _buffer_length = size(remapper_cc._interpolated_values)
    scratch_arr1 = ArrayType(zeros(FT, interpolated_values_dim...))
    scratch_arr2 = ArrayType(zeros(FT, interpolated_values_dim...))
    scratch_arr3 = ArrayType(zeros(FT, interpolated_values_dim...))

    remapping =
        (; remapper_cc, scratch_cc1, scratch_cc2, scratch_arr1, scratch_arr2, scratch_arr3)

    ocean_properties =
        (; reference_density = 1020, heat_capacity = 3991, fresh_water_density = 999.8)

    # Before version 0.96.22, the NetCDFWriter was broken on GPU
    if arch isa OC.CPU || pkgversion(OC) >= v"0.96.22"
        # Save all tracers and velocities to a NetCDF file at daily frequency
        outputs = merge(ocean.model.tracers, ocean.model.velocities)
        netcdf_writer = OC.NetCDFWriter(
            ocean.model,
            outputs;
            schedule = OC.TimeInterval(86400), # Daily output
            filename = joinpath(output_dir, "ocean_diagnostics.nc"),
            indices = (:, :, grid.Nz),
            overwrite_existing = true,
            array_type = Array{Float32},
        )
        ocean.output_writers[:diagnostics] = netcdf_writer
    end

    # Initialize with 0 ice concentration; this will be updated in `resolve_ocean_ice_fractions!`
    # if the ocean is coupled to a non-prescribed sea ice model.
    ice_concentration = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    return OceananigansSimulation(
        ocean,
        area_fraction,
        ocean_properties,
        remapping,
        ice_concentration,
    )
end

"""
    FieldExchanger.resolve_ocean_ice_fractions!(ocean_sim, ice_sim, land_fraction)

Ensure the ocean and ice area fractions are consistent with each other.
This matters in the case of a LatitudeLongitudeGrid, which is only
defined between -80 and 80 degrees latitude. In this case, we want to
set the ice fraction to `1 - land_fraction` on [-90, -80] and [80, 90]
degrees latitude, and make sure the ocean fraction is 0 there.

The land fraction is expected to be set to 1 at the poles before calling this function,
and doesn't need to be set again since its fraction is static.

This function also updates the ice concentration field in the ocean simulation
so that it can be used for weighting flux updates.
"""
function FieldExchanger.resolve_ocean_ice_fractions!(
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
        polar_mask .= abs.(lat) .>= FT(80)

        # TODO do we want both to be 0 since we use capped lat/lon?
        # Set ice fraction to 1 - land_fraction and ocean fraction to 0 where polar_mask is 1
        @. ice_fraction = ifelse.(polar_mask == FT(1), FT(1) - land_fraction, ice_fraction)
        @. ocean_fraction = ifelse.(polar_mask == FT(1), FT(0), ocean_fraction)
    end

    # Update the ice concentration field in the ocean simulation
    ocean_sim.ice_concentration .= Interfacer.get_field(ice_sim, Val(:ice_concentration))
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
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:beta}) = Float32(1)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:emissivity}) = Float32(0.97)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_direct_albedo}) =
    Float32(0.011)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_diffuse_albedo}) =
    Float32(0.069)

# NOTE: This is 3D, but it will be remapped to 2D
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_temperature}) =
    273.15 + sim.ocean.model.tracers.T

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
    grid = sim.ocean.model.grid
    ice_concentration = sim.ice_concentration

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

    # Weight by (1 - sea ice concentration)
    # TODO does this work with OC fields?
    OC.interior(F_turb_ρτxz_cc, :, :, 1) .=
        OC.interior(F_turb_ρτxz_cc, :, :, 1) .* (1.0 .- ice_concentration)
    OC.interior(F_turb_ρτyz_cc, :, :, 1) .=
        OC.interior(F_turb_ρτyz_cc, :, :, 1) .* (1.0 .- ice_concentration)

    # Set the momentum flux BCs at the correct locations using the remapped scratch fields
    oc_flux_u = surface_flux(sim.ocean.model.velocities.u)
    oc_flux_v = surface_flux(sim.ocean.model.velocities.v)
    set_from_extrinsic_vector!(
        (; u = oc_flux_u, v = oc_flux_v),
        grid,
        F_turb_ρτxz_cc,
        F_turb_ρτyz_cc,
    )

    (; reference_density, heat_capacity, fresh_water_density) = sim.ocean_properties

    # Remap the latent and sensible heat fluxes using scratch arrays
    CC.Remapping.interpolate!(sim.remapping.scratch_arr1, sim.remapping.remapper_cc, F_lh) # latent heat flux
    CC.Remapping.interpolate!(sim.remapping.scratch_arr2, sim.remapping.remapper_cc, F_sh) # sensible heat flux

    # Rename for clarity; recall F_turb_energy = F_lh + F_sh
    remapped_F_lh = sim.remapping.scratch_arr1
    remapped_F_sh = sim.remapping.scratch_arr2

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
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        F_turb_moisture,
    )
    moisture_fresh_water_flux = sim.remapping.scratch_arr1 ./ fresh_water_density
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    surface_salinity = OC.interior(sim.ocean.model.tracers.S, :, :, 1)
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
    (; reference_density, heat_capacity, fresh_water_density) = sim.ocean_properties
    ice_concentration = sim.ice_concentration

    # Remap radiative flux onto scratch array; rename for clarity
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        csf.SW_d,
    )
    remapped_SW_d = sim.remapping.scratch_arr1

    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr2,
        sim.remapping.remapper_cc,
        csf.LW_d,
    )
    remapped_LW_d = sim.remapping.scratch_arr2

    # Update only the part due to radiative fluxes. For the full update, the component due
    # to latent and sensible heat is missing and will be updated in update_turbulent_fluxes.
    oc_flux_T = surface_flux(sim.ocean.model.tracers.T)
    # TODO: get sigma from parameters
    σ = 5.67e-8
    α = Interfacer.get_field(sim, Val(:surface_direct_albedo)) # scalar
    ϵ = Interfacer.get_field(sim, Val(:emissivity)) # scalar
    OC.interior(oc_flux_T, :, :, 1) .=
        (1.0 .- ice_concentration) .* (
            -(1 - α) .* remapped_SW_d .-
            ϵ * (
                remapped_LW_d .-
                σ .* (273.15 .+ OC.interior(sim.ocean.model.tracers.T, :, :, 1)) .^ 4
            )
        ) ./ (reference_density * heat_capacity)

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
    OC.interior(oc_flux_S, :, :, 1) .=
        OC.interior(oc_flux_S, :, :, 1) .-
        OC.interior(sim.ocean.model.tracers.S, :, :, 1) .* (1.0 .- ice_concentration) .*
        (remapped_P_liq .+ remapped_P_snow) ./ fresh_water_density
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
