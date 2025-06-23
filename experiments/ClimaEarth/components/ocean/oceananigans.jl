import Oceananigans as OC
import ClimaOcean as CO
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaComms
import ClimaCore as CC
import Thermodynamics as TD
import ClimaOcean.EN4: download_dataset
using KernelAbstractions: @kernel, @index, @inbounds

"""
    OceananigansSimulation{SIM, A, OPROP, REMAP}

The ClimaCoupler simulation object used to run with Oceananigans.
This type is used by the coupler to indicate that this simulation
is an surface/ocean simulation for dispatch.

It contains the following objects:
- `ocean::SIM`: The Oceananigans simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
- `ocean_properties::OPROP`: A NamedTuple of ocean properties and parameters
- `remapping::REMAP`: Objects needed to remap from the exchange (spectral) grid to Oceananigans spaces.
"""
struct OceananigansSimulation{SIM, A, OPROP, REMAP} <: Interfacer.OceanModelSimulation
    ocean::SIM
    area_fraction::A
    ocean_properties::OPROP
    remapping::REMAP
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

    # Set up tripolar ocean grid (1 degree)
    Nx = 360
    Ny = 180
    Nz = 40
    depth = 4000 # meters
    z = OC.ExponentialDiscretization(Nz, -depth, 0; scale = 0.85 * depth)

    underlying_grid = OC.TripolarGrid(arch; size = (Nx, Ny, Nz), halo = (7, 7, 4), z)
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
    ocean = CO.ocean_simulation(grid; forcing)

    # Set initial condition to EN4 state estimate at start_date
    OC.set!(ocean.model, T = en4_temperature[1], S = en4_salinity[1])

    # Construct a remapper from the exchange grid to `Center, Center` fields
    long_cc = OC.λnodes(grid, OC.Center(), OC.Center(), OC.Center())
    lat_cc = OC.φnodes(grid, OC.Center(), OC.Center(), OC.Center())

    # Create a 2D matrix containing each lat/long combination as a LatLongPoint
    # Note this must be done on CPU since the CC.Remapper module is not GPU-compatible
    target_points_cc = Array(CC.Geometry.LatLongPoint.(lat_cc, long_cc))

    if pkgversion(CC) >= v"0.14.34"
        remapper_cc = CC.Remapping.Remapper(axes(area_fraction), target_points_cc)
    else
        remapper_cc = CC.Remapping.Remapper(axes(area_fraction), target_points_cc, nothing)
    end

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

    remapping = (; remapper_cc, scratch_cc1, scratch_cc2, scratch_arr1, scratch_arr2)

    ocean_properties = (;
        ocean_reference_density = 1020,
        ocean_heat_capacity = 3991,
        ocean_fresh_water_density = 999.8,
    )

    # Before version 0.96.22, the NetCDFWriter was broken on GPU
    if arch isa OC.CPU || pkgversion(OC) >= v"0.96.22"
        # TODO: Add more diagnostics, make them dependent on simulation duration, take
        # monthly averages
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

    sim = OceananigansSimulation(ocean, area_fraction, ocean_properties, remapping)
    return sim
end

"""
    FieldExchanger.resolve_ocean_ice_fractions!(ocean_sim, ice_sim, land_fraction)

Ensure the ocean and ice area fractions are consistent with each other.
This matters in the case of a LatitudeLongitudeGrid, which is only
defined between -80 and 80 degrees latitude. In this case, we want to
set the ice fraction to `1 - land_fraction` on [-90, -80] and [80, 90]
degrees latitude, and make sure the ocean fraction is 0 there.
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

        # Set ice fraction to 1 - land_fraction and ocean fraction to 0 where polar_mask is 1
        @. ice_fraction = ifelse.(polar_mask == FT(1), FT(1) - land_fraction, ice_fraction)
        @. ocean_fraction = ifelse.(polar_mask == FT(1), FT(0), ocean_fraction)
    end
    return nothing
end

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

# Timestep the simulation forward to time `t`
Interfacer.step!(sim::OceananigansSimulation, t) =
    OC.time_step!(sim.ocean, float(t) - sim.ocean.model.clock.time)

# We always want the surface, so we always set zero(pt.lat) for z
"""
    to_node(pt::CA.ClimaCore.Geometry.LatLongPoint)

Transform `LatLongPoint` into a tuple (long, lat, 0), where the 0 is needed because we only
care about the surface.
"""
@inline to_node(pt::CA.ClimaCore.Geometry.LatLongPoint) = pt.long, pt.lat, zero(pt.lat)
# This next one is needed if we have "LevelGrid"
@inline to_node(pt::CA.ClimaCore.Geometry.LatLongZPoint) = pt.long, pt.lat, zero(pt.lat)

"""
    map_interpolate(points, oc_field::OC.Field)

Interpolate the given 3D field onto the target points.

If the underlying grid does not contain a given point, return 0 instead.

TODO: Use a non-allocating version of this function (simply replace `map` with `map!`)
"""
function map_interpolate(points, oc_field::OC.Field)
    loc = map(L -> L(), OC.Fields.location(oc_field))
    grid = oc_field.grid
    data = oc_field.data

    # TODO: There has to be a better way
    min_lat, max_lat = extrema(OC.φnodes(grid, OC.Center(), OC.Center(), OC.Center()))

    map(points) do pt
        FT = eltype(pt)

        # The oceananigans grid does not cover the entire globe, so we should not
        # interpolate outside of its latitude bounds. Instead we return 0
        min_lat < pt.lat < max_lat || return FT(0)

        fᵢ = OC.Fields.interpolate(to_node(pt), data, loc, grid)
        convert(FT, fᵢ)::FT
    end
end

"""
    surface_flux(f::OC.AbstractField)

Extract the top boundary conditions for the given field.
"""
function surface_flux(f::OC.AbstractField)
    top_bc = f.boundary_conditions.top
    if top_bc isa OC.BoundaryCondition{<:OC.BoundaryConditions.Flux}
        return top_bc.condition
    else
        return nothing
    end
end

function Interfacer.remap(field::OC.Field, target_space)
    return map_interpolate(CC.Fields.coordinate_field(target_space), field)
end

function Interfacer.remap(operation::OC.AbstractOperations.AbstractOperation, target_space)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    return Interfacer.remap(evaluated_field, target_space)
end

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
    Float32(0.06)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_diffuse_albedo}) =
    Float32(0.06)

# NOTE: This is 3D, but it will be remapped to 2D
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_temperature}) =
    273.15 + sim.ocean.model.tracers.T

"""
    FluxCalculator.update_turbulent_fluxes!(sim::OceananigansSimulation, fields)

Update the turbulent fluxes in the simulation using the values stored in the coupler fields.
These include latent heat flux, sensible heat flux, momentum fluxes, and moisture flux.

A note on sign conventions:
SurfaceFluxes and Oceananigans both use the convention that a positive flux is an upward flux.
No sign change is needed during the exchange, except for moisture/salinity fluxes:
SurfaceFluxes provides moisture moving from atmosphere to ocean as a negative flux at the surface,
and Oceananigans represents moisture moving from atmosphere to ocean as a positive salinity flux,
so a sign change is needed when we convert from moisture to salinity flux.
"""
function FluxCalculator.update_turbulent_fluxes!(sim::OceananigansSimulation, fields)
    # Only LatitudeLongitudeGrid are supported because otherwise we have to rotate the vectors

    (; F_lh, F_sh, F_turb_ρτxz, F_turb_ρτyz, F_turb_moisture) = fields
    grid = sim.ocean.model.grid

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
    oc_flux_u = surface_flux(sim.ocean.model.velocities.u)
    oc_flux_v = surface_flux(sim.ocean.model.velocities.v)
    set_from_extrinsic_vectors!(
        (; u = oc_flux_u, v = oc_flux_v),
        grid,
        F_turb_ρτxz_cc,
        F_turb_ρτyz_cc,
    )

    (; ocean_reference_density, ocean_heat_capacity, ocean_fresh_water_density) =
        sim.ocean_properties

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
        (remapped_F_lh .+ remapped_F_sh) ./ (ocean_reference_density * ocean_heat_capacity)

    # Add the part of the salinity flux that comes from the moisture flux, we also need to
    # add the component due to precipitation (that was done with the radiative fluxes)
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        F_turb_moisture,
    )
    moisture_fresh_water_flux = sim.remapping.scratch_arr1 ./ ocean_fresh_water_density
    oc_flux_S = surface_flux(sim.ocean.model.tracers.S)
    surface_salinity = OC.interior(sim.ocean.model.tracers.S, :, :, 1)
    OC.interior(oc_flux_S, :, :, 1) .=
        OC.interior(oc_flux_S, :, :, 1) .- surface_salinity .* moisture_fresh_water_flux
    return nothing
end

"""
    set_from_extrinsic_vectors!(vectors, grid, u_cc, v_cc)

Given the extrinsic vector components `u_cc` and `v_cc` as `Center, Center`
fields, rotate them onto the target grid and remap to `Face, Center` and
`Center, Face` fields, respectively.
"""
function set_from_extrinsic_vectors!(vectors, grid, u_cc, v_cc)
    arch = grid.architecture

    # Rotate vectors onto the grid
    OC.Utils.launch!(arch, grid, :xy, _rotate_velocities!, u_cc, v_cc, grid)

    # Fill halo regions with the rotated vectors so we can use them to interpolate
    OC.fill_halo_regions!(u_cc)
    OC.fill_halo_regions!(v_cc)

    # Interpolate the vectors to face/center and center/face respectively
    OC.Utils.launch!(
        arch,
        grid,
        :xy,
        _interpolate_velocities!,
        vectors.u,
        vectors.v,
        grid,
        u_cc,
        v_cc,
    )
    return nothing
end

"""
    _rotate_velocities!(u, v, grid)

Rotate the velocities from the extrinsic coordinate system to the intrinsic
coordinate system.
"""
@kernel function _rotate_velocities!(u, v, grid)
    # Use `k = 1` to index into the reduced Fields
    i, j = @index(Global, NTuple)
    # Rotate u, v from extrinsic to intrinsic coordinate system
    ur, vr = OC.Operators.intrinsic_vector(i, j, 1, grid, u, v)
    @inbounds begin
        u[i, j, 1] = ur
        v[i, j, 1] = vr
    end
end

"""
    _interpolate_velocities!(u, v, grid, u_cc, v_cc)

Interpolate the input velocities `u_cc` and `v_cc`, which are Center/Center
Fields to Face/Center and Center/Face coordinates, respectively.
"""
@kernel function _interpolate_velocities!(u, v, grid, u_cc, v_cc)
    # Use `k = 1` to index into the reduced Fields
    i, j = @index(Global, NTuple)
    @inbounds begin
        u[i, j, 1] = OC.Operators.ℑxyᶠᶜᵃ(i, j, 1, grid, u_cc)
        v[i, j, 1] = OC.Operators.ℑxyᶜᶠᵃ(i, j, 1, grid, v_cc)
    end
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

A note on sign conventions:
ClimaAtmos and Oceananigans both use the convention that a positive flux is an upward flux.
No sign change is needed during the exchange, except for precipitation/salinity fluxes.
ClimaAtmos provides precipitation as a negative flux at the surface, and
Oceananigans represents precipitation as a positive salinity flux,
so a sign change is needed when we convert from precipitation to salinity flux.
"""
function FieldExchanger.update_sim!(sim::OceananigansSimulation, csf)
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
    get_model_prog_state(sim::OceananigansSimulation)

Returns the model state of a simulation as a `ClimaCore.FieldVector`.
It's okay to leave this unimplemented for now, but we won't be able to use the
restart system.

TODO extend this for non-ClimaCore states.
"""
function Checkpointer.get_model_prog_state(sim::OceananigansSimulation)
    @warn "get_model_prog_state not implemented for OceananigansSimulation"
end
