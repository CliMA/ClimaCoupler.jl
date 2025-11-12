import Oceananigans as OC
import ClimaOcean as CO
import ClimaAtmos as CA
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaComms
import ClimaCore as CC
import Thermodynamics as TD
import ClimaOcean.EN4: download_dataset
using KernelAbstractions: @kernel, @index, @inbounds
using XESMF # to load Oceananigans regridding extension

OceananigansXESMFExt = Base.get_extension(OC, :OceananigansXESMFExt).OceananigansXESMFExt;

"""
    OceananigansSimulation{SIM, A, OPROP, REMAP}

The ClimaCoupler simulation object used to run with Oceananigans.
This type is used by the coupler to indicate that this simulation
is an surface/ocean simulation for dispatch.

It contains the following objects:
- `ocean::SIM`: The Oceananigans simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
- `ocean_properties::OPROP`: A NamedTuple of ocean properties and parameters
- `remapper_oc_to_cc::REMAP1`: Objects needed to remap from the Oceananigans space to the exchange (spectral) grid.
- `remapper_cc_to_oc::REMAP2`: Objects needed to remap from the exchange (spectral) grid to the Oceananigans space.
"""
struct OceananigansSimulation{SIM, A, OPROP, REMAP1, REMAP2} <:
       Interfacer.OceanModelSimulation
    ocean::SIM
    area_fraction::A
    ocean_properties::OPROP
    remapper_oc_to_cc::REMAP1
    remapper_cc_to_oc::REMAP2
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

    # Get the remapper objects to go in both directions between the Oceananigans and Cubed sphere grids
    boundary_space = axes(area_fraction)
    remapper_oc_to_cc, remapper_cc_to_oc = construct_remappers(grid, boundary_space)

    ocean_properties = (;
        ocean_reference_density = 1020,
        ocean_heat_capacity = 3991,
        ocean_fresh_water_density = 999.8,
    )

    # Save all tracers and velocities to a JLD2 file at daily frequency
    outputs = merge(ocean.model.tracers, ocean.model.velocities)
    jld2_writer = OC.JLD2Writer(
        ocean.model,
        outputs;
        schedule = OC.TimeInterval(86400), # Daily output
        filename = joinpath(output_dir, "ocean_diagnostics"),
        indices = (:, :, grid.Nz),
        overwrite_existing = true,
        array_type = Array{Float32},
    )
    ocean.output_writers[:diagnostics] = jld2_writer

    sim = OceananigansSimulation(
        ocean,
        area_fraction,
        ocean_properties,
        remapper_oc_to_cc,
        remapper_cc_to_oc,
    )
    return sim
end

"""
    construct_remappers(grid, boundary_space)

Construct the remapper objects to go in both directions between the Oceananigans and Cubed sphere grids.
Both objects contain a remapper object and relevant scratch space.

- Oceananigans to ClimaCore
In this direction we use XESMF bilinear interpolation.
Note that we assume the Oceananigans Field is on Center, Center, Center.

    Example: remap the Oceananigans field `T` to the ClimaCore field `T_climacore`

    # Convert the Oceananigans field to a flat vector
    remapper_oc_to_cc.src_vec .= Array(vec(OC.interior(T, :, :, Nz))) # Oceananigans source vector

    # Apply the XESMF regridder (matrix multiply)
    remapper_oc_to_cc.remapper(remapper_oc_to_cc.dest_vec, remapper_oc_to_cc.src_vec)

    # Convert the output vector to a 2D array
    remapper_oc_to_cc.dest_arr .= reshape(remapper_oc_to_cc.dest_vec, size(remapper_oc_to_cc.dest_arr)...)

    # Copy the remapped data to the ClimaCore field
    # Note: in general we avoid accessing the parent of a field, but here we make an exception
    # since we're using the underlying array to remap.
    parent(T_climacore) .= remapper_oc_to_cc.dest_arr

- ClimaCore to Oceananigans
In this direction we use the ClimaCore remapper for this because it uses all the information
we have about the spectral element cubed sphere grid, which XESMF does not support.

    Example: remap the ClimaCore field `F_turb_ρτxz_cc` to the Oceananigans field `F_turb_ρτxz_oc`:

    # Remap the ClimaCore momentum flux to a Oceananigans field using scratch arrays and fields
    CC.Remapping.interpolate!(
        remapper_cc_to_oc.scratch_arr1,
        remapper_cc_to_oc.remapper,
        F_turb_ρτxz_cc, # ClimaCore field
    )
    OC.set!(remapper_cc_to_oc.scratch_oc1, remapper_cc_to_oc.scratch_arr1) # zonal momentum flux
    F_turb_ρτxz_oc = remapper_cc_to_oc.scratch_oc1 # Oceananigans field

Arguments:
- `grid`: The Oceananigans grid (TripolarGrid or LatitudeLongitudeGrid).
- `boundary_space`: The boundary space (ClimaCore SpectralElementSpace2D).

Returns:
- `remapper_oc_to_cc`: The remapper object to go from the Oceananigans grid to the Cubed sphere nodes.
"""
function construct_remappers(grid, boundary_space)
    ## Remapper: Oceananigans `Center, Center` to Cubed sphere nodes
    # Get the Oceananigans coordinates and put them on CPU
    coords_oc =
        OceananigansXESMFExt.xesmf_coordinates(grid, OC.Center(), OC.Center(), OC.Center())
    coords_oc = Dict(k => Array(v) for (k, v) in coords_oc)
    if grid.underlying_grid isa OC.TripolarGrid
        # TripolarGrid is defined on [0, 360], so we to convert to [-180, 180] to match ClimaCore
        coords_oc["lon"] .-= 180
    end

    # Get the latitude and longitude of each node on the boundary space
    climacore_coords = CC.Fields.coordinate_field(boundary_space)

    # Get the cubed sphere latitude and longitude, each as an Nx1 Matrix
    climacore_lat = Array(reshape(vec(parent(climacore_coords.lat)), :, 1))
    climacore_lon = Array(reshape(vec(parent(climacore_coords.long)), :, 1))

    coords_climacore = Dict("lat" => climacore_lat, "lon" => climacore_lon)

    # Construct the XESMF regridder object
    regridder_oceananigans_to_climacore =
        XESMF.Regridder(coords_oc, coords_climacore; method = "bilinear")

    # Allocate space for source an destination vectors to use as intermediate storage
    src_vec_oc = Array(vec(OC.Field{OC.Center, OC.Center, Nothing}(grid))) # 2D field on Center/Center
    field_climacore = CC.Fields.zeros(boundary_space) # 2D field on boundary space (cubed sphere)
    dest_vec_climacore = vec(parent(field_climacore))
    dest_arr_climacore =
        deepcopy(reshape(dest_vec_climacore, size(parent(field_climacore))...))

    remapper_oc_to_cc = (;
        remapper = regridder_oceananigans_to_climacore,
        src_vec = src_vec_oc,
        dest_vec = dest_vec_climacore,
        dest_arr = dest_arr_climacore,
    )

    ## Remapper: Cubed sphere nodes to Oceananigans grid `Center, Center`
    # For a TripolarGrid, latitude and longitude are already 2D arrays,
    # so we broadcast over them directly.
    long_oc = OC.λnodes(grid, OC.Center(), OC.Center(), OC.Center())
    lat_oc = OC.φnodes(grid, OC.Center(), OC.Center(), OC.Center())

    # For a LatitudeLongitudeGrid, latitude and longitude are 1D collections,
    # so we reshape them to 2D arrays.
    if grid.underlying_grid isa OC.LatitudeLongitudeGrid
        long_oc = reshape(long_oc, length(long_oc), 1)
        lat_oc = reshape(lat_oc, 1, length(lat_oc))
    else
        # TripolarGrid is defined on [0, 360], so we to convert to [-180, 180] to match ClimaCore
        long_oc .-= 180
    end
    target_points_oc = @. CC.Geometry.LatLongPoint(lat_oc, long_oc)

    # Construct the ClimaCore remapper object
    remapper_climacore_to_oceananigans =
        CC.Remapping.Remapper(boundary_space, target_points_oc)

    # Construct two 2D Center/Center fields to use as scratch space while remapping
    scratch_oc1 = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    scratch_oc2 = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    # Construct two scratch arrays to use while remapping
    # We get the array type, float type, and dimensions from the remapper object to maintain consistency
    ArrayType = ClimaComms.array_type(boundary_space)
    FT = CC.Spaces.undertype(boundary_space)
    interpolated_values_dim..., _buffer_length =
        size(remapper_climacore_to_oceananigans._interpolated_values)
    scratch_arr1 = ArrayType(zeros(FT, interpolated_values_dim...))
    scratch_arr2 = ArrayType(zeros(FT, interpolated_values_dim...))

    remapper_cc_to_oc = (;
        remapper = remapper_climacore_to_oceananigans,
        scratch_oc1,
        scratch_oc2,
        scratch_arr1,
        scratch_arr2,
    )

    return remapper_oc_to_cc, remapper_cc_to_oc
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

"""
    Interfacer.remap(field::OC.Field, target_space, remapper_oc_to_cc)
    Interfacer.remap(operation::OC.AbstractOperation, target_space, remapper_oc_to_cc)

Remap the given Oceananigans field onto the target space using the remapper object.
If an operation is provided, it is evaluated and the resulting field is remapped.

Arguments:
- `field/operation`: The Oceananigans field or operation to remap.
- `target_space`: The target space (ClimaCore SpectralElementSpace2D).
- `remapper_oc_to_cc`: The remapper object to go from the Oceananigans grid to cubed sphere nodes.
"""
function Interfacer.remap(field::OC.Field, target_space, remapper_oc_to_cc)
    # Allocate a new ClimaCore field and remap the data to it
    target_field = CC.Fields.zeros(target_space)
    Interfacer.remap!(target_field, field, target_space, remapper_oc_to_cc)
    return target_field
end

function Interfacer.remap(operation::OC.AbstractOperations.AbstractOperation, target_space)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    return Interfacer.remap(evaluated_field, target_space)
end

function Interfacer.remap!(
    target_field::CC.Fields.Field,
    field::OC.Field,
    remapper_oc_to_cc,
)
    # Since we always regrid from the surface of the ocean, take only the top layer of the field
    Nz = field.grid.underlying_grid.Nz

    # Convert the Oceananigans field to a flat vector
    remapper_oc_to_cc.src_vec .= Array(vec(OC.interior(field, :, :, Nz)))

    # Apply the XESMF regridder (matrix multiply)
    remapper_oc_to_cc.remapper(remapper_oc_to_cc.dest_vec, remapper_oc_to_cc.src_vec)

    # Convert the output vector to a 2D array
    remapper_oc_to_cc.dest_arr .=
        reshape(remapper_oc_to_cc.dest_vec, size(remapper_oc_to_cc.dest_arr)...)

    # Allocate a new ClimaCore field and copy the remapped data to it
    # Note: in general we avoid accessing the parent of a field, but here we make an exception
    # since we're already remapping with the underlying array.
    parent(target_field) .= remapper_oc_to_cc.dest_arr
    return nothing
end

function Interfacer.remap!(
    target_field::CC.Fields.Field,
    operation::OC.AbstractOperations.AbstractOperation,
    target_space,
    remapper_oc_to_cc,
)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    return Interfacer.remap!(target_field, evaluated_field, target_space, remapper_oc_to_cc)
end

"""
    Interfacer.get_remapper_to_cc(sim::OceananigansSimulation)

Return the remapper object used to remap quantities from the Oceananigans grid
to the ClimaCore boundary space.
"""
Interfacer.get_remapper_to_cc(sim::OceananigansSimulation) = sim.remapper_oc_to_cc

# TODO: Better values for this

# At the moment, we return always Float32. This is because we always want to run
# Oceananingans with Float64, so we have no way to know the float type here. Sticking with
# Float32 ensures that nothing is accidentally promoted to Float64. We will need to change
# this anyway.
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:area_fraction}) = sim.area_fraction
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
        sim.remapping.remapper_cc_to_oc,
        F_turb_ρτxz,
    )
    OC.set!(sim.remapping.scratch_oc1, sim.remapping.scratch_arr1) # zonal momentum flux
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr2,
        sim.remapping.remapper_cc_to_oc,
        F_turb_ρτyz,
    )
    OC.set!(sim.remapping.scratch_oc2, sim.remapping.scratch_arr2) # meridional momentum flux

    # Rename for clarity; these are now Center, Center Oceananigans fields
    F_turb_ρτxz_cc = sim.remapping.scratch_oc1
    F_turb_ρτyz_cc = sim.remapping.scratch_oc2

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
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc_to_oc,
        F_lh,
    ) # latent heat flux
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr2,
        sim.remapping.remapper_cc_to_oc,
        F_sh,
    ) # sensible heat flux

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
        sim.remapping.remapper_cc_to_oc,
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
        sim.remapping.remapper_cc_to_oc,
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
        sim.remapping.remapper_cc_to_oc,
        csf.P_liq,
    )
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr2,
        sim.remapping.remapper_cc_to_oc,
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
