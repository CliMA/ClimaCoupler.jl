import Oceananigans as OC
import ClimaOcean as CO
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaComms
import Thermodynamics as TD

"""
    OceananigansSimulation{SIM, A, OPROP, REMAP}

The ClimaCoupler simulation object used to run with Oceananigans.
This type is used by the coupler to indicate that this simulation
is an surface/ocean simulation for dispatch.

It contains the following objects:
- `sim::SIM`: The Oceananigans simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
- `ocean_properties::OPROP`: A NamedTuple of ocean properties and parameters
- `remapping::REMAP`: Objects needed to remap from the exchange (spectral) grid to Oceananigans spaces.
"""
struct OceananigansSimulation{SIM, A, OPROP, REMAP} <: Interfacer.OceanModelSimulation
    sim::SIM
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

Specific details about the complexity of the model
can be found in the Oceananigans.jl documentation.
"""
function OceananigansSimulation(area_fraction; output_dir, comms_ctx = ClimaComms.context())
    arch = comms_ctx.device isa ClimaComms.CUDADevice ? OC.GPU() : OC.CPU()

    # Set up ocean grid (1 degree)
    resolution_points = (360, 160, 32)
    Nz = last(resolution_points)
    z_faces = CO.exponential_z_faces(; Nz, depth = 6000, h = 34)

    # Regular LatLong because we know how to do interpolation there

    # TODO: When moving to TripolarGrid, note that we need to be careful about
    # ensuring the coordinate systems align (ie, rotate vectors on the OC grid)

    underlying_grid = OC.LatitudeLongitudeGrid(
        arch;
        size = resolution_points,
        longitude = (-180, 180),
        latitude = (-80, 80),   # NOTE: Don't goo to high up when using LatLongGrid, or the cells will be too small
        z = z_faces,
        halo = (7, 7, 7),
    )

    bottom_height =
        CO.regrid_bathymetry(underlying_grid; minimum_depth = 30, interpolation_passes = 20, major_basins = 1)

    grid = OC.ImmersedBoundaryGrid(underlying_grid, OC.GridFittedBottom(bottom_height); active_cells_map = true)

    # Create ocean simulation
    ocean = CO.ocean_simulation(grid)

    # Set up initial conditions for temperature and salinity
    T_init(λ, φ, z) = 30 * (1 - tanh((abs(φ) - 30) / 5)) / 2 + rand()
    S_init(λ, φ, z) = 30 - 5e-3 * z + rand()
    OC.set!(ocean.model, T = T_init, S = S_init)


    # TODO: Handle halos directly here instead of filling them later
    long_cc = OC.λnodes(grid, OC.Center(), OC.Center(), OC.Center())
    lat_cc = OC.φnodes(grid, OC.Center(), OC.Center(), OC.Center())

    # TODO: We can remove the `nothing` after CC > 0.14.33

    # TODO: Go from 0 to Nx+1, Ny+1 (for halos) (for LatLongGrid)

    # Construct three remappers, one to Center, Center fields, one
    # to Center, Face fields, and one to Face, Center fields.
    long_cc = reshape(long_cc, length(long_cc), 1)
    lat_cc = reshape(lat_cc, 1, length(lat_cc))
    target_points_cc = @. CC.Geometry.LatLongPoint(lat_cc, long_cc)
    remapper_cc = CC.Remapping.Remapper(axes(area_fraction), target_points_cc, nothing)

    # TODO: Get rid of this and remap to cc
    long_cf = OC.λnodes(grid, OC.Center(), OC.Face(), OC.Center())
    lat_cf = OC.φnodes(grid, OC.Center(), OC.Face(), OC.Center())

    long_cf = reshape(long_cf, length(long_cf), 1)
    lat_cf = reshape(lat_cf, 1, length(lat_cf))
    target_points_cf = @. CC.Geometry.LatLongPoint(lat_cf, long_cf)
    remapper_cf = CC.Remapping.Remapper(axes(area_fraction), target_points_cf, nothing)

    long_fc = OC.λnodes(grid, OC.Face(), OC.Center(), OC.Center())
    lat_fc = OC.φnodes(grid, OC.Face(), OC.Center(), OC.Center())

    long_fc = reshape(long_fc, length(long_fc), 1)
    lat_fc = reshape(lat_fc, 1, length(lat_fc))
    target_points_fc = @. CC.Geometry.LatLongPoint(lat_fc, long_fc)
    remapper_fc = CC.Remapping.Remapper(axes(area_fraction), target_points_fc, nothing)

    ocean_properties = (; ocean_reference_density = 1020, ocean_heat_capacity = 3991, ocean_fresh_water_density = 999.8)

    remapping = (; remapper_cc, remapper_fc, remapper_cf)

    # Before version 0.96.22, the NetCDFWriter was broken on GPU
    if arch isa OC.CPU || pkgversion(OC) >= v"0.96.22"
        # TODO: Add more diagnostics, make them dependent on simulation duration, take
        # monthly averages
        diagnostics = Dict("T" => ocean.model.tracers.T)
        netcdf_writer = OC.NetCDFWriter(
            ocean.model,
            diagnostics,
            indices = (:, :, grid.Nz),
            filename = joinpath(output_dir, "ocean_diagnostics.nc"),
            schedule = OC.TimeInterval(86400),
        )
        ocean.output_writers[:diagnostics] = netcdf_writer
    end

    sim = OceananigansSimulation(ocean, area_fraction, ocean_properties, remapping)
    return sim
end

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

# Timestep the simulation forward to time `t`
Interfacer.step!(sim::OceananigansSimulation, t) = OC.time_step!(sim.sim, float(t) - sim.sim.model.clock.time)

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
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_buoyancy}) = Float32(5.8e-5)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_momentum}) = Float32(5.8e-5)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:beta}) = Float32(1)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_direct_albedo}) = Float32(0.06)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_diffuse_albedo}) = Float32(0.06)

# NOTE: This is 3D, but it will be remapped to 2D
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_temperature}) = 273.15 + sim.sim.model.tracers.T

function FluxCalculator.update_turbulent_fluxes!(sim::OceananigansSimulation, fields)
    # Only LatitudeLongitudeGrid are supported because otherwise we have to rotate the vectors

    (; F_lh, F_sh, F_turb_ρτxz, F_turb_ρτyz, F_turb_moisture) = fields

    # TODO: Starting ClimaCore > 0.14.33 we can use interpolate!
    remapped_F_turb_ρτxz = CC.Remapping.interpolate(sim.remapping.remapper_fc, F_turb_ρτxz)
    remapped_F_turb_ρτyz = CC.Remapping.interpolate(sim.remapping.remapper_cf, F_turb_ρτyz)

    # TODO: Remap to Center Center. Instead of dealing with separate fields, we can remap everything
    # onto Center, Center fields on the Oceananingans side, and then remap from ClimaCoupler to this
    # See also comment below next to the function

    # set_from_extrinsic_vectors!((; u = oc_flux_u, v = oc_flux_v),
    #                             sim.model.grid,
    #                             remapped_F_turb_ρτxz,
    #                             remapped_F_turb_ρτyz
    #                             )

    oc_flux_u = surface_flux(sim.sim.model.velocities.u)
    oc_flux_v = surface_flux(sim.sim.model.velocities.v)

    view(oc_flux_u, :, :, 1) .= remapped_F_turb_ρτxz
    view(oc_flux_v, :, :, 1) .= remapped_F_turb_ρτyz

    # TODO: Handle this directly when interpolating
    OC.fill_halo_regions!(oc_flux_u)
    OC.fill_halo_regions!(oc_flux_v)

    (; ocean_reference_density, ocean_heat_capacity, ocean_fresh_water_density) = sim.ocean_properties
    remapped_F_lh = CC.Remapping.interpolate(sim.remapping.remapper_cc, F_lh)
    remapped_F_sh = CC.Remapping.interpolate(sim.remapping.remapper_cc, F_sh)
    remapped_F_turb_energy = remapped_F_lh + remapped_F_sh

    # TODO: Note, SW radiation penetrates the surface. Right now, we just put
    # everything on the surface, but later we will need to account for this.
    # One way we can do this is using directly ClimaOcean

    oc_flux_T = surface_flux(sim.sim.model.tracers.T)
    view(oc_flux_T, :, :, 1) .=
        OC.interior(oc_flux_T, :, :, 1) .+ remapped_F_turb_energy ./ (ocean_reference_density * ocean_heat_capacity)
    OC.fill_halo_regions!(oc_flux_T)

    # Add the part of the salinity flux that comes from the moisture flux, we also need to
    # add the component due to precipitation (that was done with the radiative fluxes)
    remapped_F_turb_moisture = CC.Remapping.interpolate(sim.remapping.remapper_cc, F_turb_moisture)
    oc_flux_S = surface_flux(sim.sim.model.tracers.S)
    surface_salinity = OC.interior(sim.sim.model.tracers.S, :, :, 1)
    moisture_fresh_water_flux = remapped_F_turb_moisture ./ ocean_fresh_water_density
    view(oc_flux_S, :, :, 1) .= OC.interior(oc_flux_S, :, :, 1) .- surface_salinity .* moisture_fresh_water_flux

    OC.fill_halo_regions!(oc_flux_S)
    return nothing
end

# TODO: This commented out function works with 3D fields and allocates fields internally, we
# need preallocate the fields and make it compatible with 2D.
# function set_from_extrinsic_vectors!(vectors, grid, u, v)
#     grid = grid
#     arch = grid.architecture
#     # TODO: Move this outside
#     uᶜᶜᶜ = CenterField(grid)
#     vᶜᶜᶜ = CenterField(grid)
#     set!(uᶜᶜᶜ, u)
#     set!(vᶜᶜᶜ, v)

#     # TODO: Change these kernels to be 2D
#     launch!(arch, grid, :xyz, _rotate_vectors!, uᶜᶜᶜ, vᶜᶜᶜ, grid)
#     launch!(arch, grid, :xyz, _interpolate_vectors!,
#             vectors.u, vectors.v, grid, uᶜᶜᶜ, vᶜᶜᶜ)
#     return nothing
# end

function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:area_fraction}, field)
    sim.area_fraction .= field
    return nothing
end

"""
    FieldExchanger.update_sim!(sim::OceananigansSimulation, csf, area_fraction)

Update the ocean simulation with the provided fields, which have been filled in
by the coupler.

Update the portion of the surface_fluxes for T and S that is due to radiation and
precipitation. The rest will be updated in `update_turbulent_fluxes!`.
"""
function FieldExchanger.update_sim!(sim::OceananigansSimulation, csf, area_fraction)
    (; ocean_reference_density, ocean_heat_capacity, ocean_fresh_water_density) = sim.ocean_properties

    remapped_F_radiative = CC.Remapping.interpolate(sim.remapping.remapper_cc, csf.F_radiative)

    # Update only the part due to radiative fluxes. For the full update, the component due
    # to latent and sensible heat is missing and will be updated in update_turbulent_fluxes.
    oc_flux_T = surface_flux(sim.sim.model.tracers.T)
    oc_flux_T .= remapped_F_radiative ./ (ocean_reference_density * ocean_heat_capacity)

    remapped_P_liq = CC.Remapping.interpolate(sim.remapping.remapper_cc, csf.P_liq)
    remapped_P_snow = CC.Remapping.interpolate(sim.remapping.remapper_cc, csf.P_snow)

    # Virtual salt flux
    oc_flux_S = surface_flux(sim.sim.model.tracers.S)
    precipitating_fresh_water_flux = (remapped_P_liq .+ remapped_P_snow) ./ ocean_fresh_water_density
    surface_salinity_flux = OC.interior(sim.sim.model.tracers.S, :, :, 1) .* precipitating_fresh_water_flux
    view(oc_flux_S, :, :, 1) .= .-surface_salinity_flux
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
