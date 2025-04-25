import Oceananigans as OC
import ClimaOcean as CO
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaComms
import Thermodynamics as TD

const KA = OC.Architectures.KernelAbstractions

import ClimaOcean.OceanSeaIceModels.PrescribedAtmospheres:
    thermodynamics_parameters,
    boundary_layer_height,
    surface_layer_height

"""
    OceananigansSimulation{SIM, A}

The ClimaCoupler simulation object used to run with Oceananigans.
This type is used by the coupler to indicate that this simulation
is an surface/ocean simulation for dispatch.

It contains the following objects:
- `sim::SIM`: The Oceananigans simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
"""
struct OceananigansSimulation{SIM, A} <: Interfacer.OceanModelSimulation
    sim::SIM
    area_fraction::A
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
function OceananigansSimulation(
    atmosphere::ClimaAtmosSimulation,
    area_fraction;
    resolution_points = (90, 45, 15),
    comms_ctx = ClimaComms.context(),
)
    # TODO fill this out
    arch = if comms_ctx.device isa ClimaComms.CUDADevice
        OC.GPU()
    else
        OC.CPU()
    end

    # Set up ocean grid (1 degree)
    Nx, Ny, Nz = resolution_points
    # Nx = 256
    # Ny = 128
    # Nz = 32
    z_faces = CO.exponential_z_faces(; Nz, depth = 6000, h = 34)
    # underlying_grid = TripolarGrid(arch; size = (Nx, Ny, Nz), z = z_faces)

    grid = OC.LatitudeLongitudeGrid(arch;
                                              size = (Nx, Ny, Nz),
                                              longitude = (0, 360),
                                              latitude = (-85, 85),
                                              z = z_faces,
                                              halo = (7, 7, 7))

    # bottom_height = regrid_bathymetry(underlying_grid; minimum_depth = 30, interpolation_passes = 20, major_basins = 1)
    # view(bottom_height, 73:78, 88:89, 1) .= -1000 # open Gibraltar strait

    # TODO how to specify FT?
    # grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height); active_cells_map = true)

    # Create ocean simulation
    ocean = CO.ocean_simulation(grid)

    # TODO: Remove this. It is here right now to add interfaces to Oceananigans
    sea_ice = CO.OceanSeaIceModels.FreezingLimitedOceanTemperature(eltype(ocean.model))

    ocean_sea_ice = CO.OceanSeaIceModel(ocean, sea_ice; atmosphere)

    # Set up initial conditions for temperature and salinity
    Tᵢ(λ, φ, z) = 300
    #273 + 30 * (1 - tanh((abs(φ) - 30) / 5)) / 2 + rand()
    Sᵢ(λ, φ, z) = 30 - 5e-3 * z + rand()
    OC.set!(ocean.model, T = Tᵢ, S = Sᵢ)

    return OceananigansSimulation(ocean_sea_ice, area_fraction)
end

Interfacer.name(::OceananigansSimulation) = "OceananigansSimulation"

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

# Timestep the simulation forward to time `t`
Interfacer.step!(sim::OceananigansSimulation, t) = OC.time_step!(sim.sim, float(t) - sim.sim.clock.time)

# Reset prognostic state and current time to initial conditions
Interfacer.reinit!(sim::OceananigansSimulation) = nothing # TODO fill this out

function OC.time_step!(atmos::ClimaAtmosSimulation, Δt)
    # It is not Oceananigans job to step the model
    return nothing
end

@inline to_node(pt::CA.ClimaCore.Geometry.LatLongPoint) = pt.long, pt.lat, zero(pt.lat)

@inline to_node(pt::CA.ClimaCore.Geometry.LatLongZPoint) = pt.long, pt.lat, pt.z

instantiate(L) = L()

function map_interpolate(points, oc_field::OC.Field) #, loc, grid)
    loc = map(instantiate, OC.Fields.location(oc_field))
    grid = oc_field.grid
    data = oc_field.data

    map(points) do pt
        FT = eltype(pt)
        fᵢ = OC.Fields.interpolate(to_node(pt), data, loc, grid)
        convert(FT, fᵢ)
    end
end

function Interfacer.remap(field::OC.Field, target_space)
    return map_interpolate(CC.Fields.coordinate_field(target_space), field)
end

# """
#     Interfacer.get_field(sim::OceananigansSimulation, ::Val{:_})

# Get the value of the specified field in the ocean simulation.
# This is used to exchange properties from the ocean to the coupler and other components.
# """
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:area_fraction}) = sim.area_fraction

# FIXME! Don't hardcode this
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_buoyancy}) = Float32(1e-3)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_momentum}) = Float32(1e-3)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:beta}) = Float32(0.1)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_direct_albedo}) = Float32(0.38)
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_diffuse_albedo}) = Float32(0.38)
# Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_humidity}) = return nothing # TODO fill this out

# This is 3D, but it will be remapped to 2D
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_temperature}) = sim.sim.ocean.model.tracers.T
#

function FluxCalculator.update_turbulent_fluxes!(sim::OceananigansSimulation, fields)
    # Left for later
    return nothing
end


# # These two methods are used to track conservation of the coupled system, so it's okay to leave them empty for now
# function Interfacer.get_field(sim::OceananigansSimulation, ::Val{:energy}, level)
#     # TODO fill this out
#     return nothing
# end
# function Interfacer.get_field(sim::OceananigansSimulation, ::Val{:water}, level)
#     # TODO fill this out
#     return nothing
# end

# """
#     Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:_}, field)

# Update the specified field in the ocean simulation with the provided field, which
# is a ClimaCore Field defined on the boundary space. Remapping will be needed here.
# This is used to exchange properties from the coupler and other components to the ocean.
# """
# # Setters to update fields stored in the ocean simulation with the provided field (a ClimaCore Field defined on boundary_space)
# # These are just the `update_field!` methods corresponding to the default coupler fields,
# #  but if the ocean requires more fields, we can add them here. (maybe air temperature, pressure, humidity, CO2?)
# # If any of these aren't needed by the ocean model, we can remove them.
# function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:air_density}, field)
#     # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
# end
# function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:area_fraction}, field)
#     parent(sim.area_fraction) .= parent(field)
# end
# function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:liquid_precipitation}, field)
#     # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
# end
# function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:snow_precipitation}, field)
#     # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
# end
# function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:lw_d}, field)
#     # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
# end
# function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:sw_d}, field)
#     # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
# end
# function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:turbulent_energy_flux}, field)
#     # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
#     #  or we can update the model directly in the flux calculation, depending on where we compute fluxes
# end
# function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:turbulent_moisture_flux}, field)
#     # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
#     #  or we can update the model directly in the flux calculation, depending on where we compute fluxes
# end
# function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:turbulent_momentum_flux_x}, field)
#     # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
#     #  or we can update the model directly in the flux calculation, depending on where we compute fluxes
# end
# function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:turbulent_momentum_flux_y}, field)
#     # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
#     #  or we can update the model directly in the flux calculation, depending on where we compute fluxes
# end

# """
#     FieldExchanger.update_sim!(sim::OceananigansSimulation, csf, turbulent_fluxes, area_fraction)

# Update the ocean simulation with the provided fields, which have been filled in
# by the coupler.
# """
# function FieldExchanger.update_sim!(sim::OceananigansSimulation, csf, turbulent_fluxes, area_fraction)
#     Interfacer.update_field!(sim, Val(:air_density), air_density)
#     Interfacer.update_field!(sim, Val(:area_fraction), area_fraction)

#     # precipitation
#     Interfacer.update_field!(sim, Val(:liquid_precipitation), csf.P_liq)
#     Interfacer.update_field!(sim, Val(:snow_precipitation), csf.P_snow)

#     # update fields for radiative transfer
#     Interfacer.update_field!(sim, Val(:sw_d), csf.SW_d)
#     Interfacer.update_field!(sim, Val(:lw_d), csf.LW_d)

#     # TODO update other fields as needed
# end

# """
#     update_turbulent_fluxes!(sim::OceananigansSimulation, fields::NamedTuple)

# Update the turbulent fluxes in the ocean simulation with the provided fields,
# which were filled in during the coupler's flux calculation.

# TODO we don't need this function if we update the model directly in the flux calculation
# """
# function FluxCalculator.update_turbulent_fluxes!(sim::OceananigansSimulation, fields::NamedTuple)
#     Interfacer.update_field!(sim, Val(:F_turb_ρτxz), fields.F_turb_ρτxz)
#     Interfacer.update_field!(sim, Val(:F_turb_ρτyz), fields.F_turb_ρτyz)
#     Interfacer.update_field!(sim, Val(:F_turb_energy), fields.F_turb_energy)
#     Interfacer.update_field!(sim, Val(:F_turb_moisture), fields.F_turb_moisture)
#     return nothing
# end

# """
#     Interfacer.add_coupler_fields!(coupler_field_names, ::OceananigansSimulation)

# Extend Interfacer.add_coupler_fields! to add the fields required for
# OceananigansSimulation, if any.
# """
# function Interfacer.add_coupler_fields!(coupler_field_names, ::OceananigansSimulation)
#     # TODO check Interfacer docs for default exchange fields, fill this out if needed
#     ocean_coupler_fields = [:SW_d, :LW_d]
#     push!(coupler_field_names, ocean_coupler_fields...)
# end

# """
#     FieldExchanger.import_atmos_fields!(csf, sim::OceananigansSimulation, atmos_sim, turbulent_fluxes)

# Update the coupler fields in-place with the values from the atmosphere simulation.
# This is defined here because the coupler exchange fields are specified by the
# Oceananigans simulation `add_coupler_fields!` method, and we need to know which fields
# are required from the atmosphere model.
# """
# function FieldExchanger.import_atmos_fields!(csf, sim::OceananigansSimulation, atmos_sim, turbulent_fluxes)
#     # TODO fix remap function calls

#     # radiative fluxes
#     remap!(csf.SW_d, Interfacer.get_field(atmos_sim, Val(:sw_d)))
#     remap!(csf.LW_d, Interfacer.get_field(atmos_sim, Val(:lw_d)))
#     # precipitation
#     remap!(csf.P_liq, Interfacer.get_field(atmos_sim, Val(:liquid_precipitation)))
#     remap!(csf.P_snow, Interfacer.get_field(atmos_sim, Val(:snow_precipitation)))
#     # air density
#     remap!(csf.air_density, Interfacer.get_field(atmos_sim, Val(:air_density)))

#     # TODO import other fields as needed (according to how we extend `add_coupler_fields!`)
# end

# ## Extend functions for ocean-specific flux calculation
# """
#     compute_surface_fluxes!(csf, sim::OceananigansSimulation, atmos_sim, boundary_space, thermo_params, surface_scheme)

# This function computes surface fluxes between the Oceananigans simulation and
# the atmosphere.

# Update the input coupler surface fields `csf` in-place with the computed fluxes
# for this model. These are then summed using area-weighting across all surface
# models to get the total fluxes.

# # Arguments
# - `csf`: [CC.Fields.Field] containing a NamedTuple of turbulent flux fields: `F_turb_ρτxz`, `F_turb_ρτyz`, `F_turb_energy`, `F_turb_moisture`.
# - `sim`: [OceananigansSimulation] the ocean simulation to compute fluxes for.
# - `atmos_sim`: [Interfacer.AtmosModelSimulation] the atmosphere simulation to compute fluxes with.
# - unused arguments: `boundary_space`, `thermo_params`, `surface_scheme`
# """
# function FluxCalculator.compute_surface_fluxes!(
#     csf,
#     sim::OceananigansSimulation,
#     atmos_sim::Interfacer.AtmosModelSimulation,
#     _...,
# )
#     # TODO fill this in -
#     #  compute F_turb_ρτxz, F_turb_ρτyz, F_shf, F_lhf, F_turb_moisture using Oceananigans.jl functions
#     #  something like this:

#     fluxes = Oceananigans.compute_surface_fluxes(p, sim.model, Y, t, atmos_sim.integrator)
#     (; F_turb_ρτxz, F_turb_ρτyz, F_shf, F_lhf, F_turb_moisture) = fluxes

#     # get area fraction (min = 0, max = 1)
#     area_fraction = Interfacer.get_field(sim, Val(:area_fraction))

#     # add the flux contributing from this surface to the coupler field
#     # note that the fluxes are area-weighted, so if a surface model is
#     #  not present at this point, the fluxes are zero
#     @. csf.F_turb_ρτxz += F_turb_ρτxz * area_fraction
#     @. csf.F_turb_ρτyz += F_turb_ρτyz * area_fraction
#     @. csf.F_turb_energy += (F_shf .+ F_lhf) * area_fraction
#     @. csf.F_turb_moisture += F_turb_moisture * area_fraction
#     return nothing
# end

"""
    get_model_prog_state(sim::OceananigansSimulation)

Returns the model state of a simulation as a `ClimaCore.FieldVector`.
It's okay to leave this unimplemented for now, but we won't be able to use the
restart system.

TODO extend this for non-ClimaCore states.
"""
function Checkpointer.get_model_prog_state(sim::OceananigansSimulation)
    error("get_model_prog_state not implemented")
end



import ClimaOcean.OceanSeaIceModels.InterfaceComputations:
    atmosphere_exchanger,
    initialize!,
    StateExchanger,
    interpolate_atmosphere_state!

import ClimaOcean.OceanSeaIceModels.PrescribedAtmospheres:
    thermodynamics_parameters,
    boundary_layer_height,
    surface_layer_height

# The height of near-surface variables used in the turbulent flux solver
surface_layer_height(s::ClimaAtmosSimulation) = 10 # meters, for example

# This is a parameter that is used in the computation of the fluxes,
# It probably should not be here but in the similarity theory type.
boundary_layer_height(atmos::ClimaAtmosSimulation) = 600

# Note: possibly, can use the atmos thermodynamic parameters directly here.
thermodynamics_parameters(atmos::ClimaAtmosSimulation) =
    atmos.integrator.p.params.thermodynamics_params

"""
    interpolate_atmospheric_state!(surface_atmosphere_state,
                                        interpolated_prescribed_freshwater_flux,
                                        atmos::ClimaAtmosSimulation,
                                        grid, clock)

Interpolate the atmospheric state in `atmos` to `surface_atmospheric_state`, a
the collection of `Field`s needed to compute turbulent fluxes.
"""
function interpolate_atmosphere_state!(interfaces,
                                       atmosphere::ClimaAtmosSimulation,
                                       coupled_model)

    interpolator = interfaces.exchanger.atmosphere_exchanger.to_exchange_interp
    exchange_atmosphere_state = interfaces.exchanger.exchange_atmosphere_state

    ue  = parent(exchange_atmosphere_state.u)
    ve  = parent(exchange_atmosphere_state.v)
    Te  = parent(exchange_atmosphere_state.T)
    qe  = parent(exchange_atmosphere_state.q)
    pe  = parent(exchange_atmosphere_state.p)

    ue = dropdims(ue, dims=3)
    ve = dropdims(ve, dims=3)
    Te = dropdims(Te, dims=3)
    qe = dropdims(qe, dims=3)
    pe = dropdims(pe, dims=3)

    Uah = CC.Geometry.UVVector.(CC.Spaces.level(atmosphere.integrator.u.c.uₕ, 1))
    ua = Uah.components.data.:1
    va = Uah.components.data.:2

    # TODO: can we avoid allocating for Ta, pa, qa?
    tsa = CC.Spaces.level(atmosphere.integrator.p.precomputed.ᶜts, 1)
    ℂa = atmosphere.integrator.p.params.thermodynamics_params
    Ta = TD.air_temperature.(ℂa, tsa)
    pa = TD.air_pressure.(ℂa, tsa)
    qa = TD.total_specific_humidity.(ℂa, tsa)

    #=
    # TODO: make this work without allocation
    #       make sure that Remapper(args...; buffer_length=5)
    #       or whatever it needs to be
    exchange_fields = cat(ue, ve, Te, pe, qe, dims=3)
    atmos_fields    = [ua, va, Ta, pa, qa]
    CC.Remapping.interpolate!(exchange_fields, remapper, atmos_fields)
    =#

    CC.Remapping.interpolate!(ue, interpolator, ua)
    CC.Remapping.interpolate!(ve, interpolator, va)
    CC.Remapping.interpolate!(Te, interpolator, Ta)
    CC.Remapping.interpolate!(pe, interpolator, pa)
    CC.Remapping.interpolate!(qe, interpolator, qa)

    #=
    # This is needed, unless the above computations include the halos.
    # OC.fill_halo_regions!(exchange_atmosphere_state.u)
    # OC.fill_halo_regions!(exchange_atmosphere_state.v)
    # OC.fill_halo_regions!(exchange_atmosphere_state.T)
    # OC.fill_halo_regions!(exchange_atmosphere_state.q)
    # OC.fill_halo_regions!(exchange_atmosphere_state.p)
    =#

    return nothing
end

# Note: this just copies, for now.
KA.@kernel function _interpolate_atmosphere_state!(exchange_state, atmos_state)
    i, j = KA.@index(Global, NTuple)
    @inbounds begin
        exchange_state.u[i, j, 1] = atmos_state.u[i, j]
        exchange_state.v[i, j, 1] = atmos_state.v[i, j]
        exchange_state.T[i, j, 1] = atmos_state.T[i, j]
        exchange_state.q[i, j, 1] = atmos_state.q[i, j]
        exchange_state.p[i, j, 1] = atmos_state.p[i, j]
    end
end

#=
mutable struct AtmosphereExchanger
    atmosphere_to_exchange
    exchange_to_atmos
end
=#

function atmosphere_exchanger(atmosphere::ClimaAtmosSimulation, exchange_grid)
    λ = OC.λnodes(exchange_grid, OC.Center(), OC.Center(), OC.Center(), with_halos=true)
    φ = OC.φnodes(exchange_grid, OC.Center(), OC.Center(), OC.Center(), with_halos=true)

    if exchange_grid isa OC.LatitudeLongitudeGrid
        λ = reshape(λ, length(λ), 1)
        φ = reshape(φ, 1, length(φ))
    end

    Xh = @. CC.Geometry.LatLongPoint(φ, λ)
    space = axes(atmosphere.integrator.u.c)
    first_level = CC.Spaces.level(space, 1)

    # Note: buffer_length gives the maximum number of variables that can be remapped
    # within a single kernel.
    to_exchange_interp = CC.Remapping.Remapper(first_level, Xh, nothing, buffer_length=1)

    # Make a remapper for exchange_to_atmos regridding
    space3 = axes(atmosphere.integrator.p.precomputed.sfc_conditions.ts)
    space2 = CC.Spaces.SpectralElementSpace2D(space3.grid.full_grid.horizontal_grid)
    regridder = ClimaUtilities.Regridders.InterpolationsRegridder(space2)
    atmos_surface_points = regridder.coordinates

    if exchange_grid isa OC.Grids.OrthogonalSphericalShellGrid
        # One quick and dirty option: https://github.com/CliMA/OrthogonalSphericalShellGrids.jl/pull/29
        error("Not supported yet!")
    end

    dummy_flux = OC.Field{OC.Center, OC.Center, Nothing}(exchange_grid)
    Qc_a = map_interpolate(atmos_surface_points, dummy_flux)
    Qv_a = map_interpolate(atmos_surface_points, dummy_flux)
    Fv_a = map_interpolate(atmos_surface_points, dummy_flux)
    ρτx_a = map_interpolate(atmos_surface_points, dummy_flux)
    ρτy_a = map_interpolate(atmos_surface_points, dummy_flux)
    turbulent_atmosphere_surface_fluxes = (; Qc_a, Qv_a, Fv_a, ρτx_a, ρτy_a)

    return (; to_exchange_interp, turbulent_atmosphere_surface_fluxes, atmos_surface_points)
end

initialize!(::StateExchanger, ::ClimaAtmosSimulation) = nothing

#=
struct AtmosOceanExchanger{A2E, E2A}
    atmos_to_exchange_regridder :: A2E
    exchange_to_atmos_regridder :: E2A
end

function atmosphere_exchanger(atmosphere::ClimaAtmosSimulation, exchange_grid)
    space3 = axes(atmosphere.integrator.p.precomputed.sfc_conditions.ts)
    space2 = CC.Spaces.SpectralElementSpace2D(space3.grid.full_grid.horizontal_grid)
    regridder = ClimaUtilities.Regridders.InterpolationsRegridder(space2)

end
=#

using ClimaUtilities
using ClimaCore.Utilities: half

function map_interpolate!(cc_field, points, oc_field::OC.Field)
    loc = map(instantiate, OC.Fields.location(oc_field))
    grid = oc_field.grid
    data = oc_field.data

    map!(cc_field, points) do pt
        FT = eltype(pt)
        fᵢ = OC.Fields.interpolate(to_node(pt), data, loc, grid)
        convert(FT, fᵢ)
    end

    return nothing
end

# function compute_net_atmosphere_fluxes!(coupled_model::ClimaCoupledModel)
#     atmosphere = coupled_model.atmosphere
#     ocean = coupled_model.ocean
#     ocean_grid = ocean.model.grid
#     interfaces = coupled_model.interfaces
#     exchanger = interfaces.exchanger

#     atmos_surface_points = exchanger.atmosphere_exchanger.atmos_surface_points
#     (; Qc_a, Qv_a, Fv_a, ρτx_a, ρτy_a) = exchanger.atmosphere_exchanger.turbulent_atmosphere_surface_fluxes

#     Qv_e = interfaces.atmosphere_ocean_interface.fluxes.latent_heat
#     Qc_e = interfaces.atmosphere_ocean_interface.fluxes.sensible_heat
#     Fv_e = interfaces.atmosphere_ocean_interface.fluxes.water_vapor
#     ρτx_e = interfaces.atmosphere_ocean_interface.fluxes.x_momentum
#     ρτy_e = interfaces.atmosphere_ocean_interface.fluxes.y_momentum

#     map_interpolate!(Qc_a,  atmos_surface_points, Qc_e)
#     map_interpolate!(Qv_a,  atmos_surface_points, Qv_e)
#     map_interpolate!(Fv_a,  atmos_surface_points, Fv_e)
#     map_interpolate!(ρτx_a, atmos_surface_points, ρτx_e)
#     map_interpolate!(ρτy_a, atmos_surface_points, ρτy_e)

#     # Project onto a vector...
#     # :eyes https://github.com/CliMA/ClimaEarth.jl/pull/5/files
#     c = atmosphere.integrator.p.scratch.ᶠtemp_scalar
#     𝒢 = CC.Fields.level(CC.Fields.local_geometry_field(c), half)
#     ρwh = atmosphere.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot
#     @. ρwh = CA.SurfaceConditions.vector_from_component(Qv_a, 𝒢) +
#              CA.SurfaceConditions.vector_from_component(Qc_a, 𝒢)

#     # Mass or volume flux: check units
#     ρwq = atmosphere.integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot
#     @. ρwq = CA.SurfaceConditions.vector_from_component(Fv_a, 𝒢)

#     # TODO: validate this?
#     ρτ = atmosphere.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ
#     @. ρτ = tensor_from_uv_components(ρτx_a, ρτy_a, 𝒢)

#     return nothing
# end
