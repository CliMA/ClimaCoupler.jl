import Oceananigans
import ClimaOcean as Ocean
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaComms

"""
    OceananigansSimulation{M, I, A}

The ClimaCoupler simulation object used to run with Oceananigans.
This type is used by the coupler to indicate that this simulation
is an surface/ocean simulation for dispatch.

It contains the following objects:
- `sim::M`: The Oceananigans simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
"""
struct OceananigansSimulation{M, I, A} <: Interfacer.OceanModelSimulation
    sim::M
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
    ::Type{FT}, # TODO decide which arguments we want here
    dt::TT,
    tspan::Tuple{TT, TT},
    start_date::Dates.DateTime,
    output_dir::String,
    area_fraction,
    comms_ctx = ClimaComms.context(),
) where {FT, TT <: Union{Float64, ITime}}
    # TODO fill this out
    arch = if comms_ctx.device isa ClimaComms.CUDADevice
        Oceananigans.GPU()
    else
        Oceananigans.CPU()
    end

    # Set up ocean grid
    Nx = 360
    Ny = 170
    Nz = 30
    z_faces = Ocean.exponential_z_faces(; Nz, h = 30, depth = 6000)

    # TODO how to specify FT?
    grid = Oceananigans.LatitudeLongitudeGrid(
        arch;
        size = (Nx, Ny, Nz),
        longitude = (0, 360),
        latitude = (-85, 85),
        z = z_faces,
        halo = (7, 7, 7),
    )

    # Choose parameterizations
    momentum_advection = Oceananigans.WENOVectorInvariant(order = 5)
    tracer_advection = Oceananigans.WENO(order = 5)
    free_surface = Oceananigans.SplitExplicitFreeSurface(grid; substeps = 30)

    # Create ocean simulation
    ocean = Ocean.ocean_simulation(grid; momentum_advection, tracer_advection, free_surface, warn = false)


    # Set up initial conditions for temperature and salinity
    Tᵢ(λ, φ, z) = 30 * (1 - tanh((abs(φ) - 30) / 5)) / 2 + rand()
    Sᵢ(λ, φ, z) = 30 - 5e-3 * z + rand()
    Oceananigans.set!(ocean.model, T = Tᵢ, S = Sᵢ)

    # TODO do we want to store sim or model?
    return OceananigansSimulation(ocean, area_fraction)
end

Interfacer.name(::OceananigansSimulation) = "OceananigansSimulation"

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

# Timestep the simulation forward to time `t`
Interfacer.step!(sim::OceananigansSimulation, t) = Oceananigans.time_step!(sim.model, t - sim.model.clock.time)

# Reset prognostic state and current time to initial conditions
Interfacer.reinit!(sim::OceananigansSimulation, t) = nothing # TODO fill this out

"""
    Interfacer.get_field(sim::OceananigansSimulation, ::Val{:_})

Get the value of the specified field in the ocean simulation.
This is used to exchange properties from the ocean to the coupler and other components.
"""
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:area_fraction}) = sim.area_fraction
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_buoyancy}) = return nothing # TODO fill this out
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:roughness_momentum}) = return nothing # TODO fill this out
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_direct_albedo}) = return nothing # TODO fill this out
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_diffuse_albedo}) = return nothing # TODO fill this out
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_humidity}) = return nothing # TODO fill this out
Interfacer.get_field(sim::OceananigansSimulation, ::Val{:surface_temperature}) = return nothing # TODO fill this out

# These two methods are used to track conservation of the coupled system, so it's okay to leave them empty for now
function Interfacer.get_field(sim::OceananigansSimulation, ::Val{:energy}, level)
    # TODO fill this out
    return nothing
end
function Interfacer.get_field(sim::OceananigansSimulation, ::Val{:water}, level)
    # TODO fill this out
    return nothing
end

"""
    Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:_}, field)

Update the specified field in the ocean simulation with the provided field, which
is a ClimaCore Field defined on the boundary space. Remapping will be needed here.
This is used to exchange properties from the coupler and other components to the ocean.
"""
# Setters to update fields stored in the ocean simulation with the provided field (a ClimaCore Field defined on boundary_space)
# These are just the `update_field!` methods corresponding to the default coupler fields,
#  but if the ocean requires more fields, we can add them here. (maybe air temperature, pressure, humidity, CO2?)
# If any of these aren't needed by the ocean model, we can remove them.
function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:air_density}, field)
    # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
end
function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:area_fraction}, field)
    parent(sim.area_fraction) .= parent(field)
end
function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:liquid_precipitation}, field)
    # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
end
function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:snow_precipitation}, field)
    # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
end
function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:lw_d}, field)
    # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
end
function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:sw_d}, field)
    # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
end
function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:turbulent_energy_flux}, field)
    # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
    #  or we can update the model directly in the flux calculation, depending on where we compute fluxes
end
function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:turbulent_moisture_flux}, field)
    # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
    #  or we can update the model directly in the flux calculation, depending on where we compute fluxes
end
function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:turbulent_momentum_flux_x}, field)
    # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
    #  or we can update the model directly in the flux calculation, depending on where we compute fluxes
end
function Interfacer.update_field!(sim::OceananigansSimulation, ::Val{:turbulent_momentum_flux_y}, field)
    # TODO fill this out - remap `field` to the ocean grid and update the correct location in the ocean model
    #  or we can update the model directly in the flux calculation, depending on where we compute fluxes
end

"""
    FieldExchanger.update_sim!(sim::OceananigansSimulation, csf, turbulent_fluxes, area_fraction)

Update the ocean simulation with the provided fields, which have been filled in
by the coupler.
"""
function FieldExchanger.update_sim!(sim::OceananigansSimulation, csf, turbulent_fluxes, area_fraction)
    Interfacer.update_field!(sim, Val(:air_density), air_density)
    Interfacer.update_field!(sim, Val(:area_fraction), area_fraction)

    # precipitation
    Interfacer.update_field!(sim, Val(:liquid_precipitation), csf.P_liq)
    Interfacer.update_field!(sim, Val(:snow_precipitation), csf.P_snow)

    # update fields for radiative transfer
    Interfacer.update_field!(sim, Val(:sw_d), csf.SW_d)
    Interfacer.update_field!(sim, Val(:lw_d), csf.LW_d)

    # TODO update other fields as needed
end

"""
    update_turbulent_fluxes!(sim::OceananigansSimulation, fields::NamedTuple)

Update the turbulent fluxes in the ocean simulation with the provided fields,
which were filled in during the coupler's flux calculation.

TODO we don't need this function if we update the model directly in the flux calculation
"""
function FluxCalculator.update_turbulent_fluxes!(sim::OceananigansSimulation, fields::NamedTuple)
    Interfacer.update_field!(sim, Val(:F_turb_ρτxz), fields.F_turb_ρτxz)
    Interfacer.update_field!(sim, Val(:F_turb_ρτyz), fields.F_turb_ρτyz)
    Interfacer.update_field!(sim, Val(:F_turb_energy), fields.F_turb_energy)
    Interfacer.update_field!(sim, Val(:F_turb_moisture), fields.F_turb_moisture)
    return nothing
end

"""
    Interfacer.add_coupler_fields!(coupler_field_names, ::OceananigansSimulation)

Extend Interfacer.add_coupler_fields! to add the fields required for
OceananigansSimulation, if any.
"""
function Interfacer.add_coupler_fields!(coupler_field_names, ::OceananigansSimulation)
    # TODO check Interfacer docs for default exchange fields, fill this out if needed
    ocean_coupler_fields = [:SW_d, :LW_d]
    push!(coupler_field_names, ocean_coupler_fields...)
end

"""
    FieldExchanger.import_atmos_fields!(csf, sim::OceananigansSimulation, atmos_sim, turbulent_fluxes)

Update the coupler fields in-place with the values from the atmosphere simulation.
This is defined here because the coupler exchange fields are specified by the
Oceananigans simulation `add_coupler_fields!` method, and we need to know which fields
are required from the atmosphere model.
"""
function FieldExchanger.import_atmos_fields!(csf, sim::OceananigansSimulation, atmos_sim, turbulent_fluxes)
    # TODO fix remap function calls

    # radiative fluxes
    remap!(csf.SW_d, Interfacer.get_field(atmos_sim, Val(:sw_d)))
    remap!(csf.LW_d, Interfacer.get_field(atmos_sim, Val(:lw_d)))
    # precipitation
    remap!(csf.P_liq, Interfacer.get_field(atmos_sim, Val(:liquid_precipitation)))
    remap!(csf.P_snow, Interfacer.get_field(atmos_sim, Val(:snow_precipitation)))
    # air density
    remap!(csf.air_density, Interfacer.get_field(atmos_sim, Val(:air_density)))

    # TODO import other fields as needed (according to how we extend `add_coupler_fields!`)
end

## Extend functions for ocean-specific flux calculation
"""
    compute_surface_fluxes!(csf, sim::OceananigansSimulation, atmos_sim, boundary_space, thermo_params, surface_scheme)

This function computes surface fluxes between the Oceananigans simulation and
the atmosphere.

Update the input coupler surface fields `csf` in-place with the computed fluxes
for this model. These are then summed using area-weighting across all surface
models to get the total fluxes.

# Arguments
- `csf`: [CC.Fields.Field] containing a NamedTuple of turbulent flux fields: `F_turb_ρτxz`, `F_turb_ρτyz`, `F_turb_energy`, `F_turb_moisture`.
- `sim`: [OceananigansSimulation] the ocean simulation to compute fluxes for.
- `atmos_sim`: [Interfacer.AtmosModelSimulation] the atmosphere simulation to compute fluxes with.
- unused arguments: `boundary_space`, `thermo_params`, `surface_scheme`
"""
function FluxCalculator.compute_surface_fluxes!(
    csf,
    sim::OceananigansSimulation,
    atmos_sim::Interfacer.AtmosModelSimulation,
    _...,
)
    # TODO fill this in -
    #  compute F_turb_ρτxz, F_turb_ρτyz, F_shf, F_lhf, F_turb_moisture using Oceananigans.jl functions
    #  something like this:

    fluxes = Oceananigans.compute_surface_fluxes(p, sim.model, Y, t, atmos_sim.integrator)
    (; F_turb_ρτxz, F_turb_ρτyz, F_shf, F_lhf, F_turb_moisture) = fluxes

    # get area fraction (min = 0, max = 1)
    area_fraction = Interfacer.get_field(sim, Val(:area_fraction))

    # add the flux contributing from this surface to the coupler field
    # note that the fluxes are area-weighted, so if a surface model is
    #  not present at this point, the fluxes are zero
    @. csf.F_turb_ρτxz += F_turb_ρτxz * area_fraction
    @. csf.F_turb_ρτyz += F_turb_ρτyz * area_fraction
    @. csf.F_turb_energy += (F_shf .+ F_lhf) * area_fraction
    @. csf.F_turb_moisture += F_turb_moisture * area_fraction
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
    error("get_model_prog_state not implemented")
end
