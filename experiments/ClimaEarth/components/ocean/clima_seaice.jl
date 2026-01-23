import Oceananigans as OC
import ClimaSeaIce as CSI
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions:
    IceWaterThermalEquilibrium, MeltingConstrainedFluxBalance, get_tracer, RadiativeEmission
import ClimaOcean as CO
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaComms
import ClimaCore as CC
import Thermodynamics as TD
import ClimaOcean.EN4: download_dataset
using KernelAbstractions: @kernel, @index, @inbounds

include("climaocean_helpers.jl")

# Rename ECCO password env variable to match ClimaOcean.jl
haskey(ENV, "ECCO_PASSWORD") && (ENV["ECCO_WEBDAV_PASSWORD"] = ENV["ECCO_PASSWORD"])

"""
    ClimaSeaIceSimulation{SIM, A, REMAP, NT, IP}

The ClimaCoupler simulation object used to run with ClimaSeaIce.
This type is used by the coupler to indicate that this simulation
is a surface/sea ice simulation for dispatch.

It contains the following objects:
- `ice::SIM`: The ClimaSeaIce simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
- `remapping::REMAP`: Objects needed to remap from the exchange (spectral) grid to Oceananigans spaces.
- `ocean_ice_fluxes::NT`: A NamedTuple of fluxes between the ocean and sea ice, computed at each coupling step.
- `ice_properties::IP`: A NamedTuple of sea ice properties, including melting speed, Stefan-Boltzmann constant,
    and the Celsius to Kelvin conversion constant.
"""
struct ClimaSeaIceSimulation{SIM, A, REMAP, NT, IP} <: Interfacer.SeaIceModelSimulation
    ice::SIM
    area_fraction::A
    remapping::REMAP
    ocean_ice_fluxes::NT
    ice_properties::IP
end

"""
    ConcentrationMaskedRadiativeEmission

A heat boundary condition that emits radiation only where the sea ice
concentration is greater than zero. This is needed to prevent radiative
emission where we have no sea ice.
"""
struct ConcentrationMaskedRadiativeEmission{FT}
    emissivity::FT
    stefan_boltzmann_constant::FT
    reference_temperature::FT
end

function ConcentrationMaskedRadiativeEmission(
    FT = Float64;
    emissivity = 1,
    stefan_boltzmann_constant = 5.67e-8,
    reference_temperature = 273.15,
)

    return ConcentrationMaskedRadiativeEmission(
        convert(FT, emissivity),
        convert(FT, stefan_boltzmann_constant),
        convert(FT, reference_temperature),
    )
end

function CSI.SeaIceThermodynamics.HeatBoundaryConditions.getflux(
    emission::ConcentrationMaskedRadiativeEmission,
    i,
    j,
    grid,
    T,
    clock,
    fields,
)
    ϵ = emission.emissivity
    σ = emission.stefan_boltzmann_constant
    Tᵣ = emission.reference_temperature
    @inbounds ℵij = fields.ℵ[i, j, 1]
    return ϵ * σ * (T + Tᵣ)^4 * (ℵij > 0)
end


"""
    ClimaSeaIceSimulation()

Creates an ClimaSeaIceSimulation object containing a model, an integrator, and
a surface area fraction field.

If a start date is provided, we initialize the sea ice concentration and thickness
using the ECCO4Monthly dataset. If no start date is provided, we initialize with zero sea ice.

Since this model does not solve for prognostic temperature, we use a
prescribed heat flux boundary condition at the top, which is used to solve for
top surface temperature and longwave emission together.
The surface temperature from the last step is provided to the coupler to be
used in computing fluxes.

Specific details about the default model configuration
can be found in the documentation for `ClimaSeaIce.sea_ice_simulation`.

# Arguments
- `ocean`: [OceananigansSimulation] the ocean simulation to couple with.
- `output_dir`: [String] the directory to save output files.
- `start_date`: [Date] the start date to initialize the sea ice concentration and thickness.
- `coupled_param_dict`: [Dict{String, Any}] the coupled parameters.
- `Δt`: [Float64] the time step.
"""
function ClimaSeaIceSimulation(
    ocean;
    output_dir,
    start_date = nothing,
    coupled_param_dict = CP.create_toml_dict(eltype(ocean.area_fraction)),
    Δt = 5 * 60.0, # 5 minutes
)
    # Initialize the sea ice with the same grid as the ocean
    grid = ocean.ocean.model.grid
    arch = OC.Architectures.architecture(grid)
    advection = ocean.ocean.model.advection.T

    ice = sea_ice_simulation(grid, ocean.ocean; Δt, advection)

    # Initialize nonzero sea ice if start date provided
    if !isnothing(start_date)
        sic_metadata = CO.DataWrangling.Metadatum(
            :sea_ice_concentration,
            dataset = CO.DataWrangling.ECCO.ECCO4Monthly(),
            date = start_date,
        )
        h_metadata = CO.DataWrangling.Metadatum(
            :sea_ice_thickness,
            dataset = CO.DataWrangling.ECCO.ECCO4Monthly(),
            date = start_date,
        )

        OC.set!(ice.model.ice_concentration, sic_metadata)
        OC.set!(ice.model.ice_thickness, h_metadata)
    end

    # Get sea ice properties from coupled parameters
    ice_properties = (;
        melting_speed = 1e-4,
        σ = coupled_param_dict["stefan_boltzmann_constant"],
        C_to_K = coupled_param_dict["temperature_water_freeze"],
    )

    # Since ocean and sea ice share the same grid, we can also share the remapping objects
    remapping = ocean.remapping

    # Before version 0.96.22, the NetCDFWriter was broken on GPU
    if arch isa OC.CPU
        # Save all tracers and velocities to a NetCDF file at daily frequency
        outputs = OC.prognostic_fields(ice.model)
        jld_writer = OC.JLD2Writer(
            ice.model,
            outputs;
            schedule = OC.TimeInterval(86400), # Daily output
            filename = joinpath(output_dir, "seaice_diagnostics.jld2"),
            overwrite_existing = true,
            array_type = Array{Float32},
        )
        ice.output_writers[:diagnostics] = jld_writer
    end

    # Allocate space for the sea ice-ocean (io) fluxes
    io_bottom_heat_flux = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    io_frazil_heat_flux = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    io_salt_flux = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    x_momentum = OC.Field{OC.Face, OC.Center, Nothing}(grid)
    y_momentum = OC.Field{OC.Center, OC.Face, Nothing}(grid)

    ocean_ice_fluxes = (
        interface_heat = io_bottom_heat_flux,
        frazil_heat = io_frazil_heat_flux,
        salt = io_salt_flux,
        x_momentum = x_momentum,
        y_momentum = y_momentum,
    )

    # Get the initial area fraction from the fractional ice concentration
    boundary_space = axes(ocean.area_fraction)
    area_fraction = Interfacer.remap(boundary_space, ice.model.ice_concentration)

    sim = ClimaSeaIceSimulation(
        ice,
        area_fraction,
        remapping,
        ocean_ice_fluxes,
        ice_properties,
    )

    # Ensure ocean temperature is above freezing where there is sea ice
    CO.OceanSeaIceModels.above_freezing_ocean_temperature!(ocean.ocean, ice)
    return sim
end

function sea_ice_simulation(
    grid,
    ocean = nothing;
    Δt = 5 * 60.0, # 5 minutes
    ice_salinity = 4, # psu
    advection = nothing, # for the moment
    tracers = (),
    ice_heat_capacity = 2100, # J kg⁻¹ K⁻¹
    ice_consolidation_thickness = 0.05, # m
    ice_density = 900, # kg m⁻³
    dynamics = CO.SeaIceSimulations.sea_ice_dynamics(grid, ocean),
    phase_transitions = CSI.PhaseTransitions(; ice_heat_capacity, ice_density),
    conductivity = 2, # kg m s⁻³ K⁻¹
    internal_heat_flux = CSI.ConductiveFlux(; conductivity),
)

    top_surface_temperature = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    top_heat_boundary_condition = MeltingConstrainedFluxBalance()
    kᴺ = size(grid, 3)
    surface_ocean_salinity = OC.interior(ocean.model.tracers.S, :, :, kᴺ:kᴺ)
    bottom_heat_boundary_condition = IceWaterThermalEquilibrium(surface_ocean_salinity)

    ice_thermodynamics = CSI.SlabSeaIceThermodynamics(
        grid;
        internal_heat_flux,
        phase_transitions,
        top_heat_boundary_condition,
        bottom_heat_boundary_condition,
    )

    bottom_heat_flux = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    top_heat_flux = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    top_heat_flux = (top_heat_flux, ConcentrationMaskedRadiativeEmission())

    # Build the sea ice model
    sea_ice_model = CSI.SeaIceModel(
        grid;
        ice_salinity,
        advection,
        tracers,
        ice_consolidation_thickness,
        ice_thermodynamics,
        dynamics,
        bottom_heat_flux,
        top_heat_flux,
    )

    # Build the simulation
    return OC.Simulation(sea_ice_model; Δt)
end

###############################################################################
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###############################################################################

# Timestep the simulation forward to time `t`
Interfacer.step!(sim::ClimaSeaIceSimulation, t) =
    OC.time_step!(sim.ice, float(t) - sim.ice.model.clock.time)

Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:area_fraction}) = sim.area_fraction
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:ice_concentration}) =
    sim.ice.model.ice_concentration

# At the moment, we return always Float32. This is because we always want to run
# Oceananingans with Float64, so we have no way to know the float type here. Sticking with
# Float32 ensures that nothing is accidentally promoted to Float64. We will need to change
# this anyway.
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:roughness_buoyancy}) =
    Float32(5.8e-5)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:roughness_momentum}) =
    Float32(5.8e-5)
# Specify constant roughness model for ClimaSeaIceSimulation (sea ice)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:roughness_model}) = :constant
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:emissivity}) = Float32(1)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_direct_albedo}) =
    Float32(0.7)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_diffuse_albedo}) =
    Float32(0.7)

# Approximate the sea ice surface temperature as the temperature computed from the
#  fluxes at the previous timestep.
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_temperature}) =
    sim.ice_properties.C_to_K + sim.ice.model.ice_thermodynamics.top_surface_temperature

"""
    FluxCalculator.update_turbulent_fluxes!(sim::ClimaSeaIceSimulation, fields)

Update the turbulent fluxes in the simulation using the values stored in the coupler fields.
These include latent heat flux, sensible heat flux, momentum fluxes, and moisture flux.

The input `fields` are already area-weighted, so there's no need to weight them again.

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
    grid = sim.ice.model.grid
    ice_concentration = sim.ice.model.ice_concentration

    # Convert the momentum fluxes from contravariant to Cartesian basis
    contravariant_to_cartesian!(sim.remapping.temp_uv_vec, F_turb_ρτxz, F_turb_ρτyz)
    F_turb_ρτxz_uv = sim.remapping.temp_uv_vec.components.data.:1
    F_turb_ρτyz_uv = sim.remapping.temp_uv_vec.components.data.:2

    # Remap momentum fluxes onto reduced 2D Center, Center fields using scratch arrays and fields
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        F_turb_ρτxz_uv,
    )
    OC.set!(sim.remapping.scratch_cc1, sim.remapping.scratch_arr1) # zonal momentum flux
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr2,
        sim.remapping.remapper_cc,
        F_turb_ρτyz_uv,
    )
    OC.set!(sim.remapping.scratch_cc2, sim.remapping.scratch_arr2) # meridional momentum flux

    # Rename for clarity; these are now Center, Center Oceananigans fields
    F_turb_ρτxz_cc = sim.remapping.scratch_cc1
    F_turb_ρτyz_cc = sim.remapping.scratch_cc2

    # Set the momentum flux BCs at the correct locations using the remapped scratch fields
    # Note that this requires the sea ice model to always be run with dynamics turned on
    si_flux_u = sim.ice.model.dynamics.external_momentum_stresses.top.u
    si_flux_v = sim.ice.model.dynamics.external_momentum_stresses.top.v
    set_from_extrinsic_vector!(
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
    si_flux_heat = sim.ice.model.external_heat_fluxes.top[1]
    OC.interior(si_flux_heat, :, :, 1) .+=
        (OC.interior(ice_concentration, :, :, 1) .> 0) .* (remapped_F_lh .+ remapped_F_sh)

    return nothing
end

function Interfacer.update_field!(sim::ClimaSeaIceSimulation, ::Val{:area_fraction}, field)
    sim.area_fraction .= field
    return nothing
end

"""
    FieldExchanger.update_sim!(sim::ClimaSeaIceSimulation, csf)

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
function FieldExchanger.update_sim!(sim::ClimaSeaIceSimulation, csf)
    ice_concentration = sim.ice.model.ice_concentration

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
    si_flux_heat = sim.ice.model.external_heat_fluxes.top[1]
    α = Interfacer.get_field(sim, Val(:surface_direct_albedo)) # scalar
    ϵ = Interfacer.get_field(sim, Val(:emissivity)) # scalar

    # Update only where ice concentration is greater than zero.
    OC.interior(si_flux_heat, :, :, 1) .=
        (OC.interior(ice_concentration, :, :, 1) .> 0) .*
        (-(1 .- α) .* remapped_SW_d .- ϵ .* remapped_LW_d)
    return nothing
end

"""
    ocean_seaice_fluxes!(ocean_sim::OceananigansSimulation, ice_sim::ClimaSeaIceSimulation)

Compute the fluxes between the ocean and sea ice, storing them in the `ocean_ice_fluxes`
fields of the ocean and sea ice simulations.

This function assumes both simulations share the same grid, so no remapping is done.

Both simulations have had their atmospheric fluxes updated already in this timestep
(see `update_sim!` and `update_turbulent_fluxes!`), so we add the contributions from the
ocean-sea ice interactions to the existing fluxes, rather than overwriting all fluxes.

!!! note
    This function must be called after the turbulent fluxes have been updated in both
    simulations. Here only the contributions from the sea ice/ocean interactions
    are added to the fluxes.
"""
function FluxCalculator.ocean_seaice_fluxes!(
    ocean_sim::OceananigansSimulation,
    ice_sim::ClimaSeaIceSimulation,
)
    melting_speed = ice_sim.ice_properties.melting_speed
    ocean_properties = ocean_sim.ocean_properties
    ice_concentration = Interfacer.get_field(ice_sim, Val(:ice_concentration))

    # Update the sea ice concentration in the ocean simulation
    ocean_sim.ice_concentration .= ice_concentration

    # Compute the fluxes and store them in the both simulations
    CO.OceanSeaIceModels.InterfaceComputations.compute_sea_ice_ocean_fluxes!(
        ice_sim.ocean_ice_fluxes,
        ocean_sim.ocean,
        ice_sim.ice,
        melting_speed,
        ocean_properties,
    )

    ## Update the internals of the sea ice model
    # Set the bottom heat flux to the sum of the frazil and interface heat fluxes
    bottom_heat_flux = ice_sim.ice.model.external_heat_fluxes.bottom

    Qf = ice_sim.ocean_ice_fluxes.frazil_heat        # frazil heat flux
    Qi = ice_sim.ocean_ice_fluxes.interface_heat     # interfacial heat flux
    bottom_heat_flux .= Qf .+ Qi

    ## Update the internals of the ocean model
    ρₒ⁻¹ = 1 / ocean_sim.ocean_properties.reference_density
    cₒ = ocean_sim.ocean_properties.heat_capacity

    # Compute fluxes for u, v, T, and S from momentum, heat, and freshwater fluxes
    oc_flux_u = surface_flux(ocean_sim.ocean.model.velocities.u)
    oc_flux_v = surface_flux(ocean_sim.ocean.model.velocities.v)

    ρτxio = ice_sim.ocean_ice_fluxes.x_momentum # sea_ice - ocean zonal momentum flux
    ρτyio = ice_sim.ocean_ice_fluxes.y_momentum # sea_ice - ocean meridional momentum flux

    # Update the momentum flux contributions from ocean/sea ice fluxes
    grid = ocean_sim.ocean.model.grid
    arch = OC.Architectures.architecture(grid)
    OC.Utils.launch!(
        arch,
        grid,
        :xy,
        _add_ocean_ice_stress!,
        oc_flux_u,
        oc_flux_v,
        grid,
        ρτxio,
        ρτyio,
        ρₒ⁻¹,
        ice_concentration,
    )

    oc_flux_T = surface_flux(ocean_sim.ocean.model.tracers.T)
    OC.interior(oc_flux_T, :, :, 1) .+=
        OC.interior(ice_concentration, :, :, 1) .* OC.interior(Qi, :, :, 1) .* ρₒ⁻¹ ./ cₒ

    oc_flux_S = surface_flux(ocean_sim.ocean.model.tracers.S)
    OC.interior(oc_flux_S, :, :, 1) .+=
        OC.interior(ice_concentration, :, :, 1) .*
        OC.interior(ice_sim.ocean_ice_fluxes.salt, :, :, 1)

    return nothing
end

"""
    _add_ocean_ice_stress!(oc_flux_u, oc_flux_v, grid, ρτxio, ρτyio, ρₒ⁻¹, ice_concentration)

Add the contribution from the ocean-ice stress to the surface fluxes for each
component of the ocean velocity (u and v).

Arguments:
- `oc_flux_u`: [Field] the surface flux for the ocean zonal velocity.
- `oc_flux_v`: [Field] the surface flux for the ocean meridional velocity.
- `grid`: [Grid] the grid used by ocean and ice.
- `ρτxio`: [Field] the ice-ocean zonal momentum flux.
- `ρτyio`: [Field] the ice-ocean meridional momentum flux.
- `ρₒ⁻¹`: [Float] the inverse of the ocean reference density.
- `ice_concentration`: [Field] the sea ice concentration.
"""
@kernel function _add_ocean_ice_stress!(
    oc_flux_u,
    oc_flux_v,
    grid,
    ρτxio,
    ρτyio,
    ρₒ⁻¹,
    ice_concentration,
)
    i, j = @index(Global, NTuple)

    # ℑxᶠᵃᵃ: interpolate faces to centers
    oc_flux_u[i, j, 1] +=
        ρτxio[i, j, 1] * ρₒ⁻¹ * OC.Operators.ℑxᶠᵃᵃ(i, j, 1, grid, ice_concentration)
    oc_flux_v[i, j, 1] +=
        ρτyio[i, j, 1] * ρₒ⁻¹ * OC.Operators.ℑyᵃᶠᵃ(i, j, 1, grid, ice_concentration)
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
