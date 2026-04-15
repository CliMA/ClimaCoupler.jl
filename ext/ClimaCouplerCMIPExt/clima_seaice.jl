using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions:
    IceWaterThermalEquilibrium, PrescribedTemperature, get_tracer
import ClimaComms
import ClimaOcean.EN4: download_dataset
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Dates
import ClimaUtilities.TimeManager: ITime, date, counter, period
using StaticArrays

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
- `ocean_ice_interface::NT`: A NamedTuple containing fluxes between the ocean and sea ice, computed at each coupling step,
                             the interfacial temperature and salinity, and the flux formulation used to compute the fluxes.
- `ice_properties::IP`: A NamedTuple of sea ice properties, including melting speed, Stefan-Boltzmann constant,
    and the Celsius to Kelvin conversion constant.
"""
struct ClimaSeaIceSimulation{SIM, A, REMAP, NT, IP, MDT} <:
       Interfacer.AbstractSeaIceSimulation
    ice::SIM
    area_fraction::A
    remapping::REMAP
    ocean_ice_interface::NT
    ice_properties::IP
    model_Δt::MDT
end

"""
    Interfacer.SeaIceSimulation(::Type{FT}, ::Val{:clima_seaice}; kwargs...)

Extension of the generic SeaIceSimulation constructor for ClimaSeaIce.
"""
function Interfacer.SeaIceSimulation(::Type{FT}, ::Val{:clima_seaice}; kwargs...) where {FT}
    return ClimaSeaIceSimulation(FT; kwargs...)
end

"""
    ClimaSeaIceSimulation()

Creates a ClimaSeaIceSimulation object containing a model, an integrator, and
a surface area fraction field.

If a start date is provided, we initialize the sea ice concentration and thickness
using the ECCO4Monthly dataset. If no start date is provided, we initialize with zero sea ice.

The top heat boundary condition is `PrescribedTemperature`: the coupler iteratively
diagnoses T_sfc via the `update_T_sfc` callback and writes it back to the ice model
after each flux computation (see `compute_surface_fluxes!`).

# Arguments
- `ocean`: [OceananigansSimulation] the ocean simulation to couple with.
- `output_dir`: [String] the directory to save output files.
- `start_date`: [Date] the start date to initialize the sea ice concentration and thickness.
- `coupled_param_dict`: [Dict{String, Any}] the coupled parameters.
- `dt`: [Float64] the time step.
"""
function ClimaSeaIceSimulation(
    ::Type{FT};
    ocean,
    output_dir,
    start_date = nothing,
    coupled_param_dict = CP.create_toml_dict(FT),
    dt = 5 * 60.0, # 5 minutes
    extra_kwargs...,
) where {FT}
    # Initialize the sea ice with the same grid as the ocean
    grid = ocean.ocean.model.grid
    arch = OC.Architectures.architecture(grid)

    advection = ocean.ocean.model.advection.T
    ice = CO.SeaIces.sea_ice_simulation(
        grid,
        ocean.ocean;
        clock = deepcopy(ocean.ocean.model.clock),
        Δt = float(dt),
        advection,
    )

    ocean_ice_flux_formulation =
        CO.OceanSeaIceModels.InterfaceComputations.ThreeEquationHeatFlux(ice)
    interface_temperature = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    interface_salinity = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    # Initialize model_Δt so that time stepping works properly
    model_Δt = dt

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
            array_type = Array{FT},
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

    # `ClimaOcean.compute_sea_ice_ocean_fluxes` expects a NamedTuple containing fluxes,
    # flux_formulation, temperature, and salinity.
    ocean_ice_interface = (;
        fluxes = ocean_ice_fluxes,
        flux_formulation = ocean_ice_flux_formulation,
        temperature = interface_temperature,
        salinity = interface_salinity,
    )

    # Build the ocean - sea ice interface object

    # Get the initial area fraction from the fractional ice concentration
    boundary_space = axes(ocean.area_fraction)
    area_fraction = Interfacer.remap(boundary_space, ice.model.ice_concentration)

    sim = ClimaSeaIceSimulation(
        ice,
        area_fraction,
        remapping,
        ocean_ice_interface,
        ice_properties,
        model_Δt,
    )

    # Ensure ocean temperature is above freezing where there is sea ice
    CO.OceanSeaIceModels.above_freezing_ocean_temperature!(ocean.ocean, grid, ice)
    return sim
end

###############################################################################
### Functions required by ClimaCoupler.jl for a AbstractSurfaceSimulation
###############################################################################

# Timestep the simulation forward to time `t`
function Interfacer.step!(sim::ClimaSeaIceSimulation, t::Float64)
    Δt = t - sim.ice.model.clock.time
    if isapprox(Δt, sim.model_Δt, atol=0.125) || Δt > sim.model_Δt
        OC.time_step!(sim.ice, Δt)
    end
    return nothing
end

function Interfacer.step!(sim::ClimaSeaIceSimulation, t::ITime)
    Δt_msec = date(t) - sim.ice.model.clock.time
    model_Δt_msec = counter(sim.model_Δt) * Dates.Millisecond(period(sim.model_Δt))
    if Δt_msec >= model_Δt_msec
        OC.time_step!(sim.ice, float(sim.model_Δt))
    end
    return nothing
end

Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:area_fraction}) = sim.area_fraction
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:ice_concentration}) =
    sim.ice.model.ice_concentration
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:ice_thickness}) =
    sim.ice.model.ice_thickness
# Internal temperature in Kelvin (same convention as :surface_temperature)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:internal_temperature}) =
    sim.ocean_ice_interface.temperature + sim.ice_properties.C_to_K

# TODO better values for roughness
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:roughness_model}) = :constant
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:roughness_buoyancy}) =
    eltype(sim.ice.model)(5.8e-5)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:roughness_momentum}) =
    eltype(sim.ice.model)(5.8e-5)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:emissivity}) =
    eltype(sim.ice.model)(1)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_direct_albedo}) =
    eltype(sim.ice.model)(0.7)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_diffuse_albedo}) =
    eltype(sim.ice.model)(0.7)

# Surface temperature from last timestep (Celsius → Kelvin)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_temperature}) =
    sim.ice_properties.C_to_K + sim.ice.model.ice_thermodynamics.top_surface_temperature

"""
    FluxCalculator.compute_surface_fluxes!(csf, sim::ClimaSeaIceSimulation, atmos_sim, thermo_params)

Compute surface fluxes for `ClimaSeaIceSimulation`, iteratively diagnosing T_sfc
via the `update_T_sfc` callback to satisfy the skin-temperature flux balance.

The diagnosed T_sfc is written back to ClimaSeaIce's `top_surface_temperature`
(used by `PrescribedTemperature`) so the ice thermodynamics stays consistent.
"""
function FluxCalculator.compute_surface_fluxes!(
    csf,
    sim::ClimaSeaIceSimulation,
    atmos_sim::Interfacer.AbstractAtmosSimulation,
    thermo_params,
)
    boundary_space = axes(csf)
    FT = CC.Spaces.undertype(boundary_space)
    surface_fluxes_params = FluxCalculator.get_surface_params(atmos_sim)

    uv_int = StaticArrays.SVector.(csf.u_int, csf.v_int)

    # Sea ice parameters for the update_T_sfc callback (load into boundary-space scratch
    # Fields first to avoid GPU scalar indexing when building the callback)
    Interfacer.get_field!(csf.scalar_temp1, sim, Val(:ice_thickness))
    Interfacer.get_field!(csf.scalar_temp2, sim, Val(:internal_temperature))
    Interfacer.get_field!(csf.scalar_temp3, sim, Val(:emissivity))
    Interfacer.get_field!(csf.scalar_temp4, sim, Val(:surface_direct_albedo))
    δ = csf.scalar_temp1
    T_i = csf.scalar_temp2
    ϵ = csf.scalar_temp3
    α_albedo = csf.scalar_temp4
    internal_heat_flux = sim.ice.model.ice_thermodynamics.internal_heat_flux
    κ = if hasfield(typeof(internal_heat_flux), :conductivity)
        FT.(internal_heat_flux.conductivity)
    else
        convert(FT, 2) # default conductivity [W m⁻¹ K⁻¹]
    end
    σ = FT(sim.ice_properties.σ)
    SW_d = csf.SW_d
    LW_d = csf.LW_d
    T_melt = FT(sim.ice_properties.C_to_K) # Melting temperature (freezing point of water)

    # Build element-wise update_T_sfc callbacks (each closes over local ice parameters)
    update_T_sfc_callback =
        ClimaCouplerCMIPExt.update_T_sfc.(κ, δ, T_i, σ, ϵ, SW_d, LW_d, α_albedo, T_melt)

    # Surface temperature guess from last timestep
    Interfacer.get_field!(csf.scalar_temp1, sim, Val(:surface_temperature))
    T_sfc = csf.scalar_temp1

    # Surface humidity
    ρ_sfc =
        SF.surface_density.(
            surface_fluxes_params,
            csf.T_atmos,
            csf.ρ_atmos,
            T_sfc,
            csf.height_int .- csf.height_sfc,
            csf.q_tot_atmos,
            0,
            0,
        )
    csf.scalar_temp2 .= TD.q_vap_saturation.(thermo_params, T_sfc, ρ_sfc, 0, 0)
    q_sfc = csf.scalar_temp2

    # Roughness and gustiness configuration
    gustiness = ones(boundary_space)
    roughness_params = FluxCalculator.get_roughness_params(csf, sim)
    config = SF.SurfaceFluxConfig.(roughness_params, SF.ConstantGustinessSpec.(gustiness))

    fluxes =
        FluxCalculator.get_surface_fluxes.(
            surface_fluxes_params,
            uv_int,
            csf.T_atmos,
            csf.q_tot_atmos,
            csf.q_liq_atmos,
            csf.q_ice_atmos,
            csf.ρ_atmos,
            csf.height_int,
            uv_int .* FT(0),
            T_sfc,
            q_sfc,
            csf.height_sfc,
            FT(0),
            config,
            update_T_sfc_callback,
        )

    FluxCalculator.update_flux_fields!(csf, sim, fluxes)
    area_fraction = Interfacer.get_field(sim, Val(:area_fraction))

    # Write diagnosed T_sfc back to ClimaSeaIce (Kelvin → Celsius, only where ice exists)
    csf.scalar_temp2 .=
        ifelse.(
            area_fraction .≈ 0,
            zero(FT),
            fluxes.T_sfc_new .- FT(sim.ice_properties.C_to_K),
        )
    CC.Remapping.interpolate!(
        sim.remapping.scratch_arr1,
        sim.remapping.remapper_cc,
        csf.scalar_temp2,
    )
    ice_concentration = sim.ice.model.ice_concentration
    top_sfc_T = sim.ice.model.ice_thermodynamics.top_surface_temperature
    OC.interior(top_sfc_T, :, :, 1) .=
        ifelse.(
            OC.interior(ice_concentration, :, :, 1) .> 0,
            sim.remapping.scratch_arr1,
            OC.interior(top_sfc_T, :, :, 1),
        )

    return nothing
end

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

    # We only need to provide momentum fluxes if the sea ice model has dynamics
    if !isnothing(sim.ice.model.dynamics)
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
        F_turb_ρτxz_cell = sim.remapping.scratch_cc1
        F_turb_ρτyz_cell = sim.remapping.scratch_cc2

        # Set the momentum flux BCs at the correct locations using the remapped scratch fields
        # Note that this requires the sea ice model to always be run with dynamics turned on
        si_flux_u = sim.ice.model.dynamics.external_momentum_stresses.top.u
        si_flux_v = sim.ice.model.dynamics.external_momentum_stresses.top.v
        set_from_extrinsic_vector!(
            (; u = si_flux_u, v = si_flux_v),
            grid,
            F_turb_ρτxz_cell,
            F_turb_ρτyz_cell,
        )
    end

    # Remap the latent and sensible heat fluxes using scratch arrays
    CC.Remapping.interpolate!(sim.remapping.scratch_arr1, sim.remapping.remapper_cc, F_lh) # latent heat flux
    CC.Remapping.interpolate!(sim.remapping.scratch_arr2, sim.remapping.remapper_cc, F_sh) # sensible heat flux

    # Rename for clarity; recall F_turb_energy = F_lh + F_sh
    remapped_F_lh = sim.remapping.scratch_arr1
    remapped_F_sh = sim.remapping.scratch_arr2

    # Update the sea ice heat flux only where the concentration is greater than zero.
    # With PrescribedTemperature the top heat flux is a FluxFunction, not a Field;
    # the flux is determined from the diagnosed T_sfc so we skip writing here.
    si_flux_heat = sim.ice.model.external_heat_fluxes.top
    if si_flux_heat isa OC.Field
        OC.interior(si_flux_heat, :, :, 1) .+=
            (OC.interior(ice_concentration, :, :, 1) .> 0) .*
            (remapped_F_lh .+ remapped_F_sh)
    end

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
    # With PrescribedTemperature the top heat flux is a FluxFunction, not a Field;
    # the flux is determined from the diagnosed T_sfc so we skip writing here.
    si_flux_heat = sim.ice.model.external_heat_fluxes.top
    if si_flux_heat isa OC.Field
        α = Interfacer.get_field(sim, Val(:surface_direct_albedo)) # scalar
        ϵ = Interfacer.get_field(sim, Val(:emissivity)) # scalar

        OC.interior(si_flux_heat, :, :, 1) .=
            (OC.interior(ice_concentration, :, :, 1) .> 0) .*
            (-(1 .- α) .* remapped_SW_d .- ϵ .* remapped_LW_d)
    end
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
    grid = ocean_sim.ocean.model.grid
    ocean_properties = ocean_sim.ocean_properties
    ice_concentration = Interfacer.get_field(ice_sim, Val(:ice_concentration))

    # Update the sea ice concentration in the ocean simulation
    ocean_sim.ice_concentration .= ice_concentration

    # Compute the fluxes and store them in the both simulations
    CO.OceanSeaIceModels.InterfaceComputations.compute_sea_ice_ocean_fluxes!(
        ice_sim.ocean_ice_interface,
        ocean_sim.ocean,
        ice_sim.ice,
        ocean_properties,
    )

    ## Update the internals of the sea ice model
    # Set the bottom heat flux to the sum of the frazil and interface heat fluxes
    bottom_heat_flux = ice_sim.ice.model.external_heat_fluxes.bottom

    Qf = ice_sim.ocean_ice_interface.fluxes.frazil_heat        # frazil heat flux
    Qi = ice_sim.ocean_ice_interface.fluxes.interface_heat     # interfacial heat flux
    bottom_heat_flux .= Qf .+ Qi

    ## Update the internals of the ocean model
    ρₒ⁻¹ = 1 / ocean_sim.ocean_properties.reference_density
    cₒ = ocean_sim.ocean_properties.heat_capacity

    # Compute fluxes for u, v, T, and S from momentum, heat, and freshwater fluxes
    oc_flux_u = surface_flux(ocean_sim.ocean.model.velocities.u)
    oc_flux_v = surface_flux(ocean_sim.ocean.model.velocities.v)

    ρτxio = ice_sim.ocean_ice_interface.fluxes.x_momentum # sea_ice - ocean zonal momentum flux
    ρτyio = ice_sim.ocean_ice_interface.fluxes.y_momentum # sea_ice - ocean meridional momentum flux

    # mask out the poles
    polar_mask = ocean_sim.remapping.polar_mask
    @. ρτxio = polar_mask * ρτxio
    @. ρτyio = polar_mask * ρτyio

    # Update the momentum flux contributions from ocean/sea ice fluxes
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

    # The heat and salt fluxes already include the SIC masking, so we don't need to
    # multiply by SIC here.
    oc_flux_T = surface_flux(ocean_sim.ocean.model.tracers.T)
    heat_flux = OC.interior(ocean_sim.remapping.scratch_cc1, :, :, 1)
    heat_flux .= OC.interior(Qi, :, :, 1) .* ρₒ⁻¹ ./ cₒ
    # mask out the poles
    @. heat_flux = polar_mask * heat_flux
    OC.interior(oc_flux_T, :, :, 1) .+= heat_flux

    oc_flux_S = surface_flux(ocean_sim.ocean.model.tracers.S)
    salt_contrib = OC.interior(ocean_sim.remapping.scratch_cc2, :, :, 1)
    salt_contrib .= OC.interior(ice_sim.ocean_ice_interface.fluxes.salt, :, :, 1)
    # mask out the poles
    @. salt_contrib = polar_mask * salt_contrib
    OC.interior(oc_flux_S, :, :, 1) .+= salt_contrib

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

# Additional ClimaSeaIceSimulation getter methods for plotting debug fields
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:u}) = sim.ice.model.velocities.u
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:v}) = sim.ice.model.velocities.v

"""
    Plotting.debug_plot_fields(sim::ClimaSeaIceSimulation)

Return the fields to include in debug plots for a ClimaSeaIce simulation.
This includes the area fraction, surface temperature, ice concentration, ice thickness, and
zonal and meridional velocity fields. Note that if the sea ice model does not have dynamics,
the velocity fields will be zero.

These plots are not polished, and are intended for debugging.
"""
Plotting.debug_plot_fields(sim::ClimaSeaIceSimulation) =
    (:area_fraction, :surface_temperature, :ice_concentration, :ice_thickness, :u, :v)
