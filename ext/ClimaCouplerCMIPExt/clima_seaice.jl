using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions:
    IceWaterThermalEquilibrium, MeltingConstrainedFluxBalance, PrescribedTemperature, get_tracer, RadiativeEmission
import ClimaComms
import ClimaOcean.EN4: download_dataset
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
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
struct ClimaSeaIceSimulation{SIM, A, REMAP, NT, IP} <: Interfacer.AbstractSeaIceSimulation
    ice::SIM
    area_fraction::A
    remapping::REMAP
    ocean_ice_interface::NT
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
    FT;
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
    Interfacer.SeaIceSimulation(::Type{FT}, ::Val{:clima_seaice}; kwargs...)

Extension of the generic SeaIceSimulation constructor for ClimaSeaIce.
"""
function Interfacer.SeaIceSimulation(::Type{FT}, ::Val{:clima_seaice}; kwargs...) where {FT}
    return ClimaSeaIceSimulation(FT; kwargs...)
end

"""
    ClimaSeaIceSimulation()

Creates an ClimaSeaIceSimulation object containing a model, an integrator, and
a surface area fraction field.

If a start date is provided, we initialize the sea ice concentration and thickness
using the ECCO4Monthly dataset. If no start date is provided, we initialize with zero sea ice.

Since the coupler iteratively diagnoses T_sfc to satisfy the flux balance
(via the `update_T_sfc` callback in `compute_surface_fluxes!`), we use
`PrescribedTemperature` as the top heat boundary condition. The coupler
writes the diagnosed T_sfc back to ClimaSeaIce's `top_surface_temperature`
field after each flux computation, so the ice thermodynamics always uses
the flux-balanced surface temperature.

Specific details about the default model configuration
can be found in the documentation for `ClimaSeaIce.sea_ice_simulation`.

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
    ice = sea_ice_simulation(grid, ocean.ocean; Δt = dt, advection)

    ocean_ice_flux_formulation =
        CO.OceanSeaIceModels.InterfaceComputations.ThreeEquationHeatFlux(ice)
    interface_temperature = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    interface_salinity = OC.Field{OC.Center, OC.Center, Nothing}(grid)

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
    )

    # Ensure ocean temperature is above freezing where there is sea ice
    CO.OceanSeaIceModels.above_freezing_ocean_temperature!(ocean.ocean, grid, ice)
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
    dynamics = CO.SeaIces.sea_ice_dynamics(grid, ocean),
    phase_transitions = CSI.PhaseTransitions(; ice_heat_capacity, ice_density),
    conductivity = 2, # kg m s⁻³ K⁻¹
    internal_heat_flux = CSI.ConductiveFlux(; conductivity),
)
    FT = eltype(grid)

    # Initialize top_surface_temperature as a mutable Field (not a ConstantField)
    # so it can be remapped and updated by the coupler with the diagnosed T_sfc.
    top_sfc_temp_init = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    OC.set!(top_sfc_temp_init, FT(0)) # 0°C initial guess; overwritten by coupler
    top_heat_boundary_condition = PrescribedTemperature(top_sfc_temp_init)
    kᴺ = size(grid, 3)
    surface_ocean_salinity = OC.interior(ocean.model.tracers.S, :, :, (kᴺ:kᴺ))
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
### Functions required by ClimaCoupler.jl for a AbstractSurfaceSimulation
###############################################################################

# Timestep the simulation forward to time `t`
Interfacer.step!(sim::ClimaSeaIceSimulation, t) =
    OC.time_step!(sim.ice, float(t) - sim.ice.model.clock.time)

Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:area_fraction}) = sim.area_fraction
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:ice_concentration}) =
    sim.ice.model.ice_concentration

# Get ice thickness field (remapped to boundary space)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:ice_thickness}) =
    Interfacer.remap(axes(sim.area_fraction), sim.ice.model.ice_thickness)

# Get internal temperature (interface temperature with ocean, remapped to boundary space)
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:internal_temperature}) =
    Interfacer.remap(axes(sim.area_fraction), sim.ocean_ice_interface.temperature)

# TODO better values for this
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

# Approximate the sea ice surface temperature as the temperature computed from the
#  fluxes at the previous timestep.
Interfacer.get_field(sim::ClimaSeaIceSimulation, ::Val{:surface_temperature}) =
    sim.ice_properties.C_to_K + sim.ice.model.ice_thermodynamics.top_surface_temperature

"""
    _get_surface_fluxes_seaice(args...)

Wrapper around `SF.surface_fluxes` that extracts scalar fields from the
`SurfaceFluxConditions` struct, returning a NamedTuple instead.

This is necessary because broadcasting `SF.surface_fluxes.()` over ClimaCore
fields would try to store `SurfaceFluxConditions` in a ClimaCore field, which
fails because ClimaCore can only store flat numeric types, not nested structs.

This follows the same pattern as `FluxCalculator.get_surface_fluxes` in
`FluxCalculator.jl`, but includes the extra callback parameters needed for
the `update_T_sfc` skin temperature solver.
"""
function _get_surface_fluxes_seaice(
    surface_fluxes_params,
    T_int,
    q_tot_int,
    q_liq_int,
    q_ice_int,
    ρ_int,
    T_sfc,
    q_vap_sfc,
    Φ_sfc,
    Δz,
    d,
    u_int,
    u_sfc,
    roughness_inputs,
    config,
    scheme,
    flux_specs,
    update_q_vap_sfc,
    update_T_sfc_cb,
    update_q_vap_sfc2,
)
    outputs = SF.surface_fluxes(
        surface_fluxes_params,
        T_int,
        q_tot_int,
        q_liq_int,
        q_ice_int,
        ρ_int,
        T_sfc,
        q_vap_sfc,
        Φ_sfc,
        Δz,
        d,
        u_int,
        u_sfc,
        roughness_inputs,
        config,
        scheme,
        flux_specs,
        update_q_vap_sfc,
        update_T_sfc_cb,
        update_q_vap_sfc2,
    )

    (; shf, lhf, evaporation, ρτxz, ρτyz, T_sfc, q_vap_sfc, L_MO, ustar) = outputs

    buoyancy_flux = SF.buoyancy_flux(
        surface_fluxes_params,
        shf,
        lhf,
        T_sfc,
        ρ_int,
        q_vap_sfc,
        q_liq_int,
        q_ice_int,
        SF.MoistModel(),
    )

    return (;
        F_turb_ρτxz = ρτxz,
        F_turb_ρτyz = ρτyz,
        F_sh = shf,
        F_lh = lhf,
        F_turb_moisture = evaporation,
        L_MO,
        ustar,
        buoyancy_flux,
        T_sfc_new = T_sfc,
    )
end

"""
    FluxCalculator.compute_surface_fluxes!(csf, sim::ClimaSeaIceSimulation, atmos_sim, thermo_params)

Compute surface fluxes for ClimaSeaIceSimulation using the `update_T_sfc` callback
to iteratively diagnose the surface temperature that satisfies the flux balance equation.

This function extends the default `compute_surface_fluxes!` to use the skin temperature
update callback from `ClimaCouplerCMIPExt.update_T_sfc`.

After computing fluxes, the diagnosed surface temperature is written back to
ClimaSeaIce's `top_surface_temperature` field. Since the ice thermodynamics uses
`PrescribedTemperature`, it takes this field as-is during its time step, ensuring
the ice model evolves consistently with the flux-balanced surface temperature
diagnosed by the coupler.
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

    # Atmosphere fields are stored in coupler fields so we only regrid them once per timestep
    uv_int = StaticArrays.SVector.(csf.u_int, csf.v_int)

    # Get initial surface temperature guess
    Interfacer.get_field!(csf.scalar_temp1, sim, Val(:surface_temperature))
    T_sfc = csf.scalar_temp1

    # TODO: CHECKING IF THIS IS COMING FROM THE REGRIDDING THROUGH MAP_INTERPOLATE!
    # REMOVE THIS REMAPPER +/- 80 degree fix! 
    @. T_sfc = ifelse(T_sfc < FT(100), csf.T_atmos, T_sfc)

    # Compute surface humidity from the surface temperature
    ρ_sfc =
        SF.surface_density.(
            surface_fluxes_params,
            csf.T_atmos,
            csf.ρ_atmos,
            T_sfc,
            csf.height_int .- csf.height_sfc,
            csf.q_tot_atmos,
            0, # q_liq
            0, # q_ice
        )

    csf.scalar_temp2 .= TD.q_vap_saturation.(thermo_params, T_sfc, ρ_sfc, 0, 0)
    q_sfc = csf.scalar_temp2

    # Set gustiness
    gustiness = ones(boundary_space)

    # Get roughness model (sea ice uses constant roughness)
    roughness_model = Interfacer.get_field(sim, Val(:roughness_model))
    roughness_params = if roughness_model == :coare3
        Interfacer.get_field(sim, Val(:coare3_roughness_params))
    elseif roughness_model == :constant
        Interfacer.get_field!(csf.scalar_temp3, sim, Val(:roughness_momentum))
        z0m = csf.scalar_temp3
        Interfacer.get_field!(csf.scalar_temp4, sim, Val(:roughness_buoyancy))
        z0b = csf.scalar_temp4
        SF.ConstantRoughnessParams.(z0m, z0b)
    else
        error("Unknown roughness_model: $roughness_model. Must be :coare3 or :constant")
    end

    config = SF.SurfaceFluxConfig.(roughness_params, SF.ConstantGustinessSpec.(gustiness))

    # Get sea ice parameters for update_T_sfc callback
    # Thermal conductivity (scalar); ConductiveFlux has .conductivity, FluxFunction does not
    internal_heat_flux = sim.ice.model.ice_thermodynamics.internal_heat_flux
    κ = if hasfield(typeof(internal_heat_flux), :conductivity)
        FT.(internal_heat_flux.conductivity)
    else
        convert(FT, 2) # default conductivity [W m⁻¹ K⁻¹] when internal_heat_flux is FluxFunction etc.
    end
    
    # Ice thickness (field, remapped to boundary space)
    δ = FT.(Interfacer.get_field(sim, Val(:ice_thickness)))
    # Internal temperature (field, remapped to boundary space).
    # Oceananigans uses Celsius; convert to Kelvin to match the coupler/SurfaceFluxes convention.
    # Guard against 0 from map_interpolate! at out-of-grid points (same issue as T_sfc).
    T_i_raw = FT.(Interfacer.get_field(sim, Val(:internal_temperature))) .+ FT(sim.ice_properties.C_to_K)
    #TODO: CHECKING IF THIS IS COMING FROM THE REGRIDDING THROUGH MAP_INTERPOLATE!
    T_i = @. ifelse(T_i_raw < FT(100), csf.T_atmos, T_i_raw)
    # Stefan-Boltzmann constant (scalar)
    σ = FT(sim.ice_properties.σ)
    # Emissivity (scalar, broadcast to field)
    ϵ = FT.(Interfacer.get_field(sim, Val(:emissivity)))
    
    # Radiation and albedo (fields)
    SW_d = csf.SW_d
    LW_d = csf.LW_d
    α_albedo = FT.(Interfacer.get_field(sim, Val(:surface_direct_albedo)))

    # Create the update_T_sfc callback element-wise
    # Since update_T_sfc returns a function, we broadcast it to create a field of callbacks
    update_T_sfc_callback = ClimaCouplerCMIPExt.update_T_sfc.(
        κ, # conductivity
        δ, # ice thickness
        T_i, # internal temperature
        σ, # stefan-boltzmann constant
        ϵ, # emissivity
        SW_d, #SW↓
        LW_d, #LW↓
        α_albedo, # albedo
    )

    # Calculate surface fluxes with the callback
    Φ_sfc = SFP.grav.(surface_fluxes_params) .* csf.height_sfc
    Δz = csf.height_int .- csf.height_sfc

    # Use a wrapper that calls SF.surface_fluxes and extracts scalar fields from
    # SurfaceFluxConditions struct. Broadcasting SF.surface_fluxes.() directly would
    # try to store SurfaceFluxConditions in a ClimaCore field, which fails because
    # ClimaCore can only store flat numeric types, not nested structs.
    fluxes = _get_surface_fluxes_seaice.(
        surface_fluxes_params,
        csf.T_atmos,
        csf.q_tot_atmos,
        csf.q_liq_atmos,
        csf.q_ice_atmos,
        csf.ρ_atmos,
        T_sfc,
        q_sfc,
        Φ_sfc,
        Δz,
        FT(0), # d
        uv_int,
        uv_int .* FT(0),
        nothing, 
        config,
        Ref(SF.PointValueScheme()),
        nothing, 
        nothing,
        update_T_sfc_callback,
        nothing, # update_q_vap_sfc (not used)
    )

    (; F_turb_ρτxz, F_turb_ρτyz, F_sh, F_lh, F_turb_moisture, L_MO, ustar, buoyancy_flux, T_sfc_new) = fluxes

    # Get area fraction
    Interfacer.get_field!(csf.scalar_temp1, sim, Val(:area_fraction))
    area_fraction = csf.scalar_temp1

    # Zero out fluxes where the area fraction is zero
    @. F_turb_ρτxz = ifelse(area_fraction ≈ 0, zero(F_turb_ρτxz), F_turb_ρτxz)
    @. F_turb_ρτyz = ifelse(area_fraction ≈ 0, zero(F_turb_ρτyz), F_turb_ρτyz)
    @. F_sh = ifelse(area_fraction ≈ 0, zero(F_sh), F_sh)
    @. F_lh = ifelse(area_fraction ≈ 0, zero(F_lh), F_lh)
    @. F_turb_moisture = ifelse(area_fraction ≈ 0, zero(F_turb_moisture), F_turb_moisture)
    @. L_MO = ifelse(area_fraction ≈ 0, zero(L_MO), L_MO)
    @. ustar = ifelse(area_fraction ≈ 0, zero(ustar), ustar)
    @. buoyancy_flux = ifelse(area_fraction ≈ 0, zero(buoyancy_flux), buoyancy_flux)

    # Update the fluxes in the simulation
    fields = (; F_turb_ρτxz, F_turb_ρτyz, F_lh, F_sh, F_turb_moisture)
    FluxCalculator.update_turbulent_fluxes!(sim, fields)

    # Update fluxes in the coupler fields (area-weighted)
    @. csf.F_turb_ρτxz += F_turb_ρτxz * area_fraction
    @. csf.F_turb_ρτyz += F_turb_ρτyz * area_fraction
    @. csf.F_lh += F_lh * area_fraction
    @. csf.F_sh += F_sh * area_fraction
    @. csf.F_turb_moisture += F_turb_moisture * area_fraction

    # Handle L_MO separately (can be Inf)
    @. csf.L_MO += ifelse(isinf(L_MO), L_MO, L_MO * area_fraction)
    @. csf.ustar += ustar * area_fraction
    @. csf.buoyancy_flux += buoyancy_flux * area_fraction

    # Write the diagnosed T_sfc back to ClimaSeaIce's top_surface_temperature.
    # With PrescribedTemperature, ClimaSeaIce uses this field as-is during its
    # time step, so updating it ensures the ice thermodynamics uses the T_sfc
    # that satisfies the flux balance diagnosed by the coupler.
    # T_sfc_new is in Kelvin (coupler convention); ClimaSeaIce uses Celsius.
    csf.scalar_temp2 .= ifelse.(
        area_fraction .≈ 0,
        zero(FT),
        T_sfc_new .- FT(sim.ice_properties.C_to_K),
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
    si_flux_heat = sim.ice.model.external_heat_fluxes.top
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
    si_flux_heat = sim.ice.model.external_heat_fluxes.top
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
        OC.interior(ice_sim.ocean_ice_interface.fluxes.salt, :, :, 1)

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
