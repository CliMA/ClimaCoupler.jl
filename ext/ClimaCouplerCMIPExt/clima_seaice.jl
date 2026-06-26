using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions:
    IceWaterThermalEquilibrium, PrescribedTemperature, get_tracer
import ClimaComms
import ClimaOcean.EN4: download_dataset
import NVTX
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Dates
import ClimaUtilities.TimeManager: ITime, date, counter, period
import ClimaUtilities.ClimaArtifacts: @clima_artifact
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
- `ice_properties::IP`: A NamedTuple of sea ice properties: Stefan–Boltzmann constant `σ`
    and Celsius offset `C_to_K` (water freezing point in K).
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
diagnoses T_sfc via the `update_T_sfc_cb` closure passed to `FluxCalculator.get_surface_fluxes!`
and writes it back to the ice model after each flux computation (see `compute_surface_fluxes!`).

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
    seaice_diagnostic_interval = "1days",
    seaice_diagnostic_mode = :average,
    extra_kwargs...,
) where {FT}
    # Initialize the sea ice with the same grid as the ocean
    grid = ocean.ocean.model.grid
    arch = OC.Architectures.architecture(grid)

    advection = ocean.ocean.model.advection.T
    ice = CO.SeaIces.sea_ice_simulation(grid, ocean.ocean; Δt = float(dt), advection)

    ocean_ice_flux_formulation = CO.InterfaceComputations.ThreeEquationHeatFlux(ice)
    interface_temperature = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    interface_salinity = OC.Field{OC.Center, OC.Center, Nothing}(grid)

    # Initialize model_Δt so that time stepping works properly
    model_Δt = dt

    # Initialize nonzero sea ice if start date provided
    if !isnothing(start_date)
        # set up the `dir` keyword argument for `Metadatum`
        if start_date == Dates.Date(2010, 1, 1)
            # we have a ClimaArtifact saved for January 1, 2010 (so that CI can always run)
            dir_kw = (; dir = @clima_artifact("ecco4_SIarea_SIheff_2010_01"))
            @info "Using $(dir_kw.dir) ClimaArtifact for sea ice initialization on $(start_date)"
        else
            # otherwise, download the data
            # (or load from scratchspace; ClimaOcean will automatically handle this)
            dir_kw = (;)
        end

        sic_metadata = CO.DataWrangling.Metadatum(
            :sea_ice_concentration,
            dataset = CO.DataWrangling.ECCO.ECCO4Monthly(),
            date = start_date;
            dir_kw...,
        )
        h_metadata = CO.DataWrangling.Metadatum(
            :sea_ice_thickness,
            dataset = CO.DataWrangling.ECCO.ECCO4Monthly(),
            date = start_date;
            dir_kw...,
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
    area_fraction = Interfacer.remap(boundary_space, ice.model.ice_concentration, remapping)

    sim = ClimaSeaIceSimulation(
        ice,
        area_fraction,
        remapping,
        ocean_ice_interface,
        ice_properties,
        model_Δt,
    )

    add_seaice_diagnostics!(
        sim;
        output_dir,
        interval = TimeManager.time_to_period(seaice_diagnostic_interval),
        mode = seaice_diagnostic_mode,
    )

    # Ensure ocean temperature is above freezing where there is sea ice
    CO.EarthSystemModels.above_freezing_ocean_temperature!(ocean.ocean, grid, ice)
    return sim
end

###############################################################################
### Functions required by ClimaCoupler.jl for a AbstractSurfaceSimulation
###############################################################################

# Timestep the simulation forward to time `t`
NVTX.@annotate function Interfacer.step!(sim::ClimaSeaIceSimulation, t::Float64)
    # `round(Int, ...)` tolerates floating point drift less than `model_dt / 2`
    n_steps = round(Int, (t - sim.ice.model.clock.time) / sim.model_Δt)
    for _ in 1:n_steps
        OC.time_step!(sim.ice, sim.model_Δt)
    end
    return nothing
end

NVTX.@annotate function Interfacer.step!(sim::ClimaSeaIceSimulation, t::ITime)
    Δt_msec = date(t) - sim.ice.model.clock.time
    model_Δt_msec = counter(sim.model_Δt) * Dates.Millisecond(period(sim.model_Δt))
    n_steps = div(Δt_msec, model_Δt_msec) # integer division; exact for Millisecond periods
    for _ in 1:n_steps
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
    extract_ice_balance_cc_nodal!(balance_nodal, csf, sim, boundary_space)

Flatten ice skin-balance inputs from boundary-space fields into nodal vectors
for polygon averaging on the intersection grid.
"""
function extract_ice_balance_cc_nodal!(balance_nodal, csf, sim, boundary_space)
    CRExt = get_ConservativeRegriddingCCExt()
    scratch = CC.Fields.zeros(boundary_space)

    balance_nodal.SW_d .= CRExt.se_field_to_vec(csf.SW_d)
    balance_nodal.LW_d .= CRExt.se_field_to_vec(csf.LW_d)

    Interfacer.get_field!(scratch, sim, Val(:ice_thickness))
    balance_nodal.δ .= CRExt.se_field_to_vec(scratch)
    Interfacer.get_field!(scratch, sim, Val(:internal_temperature))
    balance_nodal.T_i .= CRExt.se_field_to_vec(scratch)
    Interfacer.get_field!(scratch, sim, Val(:emissivity))
    balance_nodal.ϵ .= CRExt.se_field_to_vec(scratch)
    Interfacer.get_field!(scratch, sim, Val(:surface_direct_albedo))
    balance_nodal.α_albedo .= CRExt.se_field_to_vec(scratch)
    return nothing
end

"""
    extract_ice_surface_state_for_intersection!(ice_surface_temp, sim)

Extract sea-ice surface state onto the flat OC-cell layout used by the
intersection grid.
"""
function extract_ice_surface_state_for_intersection!(ice_surface_temp, sim::ClimaSeaIceSimulation)
    C_to_K = sim.ice_properties.C_to_K
    top_T = sim.ice.model.ice_thermodynamics.top_surface_temperature
    ice_surface_temp.T .= vec(OC.interior(top_T, :, :, 1)) .+ C_to_K

    FT = eltype(ice_surface_temp.T)
    z0m = Interfacer.get_field(sim, Val(:roughness_momentum))
    z0b = Interfacer.get_field(sim, Val(:roughness_buoyancy))
    fill!(ice_surface_temp.z0m, FT(z0m))
    fill!(ice_surface_temp.z0b, FT(z0b))
    fill!(ice_surface_temp.h, FT(0))
    return nothing
end

"""
    compute_ice_intersection_grid_fluxes!(sim, csf, surface_fluxes_params, thermo_params)

Compute sea-ice–atmosphere fluxes on the CC×OC intersection grid with a
per-polygon [`update_T_sfc`](@ref) skin-temperature balance.
"""
function compute_ice_intersection_grid_fluxes!(
    sim::ClimaSeaIceSimulation,
    csf,
    surface_fluxes_params,
    thermo_params,
)
    boundary_space = axes(csf)
    FT = CC.Spaces.undertype(boundary_space)
    (
        intersection_grid,
        ice_intersection_flux_state,
        ice_cc_balance_nodal,
        ice_balance_at_int,
        ice_surface_temp,
        cc_atmos_temp,
    ) = sim.remapping

    extract_cc_atmos_state!(cc_atmos_temp, csf, boundary_space)
    cc_atmos_state = (
        T = cc_atmos_temp.T,
        q_tot = cc_atmos_temp.q_tot,
        q_liq = cc_atmos_temp.q_liq,
        q_ice = cc_atmos_temp.q_ice,
        ρ = cc_atmos_temp.ρ,
        u = cc_atmos_temp.u,
        v = cc_atmos_temp.v,
        h = cc_atmos_temp.h,
    )

    extract_ice_balance_cc_nodal!(ice_cc_balance_nodal, csf, sim, boundary_space)

    internal_heat_flux = sim.ice.model.ice_thermodynamics.internal_heat_flux
    κ = if hasfield(typeof(internal_heat_flux), :conductivity)
        FT(internal_heat_flux.conductivity)
    else
        FT(2)
    end
    gather_ice_balance_to_intersection!(
        ice_balance_at_int,
        intersection_grid,
        ice_cc_balance_nodal,
        κ,
    )

    extract_ice_surface_state_for_intersection!(ice_surface_temp, sim)
    ice_surface_state = (
        T = ice_surface_temp.T,
        z0m = ice_surface_temp.z0m,
        z0b = ice_surface_temp.z0b,
        h = ice_surface_temp.h,
    )

    sic_oc = vec(OC.interior(sim.ice.model.ice_concentration, :, :, 1))
    ice_active = ice_intersection_active_mask(intersection_grid, sic_oc)

    compute_surface_fluxes_on_intersection!(
        ice_intersection_flux_state,
        intersection_grid,
        cc_atmos_state,
        ice_surface_state,
        surface_fluxes_params,
        thermo_params;
        ice_balance_at_int,
        σ = FT(sim.ice_properties.σ),
        T_melt = FT(sim.ice_properties.C_to_K),
        ice_active,
    )
    return nothing
end

"""
    write_diagnosed_ice_T_sfc_to_model!(sim::ClimaSeaIceSimulation)

Scatter polygon-diagnosed surface temperatures back to the ice model's
`top_surface_temperature` field (Celsius) on cells with `SIC > 0`.
"""
function write_diagnosed_ice_T_sfc_to_model!(sim::ClimaSeaIceSimulation)
    (; intersection_grid, ice_intersection_flux_state) = sim.remapping
    grid = sim.ice.model.grid
    Nx_oc, Ny_oc, _ = size(grid)
    n_oc_layout = Nx_oc * Ny_oc
    FT = eltype(ice_intersection_flux_state.surface_T)

    T_oc = OC.on_architecture(OC.architecture(grid), zeros(FT, n_oc_layout))
    scatter_to_oc!(T_oc, intersection_grid, ice_intersection_flux_state.surface_T)
    T_oc_2d = reshape(T_oc, Nx_oc, Ny_oc)

    C_to_K = sim.ice_properties.C_to_K
    ice_concentration = OC.interior(sim.ice.model.ice_concentration, :, :, 1)
    top_sfc_T = sim.ice.model.ice_thermodynamics.top_surface_temperature
    OC.interior(top_sfc_T, :, :, 1) .= ifelse.(
        ice_concentration .> 0,
        T_oc_2d .- C_to_K,
        OC.interior(top_sfc_T, :, :, 1),
    )
    return nothing
end

"""
    push_intersection_fluxes_to_ice!(sim::ClimaSeaIceSimulation)

Scatter precomputed ice intersection-grid fluxes to ClimaSeaIce boundary
conditions. Assumes `compute_ice_intersection_grid_fluxes!` has already run.
"""
function push_intersection_fluxes_to_ice!(sim::ClimaSeaIceSimulation)
    (; intersection_grid, ice_intersection_flux_state) = sim.remapping
    grid = sim.ice.model.grid
    ice_concentration = sim.ice.model.ice_concentration

    Nx_oc, Ny_oc, _ = size(grid)
    n_oc_layout = Nx_oc * Ny_oc
    FT = eltype(ice_intersection_flux_state.flux_sh)

    F_sh_oc = OC.on_architecture(OC.architecture(grid), zeros(FT, n_oc_layout))
    F_lh_oc = OC.on_architecture(OC.architecture(grid), zeros(FT, n_oc_layout))
    F_τx_oc = OC.on_architecture(OC.architecture(grid), zeros(FT, n_oc_layout))
    F_τy_oc = OC.on_architecture(OC.architecture(grid), zeros(FT, n_oc_layout))
    scatter_to_oc!(F_sh_oc, intersection_grid, ice_intersection_flux_state.flux_sh)
    scatter_to_oc!(F_lh_oc, intersection_grid, ice_intersection_flux_state.flux_lh)
    scatter_to_oc!(F_τx_oc, intersection_grid, ice_intersection_flux_state.flux_τx)
    scatter_to_oc!(F_τy_oc, intersection_grid, ice_intersection_flux_state.flux_τy)

    F_sh_2d = reshape(F_sh_oc, Nx_oc, Ny_oc)
    F_lh_2d = reshape(F_lh_oc, Nx_oc, Ny_oc)
    F_τx_2d = reshape(F_τx_oc, Nx_oc, Ny_oc)
    F_τy_2d = reshape(F_τy_oc, Nx_oc, Ny_oc)
    has_ice = OC.interior(ice_concentration, :, :, 1) .> 0

    if !isnothing(sim.ice.model.dynamics)
        OC.interior(sim.remapping.scratch_field_oc1, :, :, 1) .= F_τx_2d
        OC.interior(sim.remapping.scratch_field_oc2, :, :, 1) .= F_τy_2d
        si_flux_u = sim.ice.model.dynamics.external_momentum_stresses.top.u
        si_flux_v = sim.ice.model.dynamics.external_momentum_stresses.top.v
        set_from_extrinsic_vector!(
            (; u = si_flux_u, v = si_flux_v),
            grid,
            sim.remapping.scratch_field_oc1,
            sim.remapping.scratch_field_oc2,
        )
    end

    si_flux_heat = sim.ice.model.external_heat_fluxes.top
    if si_flux_heat isa OC.Field
        OC.interior(si_flux_heat, :, :, 1) .+=
            has_ice .* (F_lh_2d .+ F_sh_2d)
    end
    return nothing
end

"""
    FluxCalculator.compute_surface_fluxes!(csf, sim::ClimaSeaIceSimulation, atmos_sim, thermo_params)

Compute surface fluxes for `ClimaSeaIceSimulation` on the CC×OC intersection grid,
diagnosing `T_sfc` per polygon via [`update_T_sfc`](@ref).
"""
NVTX.@annotate function FluxCalculator.compute_surface_fluxes!(
    csf,
    sim::ClimaSeaIceSimulation,
    atmos_sim::Interfacer.AbstractAtmosSimulation,
    thermo_params,
    accumulator = nothing,
)
    boundary_space = axes(csf)
    surface_fluxes_params = FluxCalculator.get_surface_params(atmos_sim)

    compute_ice_intersection_grid_fluxes!(sim, csf, surface_fluxes_params, thermo_params)

    fluxes = intersection_fluxes_to_boundary_fields(
        boundary_space,
        sim.remapping.intersection_grid,
        sim.remapping.ice_intersection_flux_state,
    )

    FluxCalculator.update_flux_fields!(csf, sim, fluxes, accumulator)
    write_diagnosed_ice_T_sfc_to_model!(sim)
    return nothing
end

"""
    FluxCalculator.update_turbulent_fluxes!(sim::ClimaSeaIceSimulation, fields)

Push intersection-grid turbulent fluxes from `ice_intersection_flux_state` to the
ice model boundary conditions.
"""
function FluxCalculator.update_turbulent_fluxes!(sim::ClimaSeaIceSimulation, fields::NamedTuple)
    push_intersection_fluxes_to_ice!(sim)
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

    # Remap radiative fluxes onto scratch fields (separate buffers so SW is not overwritten by LW)
    Interfacer.remap!(sim.remapping.scratch_field_oc1, csf.SW_d, sim.remapping) # shortwave radiation
    remapped_SW_d = OC.interior(sim.remapping.scratch_field_oc1, :, :, 1)

    Interfacer.remap!(sim.remapping.scratch_field_oc2, csf.LW_d, sim.remapping) # longwave radiation
    remapped_LW_d = OC.interior(sim.remapping.scratch_field_oc2, :, :, 1)

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
    CO.InterfaceComputations.compute_sea_ice_ocean_fluxes!(
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
    heat_flux = OC.interior(ocean_sim.remapping.scratch_field_oc1, :, :, 1)
    heat_flux .= OC.interior(Qi, :, :, 1) .* ρₒ⁻¹ ./ cₒ
    OC.interior(oc_flux_T, :, :, 1) .+= heat_flux

    oc_flux_S = surface_flux(ocean_sim.ocean.model.tracers.S)
    salt_contrib = OC.interior(ocean_sim.remapping.scratch_field_oc2, :, :, 1)
    salt_contrib .= OC.interior(ice_sim.ocean_ice_interface.fluxes.salt, :, :, 1)
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

# Extend Interfacer.get_field to allow automatic remapping to the target space
# (defined here so `OceananigansSimulation` and `ClimaSeaIceSimulation` exist.)
function Interfacer.get_field!(
    target_field,
    sim::Union{OceananigansSimulation, ClimaSeaIceSimulation},
    quantity,
)
    Interfacer.remap!(target_field, Interfacer.get_field(sim, quantity), sim.remapping)
    return nothing
end
# TODO see if we can remove this allocating version
function Interfacer.get_field(
    target_space::CC.Spaces.AbstractSpace,
    sim::Union{OceananigansSimulation, ClimaSeaIceSimulation},
    quantity,
)
    return Interfacer.remap(
        target_space,
        Interfacer.get_field(sim, quantity),
        sim.remapping,
    )
end
