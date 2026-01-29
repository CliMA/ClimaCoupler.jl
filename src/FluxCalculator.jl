"""
    FluxCalculator

This modules contains abstract types and functions to calculate surface fluxes in the coupler,
or to call flux calculating functions from the component models.
"""
module FluxCalculator

import StaticArrays
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
import ClimaCore as CC
import ..Interfacer, ..Utilities

export turbulent_fluxes!,
    get_surface_params,
    update_turbulent_fluxes!,
    compute_surface_fluxes!,
    ocean_seaice_fluxes!

function turbulent_fluxes!(cs::Interfacer.CoupledSimulation)
    return turbulent_fluxes!(cs.fields, cs.model_sims, cs.thermo_params)
end

"""
    turbulent_fluxes!(cs::CoupledSimulation)
    turbulent_fluxes!(fields, model_sims, thermo_params)

Compute turbulent fluxes and associated quantities. Store the results in `fields` as
area-weighted sums.

This function uses `SurfaceFluxes.jl` under the hood.

Args:
- `csf`: [Field of NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `thermo_params`: [TD.Parameters.ThermodynamicsParameters] the thermodynamic parameters.

TODO:
- generalize interface for regridding and take land state out of atmos's integrator.p
- add flux accumulation
- add flux bounds

(NB: Radiation surface fluxes are calculated by the atmosphere.)
"""
function turbulent_fluxes!(csf, model_sims, thermo_params)
    boundary_space = axes(csf)
    atmos_sim = model_sims.atmos_sim
    FT = CC.Spaces.undertype(boundary_space)

    # Reset the coupler fields will compute. We need to do this because we will compute
    # area-weighted averages
    for p in (
        :F_turb_ρτxz,
        :F_turb_ρτyz,
        :F_lh,
        :F_sh,
        :F_turb_moisture,
        :L_MO,
        :ustar,
        :buoyancy_flux,
    )
        fill!(getproperty(csf, p), 0)
    end

    # Compute the surface fluxes for each surface model and add them to `csf`
    for sim in model_sims
        # If the simulation is an implicit flux simulation, the fluxes are computed in the
        # component model's `step!` function, so we don't need to compute them here.
        sim isa Interfacer.ImplicitFluxSimulation ||
            compute_surface_fluxes!(csf, sim, atmos_sim, thermo_params)
    end
    return nothing
end

"""
    get_surface_fluxes(inputs, surface_params::SF.Parameters.SurfaceFluxesParameters)

Uses SurfaceFluxes.jl to calculate turbulent surface fluxes. It should be atmos model agnostic, and columnwise.
Fluxes are computed over the entire surface, even where the relevant surface model is not present.

When available, it also computes ancillary quantities, such as the Monin-Obukov lengthscale.
"""
function get_surface_fluxes(
    surface_fluxes_params::SF.Parameters.SurfaceFluxesParameters,
    u_int,
    thermo_state_atmos,
    h_int,
    u_sfc,
    thermo_state_sfc,
    h_sfc,
    d,
    config,
)
    # Get inputs to compute surface fluxes
    thermo_params = SFP.thermodynamics_params(surface_fluxes_params)
    T_int = TD.air_temperature(thermo_params, thermo_state_atmos)
    q_tot_int = TD.total_specific_humidity(thermo_params, thermo_state_atmos)
    q_liq_int = TD.liquid_specific_humidity(thermo_params, thermo_state_atmos)
    q_ice_int = TD.ice_specific_humidity(thermo_params, thermo_state_atmos)
    ρ_int = TD.air_density(thermo_params, thermo_state_atmos)
    T_sfc = TD.air_temperature(thermo_params, thermo_state_sfc)
    q_vap_sfc = TD.vapor_specific_humidity(thermo_params, thermo_state_sfc)
    Φ_sfc = SFP.grav(surface_fluxes_params) * h_sfc
    Δz = h_int - h_sfc

    # Calculate surface fluxes
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
        nothing, # roughness_inputs
        config,
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
    )
end

"""
    get_surface_params(atmos_sim::Interfacer.AtmosModelSimulation)

Returns the surface parameters of type `SF.Parameters.SurfaceFluxesParameters`.

TODO: in the future this may not need to depend on the atmos sim, but
here retaining the dependency until we know how EDMF boundary conditions will
be handled (for consistency of parameters).
"""
function get_surface_params(atmos_sim::Interfacer.AtmosModelSimulation)
    return error(
        "get_surface_params is required to be dispatched on $(nameof(atmos_sim)), but no method defined",
    )
end

"""
    update_turbulent_fluxes!(sim::Interfacer.ComponentModelSimulation, fields::NamedTuple)

Updates the fluxes in the simulation `sim` with the fluxes in `fields`.

For surface models, this should be the fluxes computed between the surface model and the atmosphere.
For atmosphere models, this should be the area-weighted sum of fluxes across all surface models.
"""
function update_turbulent_fluxes!(sim::Interfacer.ComponentModelSimulation, fields)
    return error(
        "update_turbulent_fluxes! is required to be dispatched on $(nameof(sim)), but no method defined",
    )
end

update_turbulent_fluxes!(sim::Interfacer.AbstractSurfaceStub, fields::NamedTuple) = nothing

"""
    compute_surface_humidity!(q_sfc, T_sfc, ρ_sfc, thermo_params)

Computes the surface specific humidity based on the atmospheric state and surface temperature.
The phase of the surface is determined by the surface temperature, and the saturation
specific humidity is computed accordingly.

All fields should be on the exchange grid.

# Arguments
- `q_sfc`: [CC.Fields.Field] output field for surface specific humidity.
- `T_sfc`: [CC.Fields.Field] surface temperature.
- `ρ_sfc`: [CC.Fields.Field] surface air density.
- `thermo_params`: [TD.Parameters.ThermodynamicsParameters] the thermodynamic parameters.
"""
# TODO: use TD.q_vap_saturation after updating to the latest Thermodynamics
function compute_surface_humidity!(q_sfc, T_sfc, ρ_sfc, thermo_params)
    T_freeze = TDP.T_freeze(thermo_params)
    @. q_sfc = ifelse(
        T_sfc .> T_freeze,
        TD.q_vap_saturation_generic(thermo_params, T_sfc, ρ_sfc, TD.Liquid()),
        TD.q_vap_saturation_generic(thermo_params, T_sfc, ρ_sfc, TD.Ice()),
    )
    return nothing
end

"""
    compute_surface_fluxes!(csf, sim, atmos_sim, thermo_params)

This function computes surface fluxes between the input component model
simulation and the atmosphere.

Update the input coupler surface fields `csf` in-place with the computed fluxes
for this model. These are then summed using area-weighting across all surface
models to get the total fluxes.

Since the fluxes are computed between the input model and the atmosphere, this
function does nothing if called on an atmosphere model simulation.

The function for ImplicitFluxSimulation is a placeholder that does nothing. Currently,
the only ImplicitFluxSimulation is ClimaLandSimulation, for which compute_surface_fluxes!
is defined in the component model. We can extend this function for other ImplicitFluxSimulation
in the future.

# Arguments
- `csf`: [CC.Fields.Field] containing a NamedTuple of turbulent flux fields: `F_turb_ρτxz`, `F_turb_ρτyz`, `F_lh`, `F_sh`, `F_turb_moisture`.
- `sim`: [Interfacer.ComponentModelSimulation] the surface simulation to compute fluxes for.
- `atmos_sim`: [Interfacer.AtmosModelSimulation] the atmosphere simulation to compute fluxes with.
- `thermo_params`: [TD.Parameters.ThermodynamicsParameters] the thermodynamic parameters.

The roughness model is obtained from the simulation via `get_field(sim, Val(:roughness_model))`.
Ocean simulations return `:coare3`, while land and ice simulations return `:constant` (the default).
"""
function compute_surface_fluxes!(
    csf,
    sim::Interfacer.AtmosModelSimulation,
    atmos_sim::Interfacer.AtmosModelSimulation,
    thermo_params,
)
    # do nothing for atmos model
    return nothing
end

function compute_surface_fluxes!(
    csf,
    sim::Interfacer.ImplicitFluxSimulation,
    atmos_sim::Interfacer.AtmosModelSimulation,
    thermo_params,
)
    # do nothing for implicit flux surface model
    return nothing
end

function compute_surface_fluxes!(
    csf,
    sim::Interfacer.SurfaceModelSimulation,
    atmos_sim::Interfacer.AtmosModelSimulation,
    thermo_params,
)
    boundary_space = axes(csf)
    FT = CC.Spaces.undertype(boundary_space)
    surface_fluxes_params = FluxCalculator.get_surface_params(atmos_sim)

    # Atmosphere fields are stored in coupler fields so we only regrid them once per timestep
    # `_int` refers to atmos state of center level 1
    uv_int = StaticArrays.SVector.(csf.u_int, csf.v_int)

    # construct the atmospheric thermo state
    thermo_state_atmos =
        TD.PhaseNonEquil_ρTq.(
            thermo_params,
            csf.ρ_atmos,
            csf.T_atmos,
            TD.PhasePartition.(csf.q_atmos),
        )

    # compute surface humidity from the surface temperature, surface density, and phase
    Interfacer.get_field!(csf.scalar_temp1, sim, Val(:surface_temperature))
    T_sfc = csf.scalar_temp1

    # TODO: This is not accurate - we shouldn't assume condensate is 0.
    ρ_sfc =
        SF.surface_density.(
            surface_fluxes_params,
            csf.T_atmos,
            csf.ρ_atmos,
            T_sfc,
            csf.height_int .- csf.height_sfc,
            csf.q_atmos,
            0, # q_liq
            0, # q_ice
        )

    compute_surface_humidity!(csf.scalar_temp2, T_sfc, ρ_sfc, thermo_params)
    q_sfc = csf.scalar_temp2

    # construct the surface thermo state
    # after this we can reuse `scalar_temp1` and `scalar_temp2` again
    thermo_state_sfc =
        TD.PhaseNonEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, TD.PhasePartition.(q_sfc))

    # get area fraction (min = 0, max = 1)
    Interfacer.get_field!(csf.scalar_temp1, sim, Val(:area_fraction))
    area_fraction = csf.scalar_temp1

    # Set some scalars that we hardcode for now
    gustiness = ones(boundary_space)

    # Get roughness model from the simulation (ocean simulations return :coare3, others return :constant)
    roughness_model = Interfacer.get_field(sim, Val(:roughness_model))
    # Get roughness lengths for constant roughness model
    if roughness_model == :constant
        Interfacer.get_field!(csf.scalar_temp2, sim, Val(:roughness_momentum))
        z0m = csf.scalar_temp2
        Interfacer.get_field!(csf.scalar_temp3, sim, Val(:roughness_buoyancy))
        z0b = csf.scalar_temp3
    end

    # Set SurfaceFluxConfig containing models for roughness and gustiness
    # Reuse cached roughness params when available to avoid allocations
    roughness_params = if roughness_model == :coare3
        # All COARE3-using simulations cache coare3_roughness_params, so no allocation needed
        Interfacer.get_field(sim, Val(:coare3_roughness_params))
    elseif roughness_model == :constant
        SF.ConstantRoughnessParams.(z0m, z0b)
    else
        error("Unknown roughness_model: $roughness_model. Must be :coare3 or :constant")
    end

    config = SF.SurfaceFluxConfig.(roughness_params, SF.ConstantGustinessSpec.(gustiness))

    # calculate the surface fluxes
    fluxes =
        FluxCalculator.get_surface_fluxes.(
            surface_fluxes_params,
            uv_int,
            thermo_state_atmos,
            csf.height_int,
            StaticArrays.SVector.(0, 0), # uv_sfc
            thermo_state_sfc,
            csf.height_sfc,
            0, # d
            config,
        )
    (; F_turb_ρτxz, F_turb_ρτyz, F_sh, F_lh, F_turb_moisture, L_MO, ustar, buoyancy_flux) =
        fluxes

    # Zero out fluxes where the area fraction is zero
    # Multiplying by `area_fraction` is not sufficient because the fluxes may
    # be NaN where the area fraction is zero.
    @. F_turb_ρτxz = ifelse(area_fraction ≈ 0, zero(F_turb_ρτxz), F_turb_ρτxz)
    @. F_turb_ρτyz = ifelse(area_fraction ≈ 0, zero(F_turb_ρτyz), F_turb_ρτyz)
    @. F_sh = ifelse(area_fraction ≈ 0, zero(F_sh), F_sh)
    @. F_lh = ifelse(area_fraction ≈ 0, zero(F_lh), F_lh)
    @. F_turb_moisture = ifelse(area_fraction ≈ 0, zero(F_turb_moisture), F_turb_moisture)
    @. L_MO = ifelse(area_fraction ≈ 0, zero(L_MO), L_MO)
    @. ustar = ifelse(area_fraction ≈ 0, zero(ustar), ustar)
    @. buoyancy_flux = ifelse(area_fraction ≈ 0, zero(buoyancy_flux), buoyancy_flux)

    # update the fluxes, which are now area-weighted, of this surface model
    fields = (; F_turb_ρτxz, F_turb_ρτyz, F_lh, F_sh, F_turb_moisture)
    FluxCalculator.update_turbulent_fluxes!(sim, fields)

    # update fluxes in the coupler fields
    # add the flux contributing from this surface to the coupler field
    # note that the fluxes are area-weighted here when we provide them to the atmosphere
    @. csf.F_turb_ρτxz += F_turb_ρτxz * area_fraction
    @. csf.F_turb_ρτyz += F_turb_ρτyz * area_fraction
    @. csf.F_lh += F_lh * area_fraction
    @. csf.F_sh += F_sh * area_fraction
    @. csf.F_turb_moisture += F_turb_moisture * area_fraction

    # NOTE: This is still an area weighted contribution, which maybe doesn't make
    # too much sense for these quantities...

    # L_MO can be Inf. We don't want to multiply Inf * 0, so we can handle this
    # separately.
    @. csf.L_MO += ifelse(isinf(L_MO), L_MO, L_MO * area_fraction)
    @. csf.ustar += ustar * area_fraction
    @. csf.buoyancy_flux += buoyancy_flux * area_fraction
    return nothing
end

"""
    ocean_seaice_fluxes!(cs::CoupledSimulation)
    ocean_seaice_fluxes!(ocean_sim, ice_sim)

Compute the fluxes between the ocean and sea ice simulations.
This function does nothing by default - it should be extended
for any ocean and sea ice models that support flux calculations.
"""
function ocean_seaice_fluxes!(cs::Interfacer.CoupledSimulation)
    haskey(cs.model_sims, :ocean_sim) &&
        haskey(cs.model_sims, :ice_sim) &&
        ocean_seaice_fluxes!(cs.model_sims.ocean_sim, cs.model_sims.ice_sim)
    return nothing
end
function ocean_seaice_fluxes!(
    ocean_sim::Union{Interfacer.OceanModelSimulation, Interfacer.AbstractSurfaceStub},
    ice_sim::Union{Interfacer.SeaIceModelSimulation, Interfacer.AbstractSurfaceStub},
)
    return nothing
end

end # module
