"""
    FieldExchanger

This modules contains general functions for the exchange of fields between the
atmospheric and surface component models.
"""
module FieldExchanger

import ClimaCore as CC

import ..Interfacer, ..FluxCalculator, ..Utilities
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

export update_sim!, update_model_sims!, step_model_sims!, exchange!

"""
    update_surface_fractions!(cs::Interfacer.CoupledSimulation)

Updates dynamically changing area fractions.
Maintains the invariant that the sum of area fractions is 1 at all points.

If a surface model is not present, the area fraction is set to 0.

# Arguments
- `cs`: [Interfacer.CoupledSimulation] containing area fraction information.
"""
function update_surface_fractions!(cs::Interfacer.CoupledSimulation)
    boundary_space = cs.boundary_space
    FT = CC.Spaces.undertype(boundary_space)

    # land fraction is static
    if haskey(cs.model_sims, :land_sim)
        land_fraction = Interfacer.get_field(cs.model_sims.land_sim, Val(:area_fraction), boundary_space)
    else
        cs.fields.temp1 .= 0
        land_fraction = cs.fields.temp1
    end

    # ice and ocean fractions are dynamic
    if haskey(cs.model_sims, :ice_sim)
        ice_sim = cs.model_sims.ice_sim
        ice_fraction_before = Interfacer.get_field(ice_sim, Val(:area_fraction), boundary_space)
        # max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
        Interfacer.update_field!(
            ice_sim,
            Val(:area_fraction),
            max.(min.(ice_fraction_before, FT(1) .- land_fraction), FT(0)),
        )
        ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction), boundary_space)
    else
        cs.fields.temp1 .= 0
        ice_fraction = cs.fields.temp1
    end

    if haskey(cs.model_sims, :ocean_sim)
        ocean_sim = cs.model_sims.ocean_sim
        Interfacer.update_field!(ocean_sim, Val(:area_fraction), max.(FT(1) .- ice_fraction .- land_fraction, FT(0)))
        ocean_fraction = Interfacer.get_field(ocean_sim, Val(:area_fraction), boundary_space)
    else
        cs.fields.temp1 .= 0
        ocean_fraction = cs.fields.temp1
    end

    # check that the sum of area fractions is 1
    @assert minimum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ FT(1)
    @assert maximum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ FT(1)
end

"""
    import_atmos_fields!(csf, model_sims)

Update the coupler with the atmospheric fluxes. By default, this updates the coupler fields for
the surface air density, radiative fluxes, and precipitation. This function should be
extended for any model that requires additional fields from the atmosphere.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
"""
function import_atmos_fields!(csf, model_sims)
    for sim in model_sims
        import_atmos_fields!(csf, sim, model_sims.atmos_sim)
    end
end

"""
    import_atmos_fields!(csf, ::Interfacer.SurfaceModelSimulation, atmos_sim)

Updates the coupler simulation fields with atmospheric fluxes from the atmosphere simulation.
This is the default function to be used for most surface model simulations, as
    are computed by the coupler or atmosphere
and passed to the surfaces.
"""
function import_atmos_fields!(csf, ::Interfacer.SurfaceModelSimulation, atmos_sim)
    # surface density - needed for q_sat and requires atmos and sfc states, so it is calculated and saved in the coupler
    Interfacer.remap!(csf.ρ_sfc, FluxCalculator.calculate_surface_air_density(atmos_sim, csf.T_sfc)) # TODO: generalize to use individual T_sfc, (#445)

    # radiative fluxes
    Interfacer.get_field!(csf.F_radiative, atmos_sim, Val(:radiative_energy_flux_sfc))

    # precipitation
    Interfacer.get_field!(csf.P_liq, atmos_sim, Val(:liquid_precipitation))
    Interfacer.get_field!(csf.P_snow, atmos_sim, Val(:snow_precipitation))
    return nothing
end

import_atmos_fields!(csf, ::Interfacer.AtmosModelSimulation, atmos_sim) = nothing

"""
    import_combined_surface_fields!(csf, model_sims)

Updates the coupler with the surface properties. The `Interfacer.get_field`
functions for (`:surface_temperature`, `:surface_direct_albedo`,
`:surface_diffuse_albedo`) need to be specified for each surface model.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
"""
function import_combined_surface_fields!(csf, model_sims)
    combine_surfaces!(csf.T_sfc, model_sims, Val(:surface_temperature))
    combine_surfaces!(csf.surface_direct_albedo, model_sims, Val(:surface_direct_albedo))
    combine_surfaces!(csf.surface_diffuse_albedo, model_sims, Val(:surface_diffuse_albedo))
    return nothing
end

"""
    update_sim!(atmos_sim::Interfacer.AtmosModelSimulation, csf)

Updates the surface fields for temperature, roughness length, albedo, and specific humidity.

# Arguments
- `atmos_sim`: [Interfacer.AtmosModelSimulation] containing an atmospheric model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(atmos_sim::Interfacer.AtmosModelSimulation, csf)
    Interfacer.update_field!(atmos_sim, Val(:surface_direct_albedo), csf.surface_direct_albedo)
    Interfacer.update_field!(atmos_sim, Val(:surface_diffuse_albedo), csf.surface_diffuse_albedo)
    Interfacer.update_field!(atmos_sim, Val(:surface_temperature), csf)
    Interfacer.update_field!(atmos_sim, Val(:turbulent_fluxes), csf)
    return nothing
end

"""
    update_sim!(sim::SurfaceModelSimulation, csf, area_fraction)

Updates the surface component model cache with the current coupler fields besides turbulent fluxes.

# Arguments
- `sim`: [Interfacer.SurfaceModelSimulation] containing a surface model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(sim::Interfacer.SurfaceModelSimulation, csf, area_fraction)
    FT = eltype(area_fraction)

    # atmospheric surface density
    Interfacer.update_field!(sim, Val(:air_density), csf.ρ_sfc)

    # radiative fluxes
    Interfacer.update_field!(sim, Val(:radiative_energy_flux_sfc), FT.(area_fraction .* csf.F_radiative))

    # precipitation
    Interfacer.update_field!(sim, Val(:liquid_precipitation), csf.P_liq)
    Interfacer.update_field!(sim, Val(:snow_precipitation), csf.P_snow)
    return nothing
end

"""
    update_model_sims!(model_sims, csf)

Iterates `update_sim!` over all component model simulations saved in `cs.model_sims`.

# Arguments
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_model_sims!(model_sims, csf)
    boundary_space = axes(csf)
    for sim in model_sims
        if sim isa Interfacer.SurfaceModelSimulation
            update_sim!(sim, csf, Interfacer.get_field(sim, Val(:area_fraction), boundary_space))
        else
            update_sim!(sim, csf)
        end
    end
end

"""
    step_model_sims!(model_sims, t)

Iterates `step!` over all component model simulations saved in `cs.model_sims`.

# Arguments
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `t`: [AbstractFloat] denoting the simulation time.
"""
function step_model_sims!(model_sims, t)
    for sim in model_sims
        Interfacer.step!(sim, t)
    end
end

"""
    combine_surfaces!(combined_field::CC.Fields.Field, sims, field_name::Val)

Sums the fields, specified by `field_name`, weighted by the respective area fractions of all
surface simulations. THe result is saved in `combined_field`.

# Arguments
- `combined_field`: [CC.Fields.Field] output object containing weighted values.
- `sims`: [NamedTuple] containing simulations .
- `field_name`: [Val] containing the name Symbol of the field t be extracted by the `Interfacer.get_field` functions.

# Example
- `combine_surfaces!(temp_field, cs.model_sims, Val(:surface_temperature))`
"""
function combine_surfaces!(combined_field, sims, field_name)
    boundary_space = axes(combined_field)
    combined_field .= 0
    for sim in sims
        if sim isa Interfacer.SurfaceModelSimulation
            # Zero out the contribution from this surface if the area fraction is zero
            # Note that multiplying by `area_fraction` is not sufficient in the case of NaNs
            area_fraction = Interfacer.get_field(sim, Val(:area_fraction), boundary_space)
            combined_field .+=
                area_fraction .*
                ifelse.(area_fraction .≈ 0, zero(combined_field), Interfacer.get_field(sim, field_name, boundary_space))
        end
    end
end

"""
    compute_surface_humidity!(q_sfc, T_atmos, q_atmos, ρ_atmos, T_sfc, thermo_params)

Compute the surface specific humidity based on the atmospheric state and surface temperature.
The phase of the surface is determined by the surface temperature, and the saturation
specific humidity is computed accordingly.
All fields should be on the exchange grid.
# Arguments
- `q_sfc`: [CC.Fields.Field] output field for surface specific humidity.
- `T_atmos`: [CC.Fields.Field] atmospheric temperature.
- `q_atmos`: [CC.Fields.Field] atmospheric specific humidity.
- `ρ_atmos`: [CC.Fields.Field] atmospheric air density.
- `T_sfc`: [CC.Fields.Field] surface temperature.
- `thermo_params`: [TD.Parameters.ThermodynamicsParameters] the thermodynamic parameters.
"""
function compute_surface_humidity!(q_sfc, T_atmos, q_atmos, ρ_atmos, T_sfc, thermo_params)
    thermo_state_atmos = TD.PhaseEquil_ρTq.(thermo_params, ρ_atmos, T_atmos, q_atmos)
    ρ_sfc = FluxCalculator.extrapolate_ρ_to_sfc.(thermo_params, thermo_state_atmos, T_sfc)

    T_freeze = TDP.T_freeze(thermo_params)
    @. q_sfc = ifelse(
        T_sfc .> T_freeze,
        TD.q_vap_saturation_generic(thermo_params, T_sfc, ρ_sfc, TD.Liquid()),
        TD.q_vap_saturation_generic(thermo_params, T_sfc, ρ_sfc, TD.Ice()),
    )
    return nothing
end

"""
    exchange!(cs::Interfacer.CoupledSimulation)

Exchange fields between the surface and atmosphere models.
This is done in 2 steps:
1. Import the atmosphere fields and surface fields into the coupler.
2. Update the component model simulations with the coupler fields.

The order of these steps is important, as importing the surface fields requires
the atmosphere fields to be updated so that surface humidity can be computed.
"""
function exchange!(cs::Interfacer.CoupledSimulation)
    # Import the atmosphere fields and surface fields into the coupler
    import_atmos_fields!(cs.fields, cs.model_sims)
    import_combined_surface_fields!(cs.fields, cs.model_sims)

    # Update the component model simulations with the coupler fields
    update_model_sims!(cs.model_sims, cs.fields)
    return nothing
end


end # module
