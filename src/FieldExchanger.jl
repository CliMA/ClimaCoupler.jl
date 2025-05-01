"""
    FieldExchanger

This modules contains general functions for the exchange of fields between the
atmospheric and surface component models.
"""
module FieldExchanger

import ..Interfacer, ..FluxCalculator, ..Utilities
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

export import_atmos_fields!,
    update_surface_fractions!,
    import_combined_surface_fields!,
    update_sim!,
    update_model_sims!,
    reinit_model_sims!,
    step_model_sims!

"""
    update_surface_fractions!(cs::Interfacer.CoupledSimulation)

Updates dynamically changing area fractions.
Maintains the invariant that the sum of area fractions is 1 at all points.

# Arguments
- `cs`: [Interfacer.CoupledSimulation] containing area fraction information.
"""
function update_surface_fractions!(cs::Interfacer.CoupledSimulation)
    FT = Interfacer.float_type(cs)

    # land fraction is static
    land_fraction = Interfacer.get_field(cs.model_sims.land_sim, Val(:area_fraction))

    # ice and ocean fractions are dynamic
    ice_fraction_before = Interfacer.get_field(cs.model_sims.ice_sim, Val(:area_fraction))
    # max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
    Interfacer.update_field!(
        cs.model_sims.ice_sim,
        Val(:area_fraction),
        max.(min.(ice_fraction_before, FT(1) .- land_fraction), FT(0)),
    )
    ice_fraction = Interfacer.get_field(cs.model_sims.ice_sim, Val(:area_fraction))

    Interfacer.update_field!(
        cs.model_sims.ocean_sim,
        Val(:area_fraction),
        max.(FT(1) .- ice_fraction .- land_fraction, FT(0)),
    )
    ocean_fraction = Interfacer.get_field(cs.model_sims.ocean_sim, Val(:area_fraction))

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
    dummmy_remap!(csf.ρ_sfc, FluxCalculator.calculate_surface_air_density(atmos_sim, csf.T_sfc)) # TODO: generalize to use individual T_sfc, (#445)

    # radiative fluxes
    dummmy_remap!(csf.F_radiative, Interfacer.get_field(atmos_sim, Val(:radiative_energy_flux_sfc)))

    # precipitation
    dummmy_remap!(csf.P_liq, Interfacer.get_field(atmos_sim, Val(:liquid_precipitation)))
    dummmy_remap!(csf.P_snow, Interfacer.get_field(atmos_sim, Val(:snow_precipitation)))
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
    combined_field = csf.temp1

    # surface fields
    combine_surfaces!(combined_field, model_sims, Val(:surface_temperature))
    dummmy_remap!(csf.T_sfc, combined_field)

    combine_surfaces!(combined_field, model_sims, Val(:surface_direct_albedo))
    dummmy_remap!(csf.surface_direct_albedo, combined_field)

    combine_surfaces!(combined_field, model_sims, Val(:surface_diffuse_albedo))
    dummmy_remap!(csf.surface_diffuse_albedo, combined_field)
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
end

"""
    update_model_sims!(model_sims, csf)

Iterates `update_sim!` over all component model simulations saved in `cs.model_sims`.

# Arguments
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_model_sims!(model_sims, csf)
    for sim in model_sims
        if sim isa Interfacer.SurfaceModelSimulation
            update_sim!(sim, csf, Interfacer.get_field(sim, Val(:area_fraction)))
        else
            update_sim!(sim, csf)
        end
    end
end

"""
    reinit_model_sims!(model_sims)

Iterates `reinit!` over all component model simulations saved in `cs.model_sims`.

# Arguments
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.

"""
function reinit_model_sims!(model_sims)
    for sim in model_sims
        Interfacer.reinit!(sim)
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
    dummmy_remap!(target, source)

Simple stand-in function for remapping.
For AMIP we don't need regridding of surface model CC.Fields.
When we do, we re-introduce the ClimaCoreTempestRemap remapping functions.

# Arguments
- `target`: [CC.Fields.Field] destination of remapping.
- `source`: [CC.Fields.Field] source of remapping.
"""
function dummmy_remap!(target, source)
    parent(target) .= parent(source)
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
    combined_field .= 0
    for sim in sims
        if sim isa Interfacer.SurfaceModelSimulation
            # Zero out the contribution from this surface if the area fraction is zero
            # Note that multiplying by `area_fraction` is not sufficient in the case of NaNs
            area_fraction = Interfacer.get_field(sim, Val(:area_fraction))
            combined_field .+=
                area_fraction .*
                ifelse.(area_fraction .≈ 0, zero(combined_field), Interfacer.get_field(sim, field_name))
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


end # module
