"""
    FieldExchanger

This modules contains general functions for the exchange of fields between the
atmospheric and surface component models.
"""
module FieldExchanger

import ClimaCore as CC

import ..Interfacer, ..FluxCalculator, ..Utilities

export update_sim!,
    update_model_sims!,
    step_model_sims!,
    exchange!,
    set_caches!,
    update_surface_fractions!,
    resolve_area_fractions!

"""
    update_surface_fractions!(cs::Interfacer.CoupledSimulation)

Updates dynamically changing area fractions.
Maintains the invariant that the sum of area fractions is 1 at all points.
Area fractions are expected to be defined on the boundary space of the coupled simulation,
since they are used by the coupler.

If a surface model is not present, the area fraction is set to 0.

# Arguments
- `cs`: [Interfacer.CoupledSimulation] containing area fraction information.
"""
function update_surface_fractions!(cs::Interfacer.CoupledSimulation)
    FT = CC.Spaces.undertype(Interfacer.boundary_space(cs))

    # land fraction is static
    if haskey(cs.model_sims, :land_sim)
        land_fraction = Interfacer.get_field(cs.model_sims.land_sim, Val(:area_fraction))
    else
        cs.fields.scalar_temp1 .= 0
        land_fraction = cs.fields.scalar_temp1
    end

    # ice and ocean fractions are dynamic
    if haskey(cs.model_sims, :ice_sim)
        ice_sim = cs.model_sims.ice_sim
        Interfacer.get_field!(cs.fields.scalar_temp1, ice_sim, Val(:ice_concentration))
        ice_concentration = cs.fields.scalar_temp1

        # max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
        Interfacer.update_field!(
            ice_sim,
            Val(:area_fraction),
            max.(min.(ice_concentration, FT(1) .- land_fraction), FT(0)),
        )
        ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction))
    else
        cs.fields.scalar_temp1 .= 0
        ice_fraction = cs.fields.scalar_temp1
    end

    if haskey(cs.model_sims, :ocean_sim)
        ocean_sim = cs.model_sims.ocean_sim
        Interfacer.update_field!(
            ocean_sim,
            Val(:area_fraction),
            max.(FT(1) .- ice_fraction .- land_fraction, FT(0)),
        )
        ocean_fraction = Interfacer.get_field(ocean_sim, Val(:area_fraction))

        # ensure that ocean and ice fractions are consistent
        if haskey(cs.model_sims, :ice_sim)
            resolve_area_fractions!(ocean_sim, cs.model_sims.ice_sim, land_fraction)
        end
    else
        cs.fields.scalar_temp1 .= 0
        ocean_fraction = cs.fields.scalar_temp1
    end

    # check that the sum of area fractions is 1
    @assert minimum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ FT(1)
    @assert maximum(ice_fraction .+ land_fraction .+ ocean_fraction) ≈ FT(1)
end

"""
    resolve_area_fractions!(ocean_sim, ice_sim, land_fraction)

Ensure that the ocean and ice fractions are consistent with each other.
For most ocean and ice models, this does nothing since the ocean fraction is
defined as `1 - ice_fraction - land_fraction`. However, some models may have
additional constraints on the ice and ocean fractions that need to be enforced.
This function can be extended for such models.
"""
function resolve_area_fractions!(ocean_sim, ice_sim, land_fraction)
    return nothing
end

"""
    import_atmos_fields!(csf, model_sims)

Update the coupler with quantities from the  atmosphere model. By default, this
updates the coupler fields for quantities required for turbulent flux calculations,
radiative fluxes, and precipitation.
This function should be extended for any model that requires additional fields
from the atmosphere.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
"""
function import_atmos_fields!(csf, model_sims)
    # get atmosphere properties used for flux calculations
    Interfacer.get_field!(csf.T_atmos, model_sims.atmos_sim, Val(:air_temperature))
    Interfacer.get_field!(csf.q_atmos, model_sims.atmos_sim, Val(:specific_humidity))
    Interfacer.get_field!(csf.ρ_atmos, model_sims.atmos_sim, Val(:air_density))
    Interfacer.get_field!(csf.u_int, model_sims.atmos_sim, Val(:u_int))
    Interfacer.get_field!(csf.v_int, model_sims.atmos_sim, Val(:v_int))

    # radiative fluxes
    Interfacer.get_field!(csf.SW_d, model_sims.atmos_sim, Val(:SW_d))
    Interfacer.get_field!(csf.LW_d, model_sims.atmos_sim, Val(:LW_d))

    # precipitation
    Interfacer.get_field!(csf.P_liq, model_sims.atmos_sim, Val(:liquid_precipitation))
    Interfacer.get_field!(csf.P_snow, model_sims.atmos_sim, Val(:snow_precipitation))

    for sim in model_sims
        import_atmos_fields!(csf, sim, model_sims.atmos_sim)
    end
end

"""
    import_atmos_fields!(csf, ::Interfacer.ComponentModelSimulation, atmos_sim)

Updates the coupler simulation fields with atmospheric fluxes from the atmosphere simulation.
This function should be extended for any surface model that requires additional fields
from the atmosphere. Any fields added in a method of this function should also be added
in the corresponding method of `Interfacer.add_coupler_fields!`. The combination of these
two functions defines any extra atmosphere fields provided to the surface.
"""
import_atmos_fields!(csf, ::Interfacer.ComponentModelSimulation, atmos_sim) = nothing

"""
    import_combined_surface_fields!(csf, model_sims)

Updates the coupler with the surface properties. The `Interfacer.get_field`
functions for (`:emissivity`, `:surface_temperature`, `:surface_direct_albedo`,
`:surface_diffuse_albedo`) need to be specified for each surface model.

Note: Not all surface fields are imported here. Some quantities are retrieved
from each surface model when surface fluxes are computed, in `compute_surface_fluxes!`.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
"""
function import_combined_surface_fields!(csf, model_sims)
    combine_surfaces!(csf, model_sims, Val(:emissivity))
    combine_surfaces!(csf, model_sims, Val(:surface_temperature))
    combine_surfaces!(csf, model_sims, Val(:surface_direct_albedo))
    combine_surfaces!(csf, model_sims, Val(:surface_diffuse_albedo))
    return nothing
end

"""
    import_static_fields!(csf, model_sims)

Import static fields into the coupler fields.
This is used to import fields that are not updated during a simulation,
so it is only called at initialization.

Fields imported here are:
- the bottom cell center and face heights of the atmosphere

Any fields imported into the coupler fields here that need to be sent to
component models should be updated in the `set_cache!` function
of the receiving component model.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
"""
function import_static_fields!(csf, model_sims)
    Interfacer.get_field!(csf.height_int, model_sims.atmos_sim, Val(:height_int))
    Interfacer.get_field!(csf.height_sfc, model_sims.atmos_sim, Val(:height_sfc))
    csf.height_delta .= Interfacer.get_atmos_height_delta(csf.height_int, csf.height_sfc)

    return nothing
end

"""
    update_sim!(atmos_sim::Interfacer.AtmosModelSimulation, csf)

Updates the atmosphere's fields for surface direct and diffuse albedos, emissivity,and temperature.

# Arguments
- `atmos_sim`: [Interfacer.AtmosModelSimulation] containing an atmospheric model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(atmos_sim::Interfacer.AtmosModelSimulation, csf)
    Interfacer.update_field!(
        atmos_sim,
        Val(:surface_direct_albedo),
        csf.surface_direct_albedo,
    )
    Interfacer.update_field!(
        atmos_sim,
        Val(:surface_diffuse_albedo),
        csf.surface_diffuse_albedo,
    )
    Interfacer.update_field!(atmos_sim, Val(:emissivity), csf.emissivity)
    Interfacer.update_field!(atmos_sim, Val(:surface_temperature), csf)
    return nothing
end

"""
    update_sim!(sim::SurfaceModelSimulation, csf)

Updates the surface component model cache with the current coupler fields
*besides turbulent fluxes*, which are updated in `update_turbulent_fluxes`.

# Arguments
- `sim`: [Interfacer.SurfaceModelSimulation] containing a surface model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(sim::Interfacer.SurfaceModelSimulation, csf)
    # radiative fluxes
    Interfacer.update_field!(sim, Val(:SW_d), csf.SW_d)
    Interfacer.update_field!(sim, Val(:LW_d), csf.LW_d)

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
    for sim in model_sims
        update_sim!(sim, csf)
    end
end

"""
    step_model_sims!(model_sims, t, coupler_fields, thermo_params)
    step_model_sims!(cs::CoupledSimulation)

Iterates `step!` over all component model simulations saved in `cs.model_sims`.

# Arguments
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `t`: [AbstractFloat or ITime] denoting the simulation time.
"""
function step_model_sims!(model_sims, t, coupler_fields, thermo_params)
    # Step all surface models (the ordering doesn't matter here)
    for sim in model_sims
        sim isa Interfacer.SurfaceModelSimulation && Interfacer.step!(sim, t)

        # For an implicit flux simulation, `compute_surface_fluxes!` reads in the precomputed
        # fluxes, and puts them into the coupler fields.
        sim isa Interfacer.ImplicitFluxSimulation && FluxCalculator.compute_surface_fluxes!(
            coupler_fields,
            sim,
            model_sims.atmos_sim,
            thermo_params,
        )
    end

    # Update the atmosphere with the fluxes across all surface models
    # The surface models have already been updated with the fluxes in `compute_surface_fluxes!`,
    # or internally within the step in the case of the integrated land model.
    FluxCalculator.update_turbulent_fluxes!(model_sims.atmos_sim, coupler_fields)

    # Step the atmosphere model
    Interfacer.step!(model_sims.atmos_sim, t)
    return nothing
end

function step_model_sims!(cs::Interfacer.CoupledSimulation)
    step_model_sims!(cs.model_sims, cs.t[], cs.fields, cs.thermo_params)
end

"""
    combine_surfaces!(csf, sims, field_name_val::Val{field_name}) where {field_name}

Sums the surface fields specified by `field_name_val`, weighted by the respective area fractions
of all surface simulations. The result is saved in the coupler field specified by `field_name_val`.

For surface temperature, upward longwave radiation is computed from the temperatures
of each surface, weighted by their area fractions, and then the combined temperature
is computed from the combined upward longwave radiation.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
    Note: For the surface temperature, all coupler fields are passed in a NamedTuple.
- `sims`: [NamedTuple] containing simulations.
- `field_name_val`: [Val] containing the name Symbol of the field to be extracted by the `Interfacer.get_field` functions.

# Example
- `combine_surfaces!(temp_field, cs.model_sims, Val(:emissivity))`
"""
function combine_surfaces!(csf, sims, field_name_val::Val{field_name}) where {field_name}
    # Extract the coupler field we are updating
    combined_field = getproperty(csf, field_name)
    FT = eltype(combined_field)
    combined_field .= zero(FT)

    for sim in sims
        if sim isa Interfacer.SurfaceModelSimulation
            # Store the area fraction of this simulation in `scalar_temp` and rename for clarity
            Interfacer.get_field!(csf.scalar_temp1, sim, Val(:area_fraction))
            area_fraction = csf.scalar_temp1

            # Remap the surface field onto a coupler temporary field to avoid allocation
            Interfacer.get_field!(csf.scalar_temp2, sim, field_name_val)
            surface_field = csf.scalar_temp2

            # Zero out the contribution from this surface if the area fraction is zero.
            # Note that multiplying by `area_fraction` is not sufficient in the case of NaNs
            combined_field .+=
                area_fraction .* ifelse.(area_fraction .≈ 0, zero(FT), surface_field)
        end
    end
    return nothing
end
function combine_surfaces!(csf, sims, ::Val{:surface_temperature})
    # extract the coupler fields we need to get the surface temperature
    T_sfc = csf.T_sfc
    emissivity_sfc = csf.emissivity

    FT = eltype(T_sfc)
    T_sfc .= zero(FT)
    for sim in sims
        if sim isa Interfacer.SurfaceModelSimulation
            # Store the area fraction and emissivity of this simulation in temp fields
            Interfacer.get_field!(csf.scalar_temp1, sim, Val(:area_fraction))
            area_fraction = csf.scalar_temp1
            Interfacer.get_field!(csf.scalar_temp2, sim, Val(:emissivity))
            emissivity_sim = csf.scalar_temp2

            # Remap the surface field onto a coupler temporary field to avoid allocation
            Interfacer.get_field!(csf.scalar_temp3, sim, Val(:surface_temperature))
            T_sfc_sim = csf.scalar_temp3

            # Zero out the contribution from this surface if the area fraction is zero.
            # Note that multiplying by `area_fraction` is not sufficient in the case of NaNs
            # Compute upward longwave radiation from surface temperature for this simulation
            T_sfc .+=
                area_fraction .*
                ifelse.(area_fraction .≈ 0, zero(FT), emissivity_sim .* T_sfc_sim .^ FT(4))
        end
    end
    # Convert the combined upward longwave radiation into a surface temperature
    @. T_sfc = (T_sfc / emissivity_sfc)^FT(1 / 4)
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

"""
    set_caches!(cs::Interfacer.CoupledSimulation)

Perform any initialization of the component model caches that cannot be
done before the initial exchange. This is useful in handling cache interdependencies
between component models.

For example, the radiation callback in the atmosphere model needs to be
initialized with the surface temperatures, which are only available after the
initial exchange. The integrated land, in turn, requires its drivers in the
cache to be filled with the initial radiation fluxes, so that it can propagate
these to the rest of its cache (e.g. in canopy radative transfer).

This function can also be used to set exchanged fields that are static over the
simulation, since it is only called at initialization.
"""
function set_caches!(cs::Interfacer.CoupledSimulation)
    Interfacer.set_cache!(cs.model_sims.atmos_sim, cs.fields)
    exchange!(cs)
    for sim in cs.model_sims
        sim isa Interfacer.SurfaceModelSimulation && Interfacer.set_cache!(sim, cs.fields)
    end
    return nothing
end

end # module
