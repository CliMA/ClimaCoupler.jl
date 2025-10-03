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

export update_sim!, update_model_sims!, step_model_sims!, exchange!, set_caches!

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
        cs.fields.temp1 .= 0
        land_fraction = cs.fields.temp1
    end

    # ice and ocean fractions are dynamic
    if haskey(cs.model_sims, :ice_sim)
        ice_sim = cs.model_sims.ice_sim
        ice_fraction_before = Interfacer.get_field(ice_sim, Val(:area_fraction))
        # max needed to avoid Float32 errors (see issue #271; Heisenbug on HPC)
        Interfacer.update_field!(
            ice_sim,
            Val(:area_fraction),
            max.(min.(ice_fraction_before, FT(1) .- land_fraction), FT(0)),
        )
        ice_fraction = Interfacer.get_field(ice_sim, Val(:area_fraction))
    else
        cs.fields.temp1 .= 0
        ice_fraction = cs.fields.temp1
    end

    if haskey(cs.model_sims, :ocean_sim)
        ocean_sim = cs.model_sims.ocean_sim
        Interfacer.update_field!(
            ocean_sim,
            Val(:area_fraction),
            max.(FT(1) .- ice_fraction .- land_fraction, FT(0)),
        )
        ocean_fraction = Interfacer.get_field(ocean_sim, Val(:area_fraction))
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
    # get atmosphere properties used for flux calculations
    Interfacer.get_field!(csf.T_atmos, atmos_sim, Val(:air_temperature))
    Interfacer.get_field!(csf.q_atmos, atmos_sim, Val(:specific_humidity))
    Interfacer.get_field!(csf.ρ_atmos, atmos_sim, Val(:air_density))

    # radiative fluxes
    Interfacer.get_field!(csf.F_radiative, atmos_sim, Val(:radiative_energy_flux_sfc))

    # precipitation
    Interfacer.get_field!(csf.P_liq, atmos_sim, Val(:liquid_precipitation))
    Interfacer.get_field!(csf.P_snow, atmos_sim, Val(:snow_precipitation))
    return nothing
end

import_atmos_fields!(csf, ::Interfacer.AtmosModelSimulation, atmos_sim) = nothing

"""
    import_combined_surface_fields!(csf, model_sims, thermo_params)

Updates the coupler with the surface properties. The `Interfacer.get_field`
functions for (`:surface_temperature`, `:surface_direct_albedo`,
`:surface_diffuse_albedo`) need to be specified for each surface model.

Note: The calculation of surface humidity done here uses atmospheric properties stored in
the coupled fields. For these values to be correct, this function should be called
after `import_atmos_fields!` in a timestep.

Note 2: Not all surface fields are imported here. Some quantities are retrieved
from each surface model when surface fluxes are computed, in `compute_surface_fluxes!`.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `thermo_params`: [TD.Parameters.ThermodynamicsParameters] the thermodynamic parameters.
"""
function import_combined_surface_fields!(csf, model_sims, thermo_params)
    combine_surfaces!(csf.emissivity, model_sims, Val(:emissivity), csf.temp1)
    combine_surfaces!(
        csf.surface_direct_albedo,
        model_sims,
        Val(:surface_direct_albedo),
        csf.temp1,
    )
    combine_surfaces!(
        csf.surface_diffuse_albedo,
        model_sims,
        Val(:surface_diffuse_albedo),
        csf.temp1,
    )
    # Temperature requires emissivity, so we provide all the coupler fields
    combine_surfaces!(csf, model_sims, Val(:surface_temperature))

    # q_sfc is computed from the atmosphere state and surface temperature, so it's handled differently
    # This is computed on the exchange grid, so there's no need to remap
    compute_surface_humidity!(
        csf.q_sfc,
        csf.T_atmos,
        csf.q_atmos,
        csf.ρ_atmos,
        csf.T_sfc,
        thermo_params,
    )
    return nothing
end

"""
    compute_surface_humidity!(q_sfc, T_atmos, q_atmos, ρ_atmos, T_sfc, thermo_params)

Computes the surface specific humidity based on the atmospheric state and surface temperature.
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
    update_sim!(atmos_sim::Interfacer.AtmosModelSimulation, csf)

Updates the atmosphere's fields for surface direct and diffuse albedos, emissivity, and temperature,
as well as the turbulent fluxes.

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

    # radiative fluxes
    Interfacer.update_field!(
        sim,
        Val(:radiative_energy_flux_sfc),
        FT.(area_fraction .* csf.F_radiative),
    )

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
            update_sim!(sim, csf, Interfacer.get_field(sim, Val(:area_fraction)))
        else
            update_sim!(sim, csf)
        end
    end
end

"""
    step_model_sims!(model_sims, t)
    step_model_sims!(cs::CoupledSimulation)

Iterates `step!` over all component model simulations saved in `cs.model_sims`.

# Arguments
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `t`: [AbstractFloat or ITime] denoting the simulation time.
"""
function step_model_sims!(model_sims, t)
    for sim in model_sims
        Interfacer.step!(sim, t)
    end
    return nothing
end

function step_model_sims!(cs::Interfacer.CoupledSimulation)
    step_model_sims!(cs.model_sims, cs.t[])
end

"""
    combine_surfaces!(combined_field::CC.Fields.Field, sims, field_name::Val, temp1)

Sums the fields, specified by `field_name`, weighted by the respective area fractions of all
surface simulations. THe result is saved in `combined_field`.

For surface temperature, upward longwave radiation is computed from the temperatures
of each surface, weighted by their area fractions, and then the combined temperature
is computed from the combined upward longwave radiation.

# Arguments
- `combined_field`: [CC.Fields.Field] output object containing weighted values.
    Note: For the surface temperature, all coupler fields are passed in a NamedTuple.
- `sims`: [NamedTuple] containing simulations .
- `field_name`: [Val] containing the name Symbol of the field t be extracted by the `Interfacer.get_field` functions.
- `temp1`: [CC.Fields.Field] temporary field for intermediate calculations.
    Omitted for surface temperature method.

# Example
- `combine_surfaces!(temp_field, cs.model_sims, Val(:surface_temperature))`
"""
function combine_surfaces!(combined_field, sims, field_name, temp1)
    boundary_space = axes(combined_field)
    combined_field .= 0
    for sim in sims
        if sim isa Interfacer.SurfaceModelSimulation
            # Store the area fraction of this simulation in `temp1`
            Interfacer.get_field!(temp1, sim, Val(:area_fraction))
            # Zero out the contribution from this surface if the area fraction is zero.
            # Note that multiplying by `area_fraction` is not sufficient in the case of NaNs
            combined_field .+=
                temp1 .*
                ifelse.(
                    temp1 .≈ 0,
                    zero(combined_field),
                    Interfacer.get_field(sim, field_name, boundary_space),
                )
        end
    end
    return nothing
end
function combine_surfaces!(csf, sims, field_name::Val{:surface_temperature})
    # extract the coupler fields we need to get the surface temperature
    T_sfc = csf.T_sfc
    emissivity_sfc = csf.emissivity

    boundary_space = axes(T_sfc)
    FT = CC.Spaces.undertype(boundary_space)

    T_sfc .= FT(0)
    for sim in sims
        if sim isa Interfacer.SurfaceModelSimulation
            # Store the area fraction and emissivity of this simulation in temp fields
            Interfacer.get_field!(csf.temp1, sim, Val(:area_fraction))
            area_fraction = csf.temp1
            Interfacer.get_field!(csf.temp2, sim, Val(:emissivity))
            emissivity_sim = csf.temp2

            # Zero out the contribution from this surface if the area fraction is zero.
            # Note that multiplying by `area_fraction` is not sufficient in the case of NaNs
            # Compute upward longwave radiation from surface temperature for this simulation
            T_sfc .+=
                area_fraction .*
                ifelse.(
                    area_fraction .≈ 0,
                    zero(T_sfc),
                    emissivity_sim .*
                    Interfacer.get_field(sim, field_name, boundary_space) .^ FT(4),
                )
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
    import_combined_surface_fields!(cs.fields, cs.model_sims, cs.thermo_params)

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
"""
function set_caches!(cs::Interfacer.CoupledSimulation)
    Interfacer.set_cache!(cs.model_sims.atmos_sim)
    exchange!(cs)
    for sim in cs.model_sims
        sim isa Interfacer.SurfaceModelSimulation && Interfacer.set_cache!(sim)
    end
    return nothing
end

end # module
