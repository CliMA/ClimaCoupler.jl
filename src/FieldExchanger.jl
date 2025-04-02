"""
    FieldExchanger

This modules contains general functions for the exchange of fields between the
atmospheric and surface component models.
"""
module FieldExchanger

import ..Interfacer, ..FluxCalculator, ..Utilities

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
    import_atmos_fields!(csf, model_sims, boundary_space, turbulent_fluxes)

Update the coupler with the atmospheric fluxes. The `Interfacer.get_field` functions
(`:turbulent_energy_flux`, `:turbulent_moisture_flux`, `:radiative_energy_flux_sfc`, `:liquid_precipitation`, `:snow_precipitation`)
have to be defined for the amtospheric component model type.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `boundary_space`: [Spaces.AbstractSpace] the space of the coupler surface.
- `turbulent_fluxes`: [TurbulentFluxPartition] denotes a flag for turbulent flux calculation.
"""
function import_atmos_fields!(csf, model_sims, boundary_space, turbulent_fluxes)
    for sim in model_sims
        import_atmos_fields!(csf, sim, model_sims.atmos_sim, turbulent_fluxes)
    end
end

"""
    import_atmos_fields!(csf, ::Interfacer.SurfaceModelSimulation, atmos_sim, turbulent_fluxes)

Updates the coupler simulation fields with atmospheric fluxes from the atmosphere simulation.
This is the default function to be used for most surface model simulations, as
    are computed by the coupler or atmosphere
and passed to the surfaces.
"""
function import_atmos_fields!(csf, ::Interfacer.SurfaceModelSimulation, atmos_sim, turbulent_fluxes)
    # turbulent fluxes
    if turbulent_fluxes isa FluxCalculator.CombinedStateFluxesMOST
        dummmy_remap!(csf.F_turb_energy, Interfacer.get_field(atmos_sim, Val(:turbulent_energy_flux)))
        dummmy_remap!(csf.F_turb_moisture, Interfacer.get_field(atmos_sim, Val(:turbulent_moisture_flux)))
    end

    # surface density - needed for q_sat and requires atmos and sfc states, so it is calculated and saved in the coupler
    dummmy_remap!(csf.ρ_sfc, FluxCalculator.calculate_surface_air_density(atmos_sim, csf.T_sfc)) # TODO: generalize for PartitionedStateFluxes (#445) (use individual T_sfc)

    # radiative fluxes
    dummmy_remap!(csf.F_radiative, Interfacer.get_field(atmos_sim, Val(:radiative_energy_flux_sfc)))

    # precipitation
    dummmy_remap!(csf.P_liq, Interfacer.get_field(atmos_sim, Val(:liquid_precipitation)))
    dummmy_remap!(csf.P_snow, Interfacer.get_field(atmos_sim, Val(:snow_precipitation)))
end

import_atmos_fields!(csf, ::Interfacer.AtmosModelSimulation, atmos_sim, turbulent_fluxes) = nothing

"""
    import_combined_surface_fields!(csf, model_sims, turbulent_fluxes)

Updates the coupler with the surface properties. The `Interfacer.get_field` functions for
(`:surface_temperature`, `:surface_direct_albedo`, `:surface_diffuse_albedo`, `:roughness_momentum`, `:roughness_buoyancy`, `:beta`)
need to be specified for each surface model.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `turbulent_fluxes`: [TurbulentFluxPartition] denotes a flag for turbulent flux calculation.

"""
function import_combined_surface_fields!(csf, model_sims, turbulent_fluxes)

    combined_field = csf.temp1

    # surface fields
    combine_surfaces!(combined_field, model_sims, Val(:surface_temperature))
    dummmy_remap!(csf.T_sfc, combined_field)

    combine_surfaces!(combined_field, model_sims, Val(:surface_direct_albedo))
    dummmy_remap!(csf.surface_direct_albedo, combined_field)

    combine_surfaces!(combined_field, model_sims, Val(:surface_diffuse_albedo))
    dummmy_remap!(csf.surface_diffuse_albedo, combined_field)

    if turbulent_fluxes isa FluxCalculator.CombinedStateFluxesMOST
        combine_surfaces!(combined_field, model_sims, Val(:roughness_momentum))
        dummmy_remap!(csf.z0m_sfc, combined_field)

        combine_surfaces!(combined_field, model_sims, Val(:roughness_buoyancy))
        dummmy_remap!(csf.z0b_sfc, combined_field)

        combine_surfaces!(combined_field, model_sims, Val(:beta))
        dummmy_remap!(csf.beta, combined_field)

        combine_surfaces!(combined_field, model_sims, Val(:surface_humidity))
        dummmy_remap!(csf.q_sfc, combined_field)
    end

end

"""
    update_sim!(atmos_sim::Interfacer.AtmosModelSimulation, csf)

Updates the surface fields for temperature, roughness length, albedo, and specific humidity.

# Arguments
- `atmos_sim`: [Interfacer.AtmosModelSimulation] containing an atmospheric model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(atmos_sim::Interfacer.AtmosModelSimulation, csf, turbulent_fluxes)

    Interfacer.update_field!(atmos_sim, Val(:surface_direct_albedo), csf.surface_direct_albedo)
    Interfacer.update_field!(atmos_sim, Val(:surface_diffuse_albedo), csf.surface_diffuse_albedo)
    Interfacer.update_field!(atmos_sim, Val(:surface_temperature), csf.T_sfc)

    if turbulent_fluxes isa FluxCalculator.CombinedStateFluxesMOST
        Interfacer.update_field!(atmos_sim, Val(:roughness_momentum), csf.z0m_sfc)
        Interfacer.update_field!(atmos_sim, Val(:roughness_buoyancy), csf.z0b_sfc)
        Interfacer.update_field!(atmos_sim, Val(:beta), csf.beta) # not in this version of atmos
    end

    if turbulent_fluxes isa FluxCalculator.PartitionedStateFluxes
        Interfacer.update_field!(atmos_sim, Val(:turbulent_fluxes), csf)
    end
end

"""
    update_sim!(sim::SurfaceModelSimulation, csf, area_fraction)

Updates the surface component model cache with the current coupler fields of F_turb_energy, F_radiative, F_turb_moisture, P_liq, and ρ_sfc.

# Arguments
- `sim`: [Interfacer.SurfaceModelSimulation] containing a surface model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(sim::Interfacer.SurfaceModelSimulation, csf, turbulent_fluxes, area_fraction)
    FT = eltype(area_fraction)

    # atmospheric surface density
    Interfacer.update_field!(sim, Val(:air_density), csf.ρ_sfc)

    # turbulent fluxes
    mask = Utilities.binary_mask.(area_fraction)

    # when PartitionedStateFluxes, turbulent fluxes are updated during the flux calculation
    if turbulent_fluxes isa FluxCalculator.CombinedStateFluxesMOST
        F_turb_energy_masked = FT.(mask .* csf.F_turb_energy)
        Interfacer.update_field!(sim, Val(:turbulent_energy_flux), F_turb_energy_masked)
        Interfacer.update_field!(sim, Val(:turbulent_moisture_flux), FT.(mask .* csf.F_turb_moisture))
    end

    # radiative fluxes
    Interfacer.update_field!(sim, Val(:radiative_energy_flux_sfc), FT.(mask .* csf.F_radiative))

    # precipitation
    Interfacer.update_field!(sim, Val(:liquid_precipitation), csf.P_liq)
    Interfacer.update_field!(sim, Val(:snow_precipitation), csf.P_snow)
end

"""
    update_model_sims!(model_sims, csf, turbulent_fluxes)

Iterates `update_sim!` over all component model simulations saved in `cs.model_sims`.

# Arguments
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `csf`: [NamedTuple] containing coupler fields.
- `turbulent_fluxes`: [TurbulentFluxPartition] denotes a flag for turbulent flux calculation.

"""
function update_model_sims!(model_sims, csf, turbulent_fluxes)
    for sim in model_sims
        if sim isa Interfacer.SurfaceModelSimulation
            update_sim!(sim, csf, turbulent_fluxes, Interfacer.get_field(sim, Val(:area_fraction)))
        else
            update_sim!(sim, csf, turbulent_fluxes)
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
    nans_to_zero(v)

Replaces NaNs with zeros, otherwise returns the value.
"""
nans_to_zero(v) = isnan(v) ? typeof(v)(0) : v

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
    combined_field .= eltype(combined_field)(0)
    for sim in sims
        if sim isa Interfacer.SurfaceModelSimulation
            combined_field .+=
                Interfacer.get_field(sim, Val(:area_fraction)) .* nans_to_zero.(Interfacer.get_field(sim, field_name)) # this ensures that unitialized (masked) areas do not affect (TODO: move to mask / remove)
        end
    end
end

end # module
