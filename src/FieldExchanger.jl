"""
    FieldExchanger

This modules contains general functions for the exchange of fields between the
atmospheric and surface component models.
"""
module FieldExchanger

import SciMLBase: step!, reinit!
export import_atmos_fields!,
    import_combined_surface_fields!, update_sim!, update_model_sims!, reinit_model_sims!, step_model_sims!

using ClimaCoupler: Interfacer, FluxCalculator, Regridder, Utilities

"""
    import_atmos_fields!(csf, model_sims, boundary_space, turbulent_fluxes)

Updates the coupler with the atmospheric fluxes. The `Interfacer.get_field` functions
(`:turbulent_energy_flux`, `:turbulent_moisture_flux`, `:radiative_energy_flux`, `:liquid_precipitation`, `:snow_precipitation`)
have to be defined for the amtospheric component model type.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `boundary_space`: [Spaces.AbstractSpace] the space of the coupler surface.
- `turbulent_fluxes`: [TurbulentFluxPartition] denotes a flag for turbulent flux calculation.
"""
function import_atmos_fields!(csf, model_sims, boundary_space, turbulent_fluxes)
    (; atmos_sim) = model_sims

    # from atmos
    if turbulent_fluxes == FluxCalculator.CombinedAtmosGrid()
        Regridder.dummmy_remap!(csf.F_turb_energy, Interfacer.get_field(atmos_sim, Val(:turbulent_energy_flux)))
        Regridder.dummmy_remap!(csf.F_turb_moisture, Interfacer.get_field(atmos_sim, Val(:turbulent_moisture_flux)))
    end

    Regridder.dummmy_remap!(csf.F_radiative, Interfacer.get_field(atmos_sim, Val(:radiative_energy_flux)))
    Regridder.dummmy_remap!(csf.P_liq, Interfacer.get_field(atmos_sim, Val(:liquid_precipitation)))
    Regridder.dummmy_remap!(csf.P_snow, Interfacer.get_field(atmos_sim, Val(:snow_precipitation))) # generalize
end

"""
    import_combined_surface_fields!(csf, model_sims, boundary_space, turbulent_fluxes)

Updates the coupler with the surface properties. The `Interfacer.get_field` functions for
(`:surface_temperature`, `:albedo`, `:roughness_momentum`, `:roughness_buoyancy`, `:beta`)
need to be specified for each surface model.

# Arguments
- `csf`: [NamedTuple] containing coupler fields.
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `boundary_space`: [Spaces.AbstractSpace] the space of the coupler surface.
- `turbulent_fluxes`: [TurbulentFluxPartition] denotes a flag for turbulent flux calculation.

"""
function import_combined_surface_fields!(csf, model_sims, boundary_space, turbulent_fluxes)

    combined_field = zeros(boundary_space)

    # surface fields
    Regridder.combine_surfaces!(combined_field, model_sims, Val(:surface_temperature))
    Regridder.dummmy_remap!(csf.T_S, combined_field)

    Regridder.combine_surfaces!(combined_field, model_sims, Val(:albedo))
    Regridder.dummmy_remap!(csf.albedo, combined_field)

    if turbulent_fluxes !== FluxCalculator.PartitionedComponentModelGrid()
        Regridder.combine_surfaces!(combined_field, model_sims, Val(:roughness_momentum))
        Regridder.dummmy_remap!(csf.z0m_S, combined_field)

        Regridder.combine_surfaces!(combined_field, model_sims, Val(:roughness_buoyancy))
        Regridder.dummmy_remap!(csf.z0b_S, combined_field)

        Regridder.combine_surfaces!(combined_field, model_sims, Val(:beta))
        Regridder.dummmy_remap!(csf.beta, combined_field)
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

    Interfacer.update_field!(atmos_sim, Val(:albedo), csf.albedo)
    Interfacer.update_field!(atmos_sim, Val(:surface_temperature), csf.T_S)

    if turbulent_fluxes == FluxCalculator.CombinedAtmosGrid()
        Interfacer.update_field!(atmos_sim, Val(:roughness_momentum), csf.z0m_S)
        Interfacer.update_field!(atmos_sim, Val(:roughness_buoyancy), csf.z0b_S)
        Interfacer.update_field!(atmos_sim, Val(:beta), csf.beta) # not in this version of atmos
    end

end

"""
    update_sim!(sim::SurfaceModelSimulation, csf, area_fraction = nothing)

Updates the surface component model cache with the current coupler fields of F_turb_energy, F_radiative, F_turb_moisture, P_liq, and ρ_sfc.

# Arguments
- `sim`: [Interfacer.SurfaceModelSimulation] containing a surface model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(sim::Interfacer.SurfaceModelSimulation, csf, area_fraction = nothing)

    Interfacer.update_field!(sim, Val(:air_density), csf.ρ_sfc)

    Interfacer.update_field!(
        sim,
        Val(:turbulent_energy_flux),
        Regridder.binary_mask.(area_fraction, threshold = eps()) .* csf.F_turb_energy,
    )

    Interfacer.update_field!(
        sim,
        Val(:turbulent_moisture_flux),
        Regridder.binary_mask.(area_fraction, threshold = eps()) .* csf.F_turb_moisture,
    )

    Interfacer.update_field!(
        sim,
        Val(:radiative_energy_flux),
        Regridder.binary_mask.(area_fraction, threshold = eps()) .* csf.F_radiative,
    )

    Interfacer.update_field!(sim, Val(:liquid_precipitation), .-csf.P_liq)
    Interfacer.update_field!(sim, Val(:snow_precipitation), .-csf.P_snow)
end

"""
    update_sim!(::SurfaceStub, csf, area_fraction)

The stub surface simulation is not updated by this function.
"""
update_sim!(::Interfacer.SurfaceStub, csf, area_fraction) = nothing

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
            update_sim!(sim, csf, Interfacer.get_field(sim, Val(:area_fraction)))
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
        reinit!(sim)
    end
end

"""
    reinit!(cs::SurfaceStub)

The stub surface simulation is not updated by this function. Extends `SciMLBase.reinit!`.
"""
reinit!(::Interfacer.SurfaceStub) = nothing

"""
    step_model_sims!(model_sims, t)

Iterates `step!` over all component model simulations saved in `cs.model_sims`.

# Arguments
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `t`: [AbstractFloat] denoting the simulation time.
"""
function step_model_sims!(model_sims, t)
    for sim in model_sims
        step!(sim, t)
    end
end

"""
    step!(::SurfaceStub, t)

The stub surface simulation is not updated by this function. Extends `SciMLBase.step!`.
"""
step!(::Interfacer.SurfaceStub, _) = nothing


end # module
