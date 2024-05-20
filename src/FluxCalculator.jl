"""
    FluxCalculator

This modules contains abstract types and functions to calculate surface fluxes in the coupler,
or to call flux calculating functions from the component models.
"""

module FluxCalculator

import StaticArrays
import SurfaceFluxes as SF
import Thermodynamics as TD
import ClimaCore as CC
import ..Interfacer, ..Regridder

export PartitionedStateFluxes,
    CombinedStateFluxesMOST,
    combined_turbulent_fluxes!,
    TurbulentFluxPartition,
    atmos_turbulent_fluxes_most!,
    calculate_surface_air_density,
    extrapolate_ρ_to_sfc,
    MoninObukhovScheme,
    BulkScheme,
    partitioned_turbulent_fluxes!,
    get_surface_params,
    update_turbulent_fluxes!,
    water_albedo_from_atmosphere!

"""
    TurbulentFluxPartition

Abstract type for flags that denote where and how to calculate tubulent fluxes.
"""
abstract type TurbulentFluxPartition end

"""
    PartitionedStateFluxes <: TurbulentFluxPartition

A flag indicating that the turbulent fluxes should be partitioned and calculated
over each surface model and then combined. This is calculated on the coupler grid.
"""
struct PartitionedStateFluxes <: TurbulentFluxPartition end

"""
    CombinedStateFluxesMOST <: TurbulentFluxPartition

A flag indicating that the turbulent fluxes (e.g. sensible and latent heat fluxes,
drag and moisture fluxes) are to be  calculated on the Atmos grid, and saved in Atmos cache.
"""
struct CombinedStateFluxesMOST <: TurbulentFluxPartition end

"""
    combined_turbulent_fluxes!(model_sims, csf, turbulent_fluxes::TurbulentFluxPartition)

Calls the method(s) which calculate turbulent surface fluxes from combined surface states in coupler fields, `csf`.

# Arguments
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `csf`: [NamedTuple] containing coupler fields.
- `turbulent_fluxes`: [TurbulentFluxPartition] denotes a flag for turbulent flux calculation.

"""
function combined_turbulent_fluxes!(model_sims, csf, turbulent_fluxes::TurbulentFluxPartition)
    if turbulent_fluxes isa CombinedStateFluxesMOST
        atmos_turbulent_fluxes_most!(model_sims.atmos_sim, csf)
    else
        nothing # TODO: may want to add CombinedCouplerGrid
    end
end

"""
    atmos_turbulent_fluxes_most!(sim::Interfacer.ComponentModelSimulation, csf)

A function to calculate turbulent surface fluxes using the combined surface states
and the Monin Obukhov Similarity Theory.
It is required that a method is defined for the given `sim` and that the fluxes are
saved in that sim's cache. `csf` refers to the coupler fields.

# Arguments
- `sim`: [Interfacer.ComponentModelSimulation] object containing the component model simulation.
- `csf`: [NamedTuple] containing coupler fields.

# Example:

```
function atmos_turbulent_fluxes_most!(atmos_sim::ClimaAtmosSimulation, csf)
    atmos_sim.cache.flux .= atmos_sim.c .* (csf.T_S .- atmos_sim.temperature)
end
```

"""
atmos_turbulent_fluxes_most!(sim::Interfacer.ComponentModelSimulation, _) =
    error("calling flux calculation in " * Interfacer.name(sim) * " but no method defined")


"""
    calculate_surface_air_density(atmos_sim::ClimaAtmosSimulation, T_S::CC.Fields.Field)

Extension for this  to to calculate surface density.
"""
function calculate_surface_air_density(atmos_sim::Interfacer.AtmosModelSimulation, T_S::CC.Fields.Field)
    error("this function is required to be dispatched on" * Interfacer.name(atmos_sim) * ", but no method defined")
end

"""
    partitioned_turbulent_fluxes!(model_sims::NamedTuple, fields::NamedTuple, boundary_space::CC.Spaces.AbstractSpace, surface_scheme, thermo_params::TD.Parameters.ThermodynamicsParameters)

The current setup calculates the aerodynamic fluxes in the coupler (assuming no regridding is needed)
using adapter function `get_surface_fluxes!`, which calls `SurfaceFluxes.jl`. The coupler saves
the area-weighted sums of the fluxes.

Args:
- `model_sims`: [NamedTuple] containing `ComponentModelSimulation`s.
- `fields`: [NamedTuple] containing coupler fields.
- `boundary_space`: [CC.Spaces.AbstractSpace] the space of the coupler surface.
- `surface_scheme`: [AbstractSurfaceFluxScheme] the surface flux scheme.
- `thermo_params`: [TD.Parameters.ThermodynamicsParameters] the thermodynamic parameters.

TODO:
- generalize interface for regridding and take land state out of atmos's integrator.p
- add flux accumulation
- add flux bounds

(NB: Radiation surface fluxes are calculated by the atmosphere.)

"""
function partitioned_turbulent_fluxes!(
    model_sims::NamedTuple,
    fields::NamedTuple,
    boundary_space::CC.Spaces.AbstractSpace,
    surface_scheme,
    thermo_params::TD.Parameters.ThermodynamicsParameters,
)

    atmos_sim = model_sims.atmos_sim
    csf = fields
    FT = eltype(csf[1])

    # reset coupler fields
    csf.F_turb_ρτxz .*= FT(0)
    csf.F_turb_ρτyz .*= FT(0)
    csf.F_turb_energy .*= FT(0)
    csf.F_turb_moisture .*= FT(0)

    # atmos state of center level 1
    z_int = Interfacer.get_field(atmos_sim, Val(:height_int))
    uₕ_int = Interfacer.get_field(atmos_sim, Val(:uv_int))
    thermo_state_int = Interfacer.get_field(atmos_sim, Val(:thermo_state_int))
    z_sfc = Interfacer.get_field(atmos_sim, Val(:height_sfc))

    for sim in model_sims
        # iterate over all surface models
        if sim isa Interfacer.SurfaceModelSimulation
            # get area fraction (min = 0, max = 1)
            area_fraction = Interfacer.get_field(sim, Val(:area_fraction))
            # get area mask [0, 1], where area_mask = 1 if area_fraction > 0
            area_mask = Regridder.binary_mask.(area_fraction)

            thermo_state_sfc = FluxCalculator.surface_thermo_state(sim, thermo_params, thermo_state_int)

            # set inputs based on whether the surface_scheme is `MoninObukhovScheme` or `BulkScheme`
            surface_params = FluxCalculator.get_surface_params(atmos_sim)
            scheme_properties = FluxCalculator.get_scheme_properties(surface_scheme, sim)

            input_args = (;
                thermo_state_sfc,
                thermo_state_int,
                uₕ_int,
                z_int,
                z_sfc,
                scheme_properties,
                surface_params,
                surface_scheme,
                boundary_space,
            )
            inputs = FluxCalculator.surface_inputs(surface_scheme, input_args)

            # calculate the surface fluxes
            fluxes = FluxCalculator.get_surface_fluxes!(inputs, surface_params)
            (; F_turb_ρτxz, F_turb_ρτyz, F_shf, F_lhf, F_turb_moisture) = fluxes

            # perform additional diagnostics if required
            FluxCalculator.differentiate_turbulent_fluxes!(sim, (thermo_params, input_args, fluxes))

            # update fluxes in the coupler
            fields = (;
                F_turb_ρτxz = F_turb_ρτxz,
                F_turb_ρτyz = F_turb_ρτyz,
                F_turb_energy = F_shf .+ F_lhf,
                F_turb_moisture = F_turb_moisture,
            )

            # update the fluxes of this surface model
            FluxCalculator.update_turbulent_fluxes!(sim, fields)

            # add the flux contributing from this surface to the coupler field
            # note that the fluxes are area-weighted, so if a surface model is
            #  not present at this point, the fluxes are zero
            @. csf.F_turb_ρτxz += F_turb_ρτxz * area_fraction * area_mask
            @. csf.F_turb_ρτyz += F_turb_ρτyz * area_fraction * area_mask
            @. csf.F_turb_energy += (F_shf .+ F_lhf) * area_fraction * area_mask
            @. csf.F_turb_moisture += F_turb_moisture * area_fraction * area_mask
        end
    end

    # TODO: add allowable bounds here, check explicitly that all fluxes are equal

end

abstract type AbstractSurfaceFluxScheme end
struct BulkScheme <: AbstractSurfaceFluxScheme end
struct MoninObukhovScheme <: AbstractSurfaceFluxScheme end

"""
    get_scheme_properties(scheme::AbstractSurfaceFluxScheme, sim::Interfacer.SurfaceModelSimulation)

Returns the scheme-specific properties for the surface model simulation `sim`.
"""
function get_scheme_properties(::BulkScheme, sim::Interfacer.SurfaceModelSimulation)
    Ch = Interfacer.get_field(sim, Val(:heat_transfer_coefficient))
    Cd = Interfacer.get_field(sim, Val(:drag_coefficient))
    beta = Interfacer.get_field(sim, Val(:beta))
    FT = eltype(Ch)
    return (; z0b = FT(0), z0m = FT(0), Ch = Ch, Cd = Cd, beta = beta, gustiness = FT(1))
end
function get_scheme_properties(::MoninObukhovScheme, sim::Interfacer.SurfaceModelSimulation)
    z0m = Interfacer.get_field(sim, Val(:roughness_momentum))
    z0b = Interfacer.get_field(sim, Val(:roughness_buoyancy))
    beta = Interfacer.get_field(sim, Val(:beta))
    FT = eltype(z0m)
    return (; z0b = z0b, z0m = z0m, Ch = FT(0), Cd = FT(0), beta = beta, gustiness = FT(1))
end

"""
    surface_inputs(scheme::AbstractSurfaceFluxScheme, thermo_state_sfc, thermo_state_int, uₕ_int, z_int, z_sfc, z0b, z0m, Ch, Cd, beta, gustiness)

Returns the inputs for the surface model simulation `sim`.
"""
function surface_inputs(::BulkScheme, input_args::NamedTuple)

    (; thermo_state_sfc, thermo_state_int, uₕ_int, z_int, z_sfc, scheme_properties, boundary_space) = input_args
    FT = CC.Spaces.undertype(axes(z_sfc))
    (; Ch, Cd, beta, gustiness) = scheme_properties

    # Extract the underlying data layouts of each field
    # Note: this is a bit "dangerous" because it circumvents ClimaCore, but
    #  it allows us to broadcast over fields on slightly different spaces
    fv = CC.Fields.field_values
    z_int_fv = fv(z_int)
    uₕ_int_fv = fv(uₕ_int)
    thermo_state_int_fv = fv(thermo_state_int)
    z_sfc_fv = fv(z_sfc)
    thermo_state_sfc_fv = fv(thermo_state_sfc)
    beta_fv = fv(beta)

    # wrap state values
    result = @. SF.Coefficients(
        SF.StateValues(z_int_fv, uₕ_int_fv, thermo_state_int_fv), # state_in
        SF.StateValues(                                   # state_sfc
            z_sfc_fv,
            StaticArrays.SVector(FT(0), FT(0)),
            thermo_state_sfc_fv,
        ),
        Cd,                                     # Cd
        Ch,                                     # Ch
        gustiness,                              # gustiness
        beta_fv,                                # beta
    )

    # Put the result data layout back onto the surface space
    return CC.Fields.Field(result, boundary_space)
end
function surface_inputs(::MoninObukhovScheme, input_args::NamedTuple)
    (; thermo_state_sfc, thermo_state_int, uₕ_int, z_int, z_sfc, scheme_properties, boundary_space) = input_args
    FT = CC.Spaces.undertype(axes(z_sfc))
    (; z0b, z0m, beta, gustiness) = scheme_properties

    # Extract the underlying data layouts of each field
    # Note: this is a bit "dangerous" because it circumvents ClimaCore, but
    #  it allows us to broadcast over fields on slightly different spaces
    fv = CC.Fields.field_values
    z_int_fv = fv(z_int)
    uₕ_int_fv = uₕ_int isa CC.Fields.Field ? fv(uₕ_int) : uₕ_int
    thermo_state_int_fv = fv(thermo_state_int)
    z_sfc_fv = fv(z_sfc)
    thermo_state_sfc_fv = fv(thermo_state_sfc)
    beta_fv = beta isa CC.Fields.Field ? fv(beta) : beta

    # Compute state values
    result = @. SF.ValuesOnly(
        SF.StateValues(z_int_fv, uₕ_int_fv, thermo_state_int_fv), # state_in
        SF.StateValues(                                  # state_sfc
            z_sfc_fv,
            StaticArrays.SVector(FT(0), FT(0)),
            thermo_state_sfc_fv,
        ),
        z0m,                                   # z0m
        z0b,                                   # z0b
        gustiness,                             # gustiness
        beta_fv,                               # beta
    )

    # Put the result data layout back onto the surface space
    return CC.Fields.Field(result, boundary_space)
end

"""
    surface_thermo_state(sim::Interfacer.SurfaceModelSimulation, thermo_params::TD.Parameters.ThermodynamicsParameters, thermo_state_int)

Returns the surface parameters for the surface model simulation `sim`. The default is assuming saturated surfaces, unless an extension is defined for the given `SurfaceModelSimulation`.
"""
function surface_thermo_state(
    sim::Interfacer.SurfaceModelSimulation,
    thermo_params::TD.Parameters.ThermodynamicsParameters,
    thermo_state_int;
    δT_sfc = 0,
)
    FT = eltype(parent(thermo_state_int))
    @warn("Simulation " * Interfacer.name(sim) * " uses the default thermo (saturated) surface state", maxlog = 10)
    # get surface temperature (or perturbed surface temperature for differentiation)
    T_sfc = Interfacer.get_field(sim, Val(:surface_temperature)) .+ FT(δT_sfc)
    ρ_sfc = extrapolate_ρ_to_sfc.(thermo_params, thermo_state_int, T_sfc) # ideally the # calculate elsewhere, here just getter...
    q_sfc = TD.q_vap_saturation_generic.(thermo_params, T_sfc, ρ_sfc, TD.Liquid()) # default: saturated liquid surface
    @. TD.PhaseEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, q_sfc)
end


# TODO: (an equivalent of) this function also lives in Atmos and Land - should move to general utilities
"""
    extrapolate_ρ_to_sfc(thermo_params, ts_int, T_sfc)

Uses the ideal gas law and hydrostatic balance to extrapolate for surface density.
"""
function extrapolate_ρ_to_sfc(thermo_params, ts_in, T_sfc)
    T_int = TD.air_temperature(thermo_params, ts_in)
    Rm_int = TD.gas_constant_air(thermo_params, ts_in)
    ρ_air = TD.air_density(thermo_params, ts_in)
    ρ_air * (T_sfc / T_int)^(TD.cv_m(thermo_params, ts_in) / Rm_int)
end

"""
    get_surface_fluxes!(inputs, surface_params::SF.Parameters.SurfaceFluxesParameters)

Uses SurfaceFluxes.jl to calculate turbulent surface fluxes. It should be atmos model agnostic, and columnwise.
"""
function get_surface_fluxes!(inputs, surface_params::SF.Parameters.SurfaceFluxesParameters)
    # calculate all fluxes (saturated surface conditions)
    outputs = SF.surface_conditions.(surface_params, inputs)

    # drag
    F_turb_ρτxz = outputs.ρτxz
    F_turb_ρτyz = outputs.ρτyz

    # energy fluxes
    F_shf = outputs.shf
    F_lhf = outputs.lhf

    # moisture
    F_turb_moisture = SF.evaporation.(surface_params, inputs, outputs.Ch)

    # At locations where this surface model is not evaluated, we get `NaN` for
    # surface fluxes. In that case, we replace the values with 0.
    @. F_turb_ρτxz = ifelse(isnan(F_turb_ρτxz), zero(F_turb_ρτxz), F_turb_ρτxz)
    @. F_turb_ρτyz = ifelse(isnan(F_turb_ρτyz), zero(F_turb_ρτyz), F_turb_ρτyz)
    @. F_shf = ifelse(isnan(F_shf), zero(F_shf), F_shf)
    @. F_lhf = ifelse(isnan(F_lhf), zero(F_lhf), F_lhf)
    @. F_turb_moisture = ifelse(isnan(F_turb_moisture), zero(F_turb_moisture), F_turb_moisture)

    return (;
        F_turb_ρτxz = F_turb_ρτxz,
        F_turb_ρτyz = F_turb_ρτyz,
        F_shf = F_shf,
        F_lhf = F_lhf,
        F_turb_moisture = F_turb_moisture,
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
        "get_surface_params is required to be dispatched on" * Interfacer.name(atmos_sim) * ", but no method defined",
    )
end

"""
    update_turbulent_fluxes!(sim::Interfacer.SurfaceModelSimulation, fields::NamedTuple)

Updates the fluxes in the surface model simulation `sim` with the fluxes in `fields`.
"""
function update_turbulent_fluxes!(sim::Interfacer.SurfaceModelSimulation, fields::NamedTuple)
    return error(
        "update_turbulent_fluxes! is required to be dispatched on" * Interfacer.name(sim) * ", but no method defined",
    )
end

update_turbulent_fluxes!(sim::Interfacer.SurfaceStub, fields::NamedTuple) = nothing

"""
    differentiate_turbulent_fluxes!(sim::Interfacer.SurfaceModelSimulation, args)

This function provides a placeholder for differentiating fluxes with respect to
surface temperature in surface energy balance calculations.
"""
differentiate_turbulent_fluxes!(::Interfacer.SurfaceModelSimulation, args) = nothing

"""
    water_albedo_from_atmosphere!(cs::Interfacer.CoupledSimulation, _)

Callback to calculate the water albedo from atmospheric state. This is a placeholder for the full radiation callback.
"""
function water_albedo_from_atmosphere!(cs::Interfacer.CoupledSimulation, _)
    atmos_sim = cs.model_sims.atmos_sim
    ocean_sim = cs.model_sims.ocean_sim
    cf = cs.fields

    # use temp fields
    water_albedo_from_atmosphere!(atmos_sim, cf.temp1, cf.temp2)

    Interfacer.update_field!(ocean_sim, Val(:surface_direct_albedo), cf.temp1)
    Interfacer.update_field!(ocean_sim, Val(:surface_diffuse_albedo), cf.temp2)

end

"""
    water_albedo_from_atmosphere!(atmos_sim::Interfacer.AtmosModelSimulation, ::CC.Fields.Field, ::CC.Fields.Field)

Placeholder for the water albedo calculation from the atmosphere. It returns an error if not extended.
"""
function water_albedo_from_atmosphere!(atmos_sim::Interfacer.AtmosModelSimulation, ::CC.Fields.Field, ::CC.Fields.Field)
    error("this function is required to be dispatched on" * Interfacer.name(atmos_sim) * ", but no method defined")
end

end # module
