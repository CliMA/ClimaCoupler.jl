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
import ..Interfacer, ..Utilities

export calculate_surface_air_density,
    extrapolate_ρ_to_sfc,
    MoninObukhovScheme,
    BulkScheme,
    turbulent_fluxes!,
    get_surface_params,
    update_turbulent_fluxes!,
    water_albedo_from_atmosphere!,
    compute_surface_fluxes!

"""
    calculate_surface_air_density(atmos_sim::ClimaAtmosSimulation, T_sfc::CC.Fields.Field)

Extension for this  to to calculate surface density.
"""
function calculate_surface_air_density(atmos_sim::Interfacer.AtmosModelSimulation, T_sfc::CC.Fields.Field)
    error("this function is required to be dispatched on" * Interfacer.name(atmos_sim) * ", but no method defined")
end

"""
    turbulent_fluxes!(model_sims::NamedTuple,
                                  csf::CC.Fields.Field,
                                  boundary_space::CC.Spaces.AbstractSpace,
                                  surface_scheme,
                                  thermo_params::TD.Parameters.ThermodynamicsParameters)

Compute turbulent fluxes and associated quantities. Store the results in `fields` as
area-weighted sums.

This function uses `SurfaceFluxes.jl` under the hood.

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
function turbulent_fluxes!(
    model_sims::NamedTuple,
    csf::CC.Fields.Field,
    boundary_space::CC.Spaces.AbstractSpace,
    surface_scheme,
    thermo_params::TD.Parameters.ThermodynamicsParameters,
)
    atmos_sim = model_sims.atmos_sim
    FT = CC.Spaces.undertype(boundary_space)

    # Reset the coupler fields will compute. We need to do this because we will compute
    # area-weighted averages
    csf.F_turb_ρτxz .*= FT(0)
    csf.F_turb_ρτyz .*= FT(0)
    csf.F_turb_energy .*= FT(0)
    csf.F_turb_moisture .*= FT(0)
    csf.z0m_sfc .*= FT(0)
    csf.z0b_sfc .*= FT(0)
    csf.beta .*= FT(0)
    csf.q_sfc .*= FT(0)
    csf.L_MO .*= FT(0)
    csf.ustar .*= FT(0)
    csf.buoyancy_flux .*= FT(0)

    # Compute the surface fluxes for each surface model and add them to `csf`
    for sim in model_sims
        compute_surface_fluxes!(csf, sim, atmos_sim, boundary_space, thermo_params, surface_scheme)
    end

    # TODO: add allowable bounds here, check explicitly that all fluxes are equal
    return nothing
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
    FT = CC.Spaces.undertype(boundary_space)
    (; z0b, z0m, beta, gustiness) = scheme_properties

    # Extract the underlying data layouts of each field
    # Note: this is a bit "dangerous" because it circumvents ClimaCore, but
    #  it allows us to broadcast over fields on slightly different spaces
    maybe_fv = (x) -> x isa CC.Fields.Field ? CC.Fields.field_values(x) : x

    z_int_fv = maybe_fv(z_int)
    uₕ_int_fv = maybe_fv(uₕ_int)
    thermo_state_int_fv = maybe_fv(thermo_state_int)
    z_sfc_fv = maybe_fv(z_sfc)
    thermo_state_sfc_fv = maybe_fv(thermo_state_sfc)
    beta_fv = maybe_fv(beta)
    z0m_fv = maybe_fv(z0m)
    z0b_fv = maybe_fv(z0b)
    gustiness_fv = maybe_fv(gustiness)

    # Compute state values
    result = @. SF.ValuesOnly(
        SF.StateValues(z_int_fv, uₕ_int_fv, thermo_state_int_fv), # state_in
        SF.StateValues(                                  # state_sfc
            z_sfc_fv,
            StaticArrays.SVector(FT(0), FT(0)),
            thermo_state_sfc_fv,
        ),
        z0m_fv,
        z0b_fv,
        gustiness_fv,
        beta_fv,
    )

    # Put the result data layout back onto the surface space
    return CC.Fields.Field(result, boundary_space)
end

"""
    surface_thermo_state(sim::Interfacer.SurfaceModelSimulation,
                         thermo_params::TD.Parameters.ThermodynamicsParameters,
                         thermo_state_int)

Return the surface thermo state the surface model simulation `sim`.

This is obtained by using the model surface temperature in `sim`, extrapolating atmospheric
density adiabatically to the surface, and using the model surface humidity (which, by
default, is computed assuming a liquid phase).
"""
function surface_thermo_state(
    sim::Interfacer.SurfaceModelSimulation,
    thermo_params::TD.Parameters.ThermodynamicsParameters,
    thermo_state_int,
    δT_sfc = 0,
)
    FT = eltype(parent(thermo_state_int))

    # get surface temperature (or perturbed surface temperature for differentiation)
    T_sfc = Interfacer.get_field(sim, Val(:surface_temperature)) .+ FT(δT_sfc)
    # Note that the surface air density, ρ_sfc, is computed using the atmospheric state at the first level and making ideal gas
    # and hydrostatic balance assumptions. The land model does not compute the surface air density so this is
    # a reasonable stand-in.
    #
    # NOTE: This allocates! Fix me!
    ρ_sfc = FluxCalculator.extrapolate_ρ_to_sfc.(thermo_params, thermo_state_int, T_sfc) # ideally the # calculate elsewhere, here just getter...

    # For SurfaceStabs, this is just liquid phase
    q_sfc = TD.q_vap_saturation_generic.(thermo_params, T_sfc, ρ_sfc, TD.Liquid()) # default: saturated liquid surface

    # NOTE: This allocates! Fix me!
    return @. TD.PhaseEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, q_sfc)
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
    return ρ_air * (T_sfc / T_int)^(TD.cv_m(thermo_params, ts_in) / Rm_int)
end

"""
    get_surface_fluxes!(inputs, surface_params::SF.Parameters.SurfaceFluxesParameters)

Uses SurfaceFluxes.jl to calculate turbulent surface fluxes. It should be atmos model agnostic, and columnwise.
Fluxes are computed over the entire surface, even where the relevant surface model is not present.

When available, it also computes ancillary quantities, such as the Monin-Obukov lengthscale.
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

    L_MO = outputs.L_MO
    ustar = outputs.ustar
    buoyancy_flux = outputs.buoy_flux

    return (; F_turb_ρτxz, F_turb_ρτyz, F_shf, F_lhf, F_turb_moisture, L_MO, ustar, buoyancy_flux)
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

update_turbulent_fluxes!(sim::Interfacer.AbstractSurfaceStub, fields::NamedTuple) = nothing

"""
    differentiate_turbulent_fluxes!(sim::Interfacer.SurfaceModelSimulation, args)

This function provides a placeholder for differentiating fluxes with respect to
surface temperature in surface energy balance calculations.
"""
differentiate_turbulent_fluxes!(::Interfacer.SurfaceModelSimulation, args) = nothing

"""
    water_albedo_from_atmosphere!(cs::Interfacer.CoupledSimulation)

Callback to calculate the water albedo from atmospheric state. This is a placeholder for the full radiation callback.
"""
function water_albedo_from_atmosphere!(cs::Interfacer.CoupledSimulation)
    atmos_sim = cs.model_sims.atmos_sim
    ocean_sim = cs.model_sims.ocean_sim
    cf = cs.fields

    # use temp fields
    water_albedo_from_atmosphere!(atmos_sim, cf.temp1, cf.temp2)

    Interfacer.update_field!(ocean_sim, Val(:surface_direct_albedo), cf.temp1)
    Interfacer.update_field!(ocean_sim, Val(:surface_diffuse_albedo), cf.temp2)
    return nothing
end

"""
    water_albedo_from_atmosphere!(atmos_sim::Interfacer.AtmosModelSimulation, ::CC.Fields.Field, ::CC.Fields.Field)

Placeholder for the water albedo calculation from the atmosphere. It returns an error if not extended.
"""
function water_albedo_from_atmosphere!(atmos_sim::Interfacer.AtmosModelSimulation, ::CC.Fields.Field, ::CC.Fields.Field)
    error("this function is required to be dispatched on" * Interfacer.name(atmos_sim) * ", but no method defined")
end

"""
    compute_surface_fluxes!(csf, sim, atmos_sim, boundary_space, thermo_params, surface_scheme)

This function computes surface fluxes between the input component model
simulation and the atmosphere.

Update the input coupler surface fields `csf` in-place with the computed fluxes
for this model. These are then summed using area-weighting across all surface
models to get the total fluxes.

Since the fluxes are computed between the input model and the atmosphere, this
function does nothing if called on an atmosphere model simulation.

# Arguments
- `csf`: [CC.Fields.Field] containing a NamedTuple of turbulent flux fields: `F_turb_ρτxz`, `F_turb_ρτyz`, `F_turb_energy`, `F_turb_moisture`.
- `sim`: [Interfacer.ComponentModelSimulation] the surface simulation to compute fluxes for.
- `atmos_sim`: [Interfacer.AtmosModelSimulation] the atmosphere simulation to compute fluxes with.
- `boundary_space`: [CC.Spaces.AbstractSpace] the space of the coupler surface.
- `thermo_params`: [TD.Parameters.ThermodynamicsParameters] the thermodynamic parameters.
- `surface_scheme`: [AbstractSurfaceFluxScheme] the surface flux scheme.
"""
function compute_surface_fluxes!(
    csf,
    sim::Interfacer.AtmosModelSimulation,
    atmos_sim::Interfacer.AtmosModelSimulation,
    boundary_space,
    thermo_params,
    surface_scheme,
)
    # do nothing for atmos model
    return nothing
end

function compute_surface_fluxes!(
    csf,
    sim::Interfacer.SurfaceModelSimulation,
    atmos_sim::Interfacer.AtmosModelSimulation,
    boundary_space,
    thermo_params,
    surface_scheme,
)
    # `_int` refers to atmos state of center level 1
    z_int = Interfacer.get_field(atmos_sim, Val(:height_int))
    uₕ_int = Interfacer.get_field(atmos_sim, Val(:uv_int))
    thermo_state_int = Interfacer.get_field(atmos_sim, Val(:thermo_state_int))
    z_sfc = Interfacer.get_field(atmos_sim, Val(:height_sfc))

    # get area fraction (min = 0, max = 1)
    area_fraction = Interfacer.get_field(sim, Val(:area_fraction))

    thermo_state_sfc = FluxCalculator.surface_thermo_state(sim, thermo_params, thermo_state_int)

    # set inputs based on whether the surface_scheme is `MoninObukhovScheme` or `BulkScheme`
    surface_params = FluxCalculator.get_surface_params(atmos_sim)
    scheme_properties = FluxCalculator.get_scheme_properties(surface_scheme, sim)

    # surface_params, surface_scheme are used by EinsenmanIceSimulation
    input_args = (;
        thermo_state_sfc,
        thermo_state_int,
        uₕ_int,
        z_int,
        z_sfc,
        scheme_properties,
        boundary_space,
        surface_params,
        surface_scheme,
    )
    inputs = FluxCalculator.surface_inputs(surface_scheme, input_args)

    # calculate the surface fluxes
    fluxes = FluxCalculator.get_surface_fluxes!(inputs, surface_params)
    (; F_turb_ρτxz, F_turb_ρτyz, F_shf, F_lhf, F_turb_moisture, L_MO, ustar, buoyancy_flux) = fluxes


    # Zero out fluxes where the area fraction is zero
    # Multiplying by `area_fraction` is not sufficient because the fluxes may
    # be NaN where the area fraction is zero.
    @. F_turb_ρτxz = ifelse(area_fraction ≈ 0, zero(F_turb_ρτxz), F_turb_ρτxz)
    @. F_turb_ρτyz = ifelse(area_fraction ≈ 0, zero(F_turb_ρτyz), F_turb_ρτyz)
    @. F_shf = ifelse(area_fraction ≈ 0, zero(F_shf), F_shf)
    @. F_lhf = ifelse(area_fraction ≈ 0, zero(F_lhf), F_lhf)
    @. F_turb_moisture = ifelse(area_fraction ≈ 0, zero(F_turb_moisture), F_turb_moisture)
    @. L_MO = ifelse(area_fraction ≈ 0, zero(L_MO), L_MO)
    @. ustar = ifelse(area_fraction ≈ 0, zero(ustar), ustar)
    @. buoyancy_flux = ifelse(area_fraction ≈ 0, zero(buoyancy_flux), buoyancy_flux)

    # multiply fluxes by area fraction
    F_turb_ρτxz .*= area_fraction
    F_turb_ρτyz .*= area_fraction
    F_shf .*= area_fraction
    F_lhf .*= area_fraction
    F_turb_moisture .*= area_fraction

    # perform additional diagnostics if required
    FluxCalculator.differentiate_turbulent_fluxes!(sim, (thermo_params, input_args, fluxes))

    # update the fluxes, which are now area-weighted, of this surface model
    fields = (;
        F_turb_ρτxz = F_turb_ρτxz,
        F_turb_ρτyz = F_turb_ρτyz,
        F_turb_energy = F_shf .+ F_lhf,
        F_turb_moisture = F_turb_moisture,
    )
    FluxCalculator.update_turbulent_fluxes!(sim, fields)

    # update fluxes in the coupler fields
    # add the flux contributing from this surface to the coupler field
    # note that the fluxes area-weighted, so if a surface model is
    #  not present at a point, the fluxes are zero
    @. csf.F_turb_ρτxz += F_turb_ρτxz
    @. csf.F_turb_ρτyz += F_turb_ρτyz
    @. csf.F_turb_energy += (F_shf .+ F_lhf)
    @. csf.F_turb_moisture += F_turb_moisture

    # NOTE: This is still an area weighted contribution, which maybe doesn't make
    # too much sense for these quantities...

    # L_MO can be Inf. We don't want to multiply Inf * 0, so we can handle this
    # separately.
    @. csf.L_MO += ifelse(isinf(L_MO), L_MO, L_MO * area_fraction)
    @. csf.ustar += ustar * area_fraction
    @. csf.buoyancy_flux += buoyancy_flux * area_fraction

    z0m = Interfacer.get_field(sim, Val(:roughness_momentum))
    z0b = Interfacer.get_field(sim, Val(:roughness_buoyancy))
    beta = Interfacer.get_field(sim, Val(:beta))

    @. csf.z0m_sfc += z0m * area_fraction
    @. csf.z0b_sfc += z0b * area_fraction
    @. csf.beta += beta * area_fraction

    # NOTE: This is essentially setting q_sfc to the Atmos q_sfc (because we compute the
    # thermo_state_sfc by extrapolating the atmos properties onto the surface)
    @. csf.q_sfc += TD.total_specific_humidity.(thermo_params, thermo_state_sfc) * area_fraction
    return nothing
end

end # module
