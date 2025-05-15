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

export extrapolate_ρ_to_sfc, turbulent_fluxes!, get_surface_params, update_turbulent_fluxes!, compute_surface_fluxes!

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
        :z0m_sfc,
        :z0b_sfc,
        :beta,
        :L_MO,
        :ustar,
        :buoyancy_flux,
    )
        fill!(getproperty(csf, p), 0)
    end

    # Compute the surface fluxes for each surface model and add them to `csf`
    for sim in model_sims
        compute_surface_fluxes!(csf, sim, atmos_sim, thermo_params)
    end

    # Update the atmosphere with the fluxes across all surface models
    # The surface models have already been updated with the fluxes in `compute_surface_fluxes!`
    # TODO this should be `update_turbulent_fluxes` to match the surface models
    Interfacer.update_field!(atmos_sim, Val(:turbulent_fluxes), csf)
    return nothing
end

function new_surface_inputs(input_args::NamedTuple)
    # Unpack stuff
    (; thermo_state_sfc, thermo_state_atmos, uₕ_int, z_int, z_sfc, scheme_properties, boundary_space) = input_args
    # Floattype
    FT = CC.Spaces.undertype(boundary_space)
    # Unpack roughness
    (; z0b, z0m, beta, gustiness) = scheme_properties
    z_int_fv = maybe_fv(z_int)
    uₕ_int_fv = maybe_fv(uₕ_int)
    thermo_state_atmos_fv = maybe_fv(thermo_state_atmos)
    z_sfc_fv = maybe_fv(z_sfc)
    thermo_state_sfc_fv = maybe_fv(thermo_state_sfc)
    beta_fv = maybe_fv(beta)
    z0m_fv = maybe_fv(z0m)
    z0b_fv = maybe_fv(z0b)
    gustiness_fv = maybe_fv(gustiness)
    return input_args
end


function surface_inputs(input_args::NamedTuple)
    (; thermo_state_sfc, thermo_state_atmos, uₕ_int, z_int, z_sfc, scheme_properties, boundary_space) = input_args
    FT = CC.Spaces.undertype(boundary_space)
    (; z0b, z0m, beta, gustiness) = scheme_properties

    # Extract the underlying data layouts of each field
    # Note: this is a bit "dangerous" because it circumvents ClimaCore, but
    #  it allows us to broadcast over fields on slightly different spaces
    maybe_fv = (x) -> x isa CC.Fields.Field ? CC.Fields.field_values(x) : x

    z_int_fv = maybe_fv(z_int)
    uₕ_int_fv = maybe_fv(uₕ_int)
    thermo_state_atmos_fv = maybe_fv(thermo_state_atmos)
    z_sfc_fv = maybe_fv(z_sfc)
    thermo_state_sfc_fv = maybe_fv(thermo_state_sfc)
    beta_fv = maybe_fv(beta)
    z0m_fv = maybe_fv(z0m)
    z0b_fv = maybe_fv(z0b)
    gustiness_fv = maybe_fv(gustiness)

    # Compute state values
    result = @. SF.ValuesOnly(
        SF.StateValues(z_int_fv, uₕ_int_fv, thermo_state_atmos_fv), # state_in
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
    get_surface_fluxes(inputs, surface_params::SF.Parameters.SurfaceFluxesParameters)

Uses SurfaceFluxes.jl to calculate turbulent surface fluxes. It should be atmos model agnostic, and columnwise.
Fluxes are computed over the entire surface, even where the relevant surface model is not present.

When available, it also computes ancillary quantities, such as the Monin-Obukov lengthscale.
"""
function get_surface_fluxes(inputs, input_args, surface_params::SF.Parameters.SurfaceFluxesParameters)
    # calculate all fluxes (saturated surface conditions)
    #outputs = SF.surface_conditions.(surface_params, inputs)

    ### TODO ######################################
    ### Attempt to use new surface flux formulation
    ts_sfc = input_args.thermo_state_sfc
    ts_atmos = input_args.thermo_state_atmos
    thermo_params = input_args.thermo_params
    u_int = input_args.u_int
    v_int = input_args.v_int
    z_int = input_args.z_int
    z_sfc = input_args.z_sfc
    ρ_sfc = input_args.ρ_sfc
    FT = eltype(u_int)
 
    ρ_atmos = TD.air_density.(thermo_params, ts_atmos)
    θ_sfc = TD.virtual_pottemp.(thermo_params, ts_sfc)
    θ_atmos = TD.virtual_pottemp.(thermo_params, ts_atmos)
    surface_params = input_args.surface_params
    
    # Testing with non-constant roughness length
    function z0test(surface_args, similarity_scales, atmos_state, surface_params)
        u★ = similarity_scales.momentum
        FT = eltype(u★)
        return FT.(0.015 .* u★ .^ 2 ./ 9.81)
    end

    # Property containers 
    atmos_state = SF.AtmosState(
                      u_int,
                      v_int,
                      TD.total_specific_humidity.(thermo_params, ts_atmos),
                      TD.virtual_pottemp.(thermo_params, ts_atmos),
                      z_int,
                      FT(1), # gustiness needs to be a function of u,v, ustar
                      FT(1000),
                      (ρ=ρ_atmos, arg𝑏=(0.01), arg𝑐=FT(0.001)),
                   )
 
    surface_state = SF.SurfaceState(
                      (𝑧0m=FT(0.01), 𝑧0θ=FT(0.01), 𝑧0q=z0test),
                      zero(u_int),
                      zero(v_int),
                      TD.total_specific_humidity.(thermo_params, ts_sfc),
                      TD.virtual_pottemp.(thermo_params, ts_sfc),
                      zero(z_int),
                      (arg𝑏=FT(0.01), arg𝑐=FT(0.001)),
                    )
 
    gustiness = atmos_state.gustiness_parameter
    (; 𝑧0m, 𝑧0θ, 𝑧0q) = surface_state.roughness_lengths

    # Initial guesses
    ζ₀ = fill!(similar(u_int), FT(10))
    var★ = fill!(similar(ζ₀), FT(1e-4))
    Σ₀ = SF.SimilarityScales(var★,var★,var★,var★)
    Δstate = SF.state_differences(surface_state, atmos_state, Σ₀, surface_params);
    L★ = Δstate.Δh ./ ζ₀
    Σ_est = Σ₀
    ΔU_est = fill!(similar(ζ₀), FT(10))

    # Similarity profiles
    UF = SF.UniversalFunctions
    param_set = SF.Parameters.SurfaceFluxesParameters(FT, UF.BusingerParams)
    similarity_theory = SF.Parameters.universal_func_type(param_set)
    sfc_params = SF.Parameters.uf_params(param_set)
    ufunc = SF.UniversalFunctions.universal_func.(Ref(similarity_theory), L★, Ref(sfc_params))
    
    outputs = SF.compute_similarity_theory_fluxes(ufunc, surface_state, atmos_state, param_set)
    ### TODO : End new surface flux formulation
    
    # drag
    F_turb_ρτxz = outputs.x_momentum
    F_turb_ρτyz = outputs.y_momentum

    # energy fluxes
    F_sh = outputs.sensible_heat
    F_lh = outputs.latent_heat

    # moisture
    #F_turb_moisture = SF.evaporation.(surface_params, inputs, outputs.Ch)
    F_turb_moisture = outputs.water_vapor

    # scale variables
    (u★, θ★, q★, b★) = outputs.scale_vars

    L_MO = outputs.l_mo
    ustar = u★
    buoyancy_flux = outputs.buoy_flux

    return (; F_turb_ρτxz, F_turb_ρτyz, F_sh, F_lh, F_turb_moisture, L_MO, ustar, buoyancy_flux)
end

"""
    get_surface_params(atmos_sim::Interfacer.AtmosModelSimulation)

Returns the surface parameters of type `SF.Parameters.SurfaceFluxesParameters`.

TODO: in the future this may not need to depend on the atmos sim, but
here retaining the dependency until we know how EDMF boundary conditions will
be handled (for consistency of parameters).
"""
function get_surface_params(atmos_sim::Interfacer.AtmosModelSimulation)
    return error("get_surface_params is required to be dispatched on $(nameof(atmos_sim)), but no method defined")
end

"""
    update_turbulent_fluxes!(sim::Interfacer.SurfaceModelSimulation, fields::NamedTuple)

Updates the fluxes in the surface model simulation `sim` with the fluxes in `fields`.
"""
function update_turbulent_fluxes!(sim::Interfacer.SurfaceModelSimulation, fields)
    return error("update_turbulent_fluxes! is required to be dispatched on $(nameof(sim)), but no method defined")
end

update_turbulent_fluxes!(sim::Interfacer.AbstractSurfaceStub, fields::NamedTuple) = nothing

"""
    compute_surface_fluxes!(csf, sim, atmos_sim, thermo_params)

This function computes surface fluxes between the input component model
simulation and the atmosphere.

Update the input coupler surface fields `csf` in-place with the computed fluxes
for this model. These are then summed using area-weighting across all surface
models to get the total fluxes.

Since the fluxes are computed between the input model and the atmosphere, this
function does nothing if called on an atmosphere model simulation.

# Arguments
- `csf`: [CC.Fields.Field] containing a NamedTuple of turbulent flux fields: `F_turb_ρτxz`, `F_turb_ρτyz`, `F_lh`, `F_sh`, `F_turb_moisture`.
- `sim`: [Interfacer.ComponentModelSimulation] the surface simulation to compute fluxes for.
- `atmos_sim`: [Interfacer.AtmosModelSimulation] the atmosphere simulation to compute fluxes with.
- `thermo_params`: [TD.Parameters.ThermodynamicsParameters] the thermodynamic parameters.
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
    sim::Interfacer.SurfaceModelSimulation,
    atmos_sim::Interfacer.AtmosModelSimulation,
    thermo_params,
)
    boundary_space = axes(csf)
    # TODO lots of allocations here
    # `_int` refers to atmos state of center level 1
    z_int = Interfacer.get_field(atmos_sim, Val(:height_int), boundary_space)
    u_int = Interfacer.get_field(atmos_sim, Val(:u_int), boundary_space)
    v_int = Interfacer.get_field(atmos_sim, Val(:v_int), boundary_space)
    uₕ_int = @. StaticArrays.SVector(u_int, v_int)
    z_sfc = Interfacer.get_field(atmos_sim, Val(:height_sfc), boundary_space)

    # construct the atmospheric thermo states
    thermo_state_atmos = TD.PhaseEquil_ρTq.(thermo_params, csf.ρ_atmos, csf.T_atmos, csf.q_atmos)

    # construct the surface thermo state
    # get surface air density by extrapolating atmospheric density to the surface
    ρ_sfc = extrapolate_ρ_to_sfc.(thermo_params, thermo_state_atmos, csf.T_sfc)

    # compute surface humidity from the surface temperature, surface density, and phase
    thermo_state_sfc = TD.PhaseEquil_ρTq.(thermo_params, ρ_sfc, csf.T_sfc, csf.q_sfc)

    # get area fraction (min = 0, max = 1)
    area_fraction = Interfacer.get_field(sim, Val(:area_fraction))

    surface_params = FluxCalculator.get_surface_params(atmos_sim)

    z0m = Interfacer.get_field(sim, Val(:roughness_momentum), boundary_space)
    z0b = Interfacer.get_field(sim, Val(:roughness_buoyancy), boundary_space)
    beta = Interfacer.get_field(sim, Val(:beta), boundary_space)
    FT = eltype(z0m)
    scheme_properties = (; z0b = z0b, z0m = z0m, Ch = FT(0), Cd = FT(0), beta = beta, gustiness = FT(1))

    input_args = (;
        thermo_state_sfc,
        thermo_state_atmos,
        uₕ_int,
        z_int,
        z_sfc,
        scheme_properties,
        boundary_space,
        surface_params,
        thermo_params,
        u_int,
        v_int,
        ρ_sfc,
    )
    inputs = FluxCalculator.surface_inputs(input_args)

    # calculate the surface fluxes
    fluxes = FluxCalculator.get_surface_fluxes(inputs, input_args, surface_params)
    (; F_turb_ρτxz, F_turb_ρτyz, F_sh, F_lh, F_turb_moisture, L_MO, ustar, buoyancy_flux) = fluxes

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

    # multiply fluxes by area fraction
    F_turb_ρτxz .*= area_fraction
    F_turb_ρτyz .*= area_fraction
    F_sh .*= area_fraction
    F_lh .*= area_fraction
    F_turb_moisture .*= area_fraction

    # update the fluxes, which are now area-weighted, of this surface model
    fields = (; F_turb_ρτxz, F_turb_ρτyz, F_lh, F_sh, F_turb_moisture)
    FluxCalculator.update_turbulent_fluxes!(sim, fields)

    # update fluxes in the coupler fields
    # add the flux contributing from this surface to the coupler field
    # note that the fluxes area-weighted, so if a surface model is
    #  not present at a point, the fluxes are zero
    @. csf.F_turb_ρτxz += F_turb_ρτxz
    @. csf.F_turb_ρτyz += F_turb_ρτyz
    @. csf.F_lh += F_lh
    @. csf.F_sh += F_sh
    @. csf.F_turb_moisture += F_turb_moisture

    # NOTE: This is still an area weighted contribution, which maybe doesn't make
    # too much sense for these quantities...

    # L_MO can be Inf. We don't want to multiply Inf * 0, so we can handle this
    # separately.
    @. csf.L_MO += ifelse(isinf(L_MO), L_MO, L_MO * area_fraction)
    @. csf.ustar += ustar * area_fraction
    @. csf.buoyancy_flux += buoyancy_flux * area_fraction

    z0m = Interfacer.get_field(sim, Val(:roughness_momentum), boundary_space)
    z0b = Interfacer.get_field(sim, Val(:roughness_buoyancy), boundary_space)
    beta = Interfacer.get_field(sim, Val(:beta), boundary_space)

    @. csf.z0m_sfc += z0m * area_fraction
    @. csf.z0b_sfc += z0b * area_fraction
    @. csf.beta += beta * area_fraction
    return nothing
end

end # module
