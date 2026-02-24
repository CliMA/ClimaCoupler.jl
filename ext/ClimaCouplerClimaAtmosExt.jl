"""
    ClimaCouplerClimaAtmosExt

This module contains code for extending the ClimaCoupler interface for
ClimaAtmos. For more information about the atmos model, please see
`experiments/ClimaEarth/README.md`
"""
module ClimaCouplerClimaAtmosExt

# atmos_init: for ClimaAtmos interface
import StaticArrays
import ClimaAtmos as CA
import ClimaAtmos.Parameters as CAP
import ClimaParams as CP
import ClimaCore as CC
import ClimaCore.Geometry: ⊗
import SurfaceFluxes as SF
import Thermodynamics as TD
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities, Plotting
import ClimaUtilities.TimeManager: ITime
import ClimaComms

###
### Functions required by ClimaCoupler.jl for an AbstractAtmosSimulation
###
struct ClimaAtmosSimulation{P, D, I, OW} <: Interfacer.AbstractAtmosSimulation
    params::P
    domain::D
    integrator::I
    output_writers::OW
end

function hasradiation(integrator)
    return !isnothing(integrator.p.atmos.radiation_mode)
end

function hasmoisture(integrator)
    return !(integrator.p.atmos.microphysics_model isa CA.DryModel)
end

"""
    Interfacer.AtmosSimulation(::Val{:climaatmos}; kwargs...)

Extension of the generic AtmosSimulation constructor for ClimaAtmos.

Note that this is currently the only atmosphere model supported by
ClimaCoupler.jl.
"""
function Interfacer.AtmosSimulation(
    ::Val{:climaatmos};
    config_dict,
    atmos_output_dir,
    coupled_param_dict,
    comms_ctx,
)
    return ClimaAtmosSimulation(
        config_dict,
        atmos_output_dir,
        coupled_param_dict,
        comms_ctx,
    )
end

function ClimaAtmosSimulation(
    config_dict::Dict,
    atmos_output_dir,
    coupled_param_dict::CP.ParamDict,
    comms_ctx,
)
    atmos_config =
        get_atmos_config_dict(config_dict, atmos_output_dir, coupled_param_dict, comms_ctx)
    return ClimaAtmosSimulation(atmos_config)
end

function ClimaAtmosSimulation(atmos_config)
    # By passing `parsed_args` to `AtmosConfig`, `parsed_args` overwrites the default atmos config
    FT = atmos_config.parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
    simulation = CA.get_simulation(atmos_config)
    (; integrator, output_writers) = simulation
    Y = integrator.u
    center_space = axes(Y.c.ρe_tot)
    face_space = axes(Y.f.u₃)
    spaces = (; center_space = center_space, face_space = face_space)

    # define shorter references for long variable names to increase readability, and set to zero
    ρ_flux_h_tot = integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot
    ρ_flux_uₕ = integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ
    @. ρ_flux_h_tot = CC.Geometry.Covariant3Vector(FT(0.0))
    ρ_flux_uₕ.components .= Ref(StaticArrays.SMatrix{1, 2}([FT(0), FT(0)]))

    if hasmoisture(integrator)
        ρ_flux_q_tot = integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot
        @. ρ_flux_q_tot = CC.Geometry.Covariant3Vector(FT(0.0))

        integrator.p.precomputed.surface_rain_flux .= FT(0)
        integrator.p.precomputed.surface_snow_flux .= FT(0)
    end

    microphysics_model = integrator.p.atmos.microphysics_model
    if microphysics_model isa CA.EquilibriumMicrophysics0M
        ᶜS_ρq_tot = integrator.p.precomputed.ᶜS_ρq_tot
        ᶜS_ρq_tot .= FT(0)
    end
    if hasradiation(integrator)
        ᶠradiation_flux = integrator.p.radiation.ᶠradiation_flux
        @. ᶠradiation_flux = CC.Geometry.WVector(FT(0))
    end

    sim = ClimaAtmosSimulation(integrator.p.params, spaces, integrator, output_writers)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

"""
    get_surface_space(sim::ClimaAtmosSimulation)

Get the surface space of the atmosphere simulation, which is the bottom
level of the face space.
"""
get_surface_space(sim::ClimaAtmosSimulation) = axes(
    CC.Spaces.level(CC.Fields.coordinate_field(sim.domain.face_space).z, CC.Utilities.half),
)

"""
    Checkpointer.get_model_prog_state(sim::ClimaAtmosSimulation)

Extension of Checkpointer.get_model_prog_state to get the model state.
"""
function Checkpointer.get_model_prog_state(sim::ClimaAtmosSimulation)
    return sim.integrator.u
end

function Checkpointer.get_model_cache(sim::ClimaAtmosSimulation)
    return sim.integrator.p
end

"""
    Checkpointer.restore_cache!(sim::ClimaAtmosSimulation, cache_vec)

Restore the cache in `sim` using `cache_vec` saved from
`checkpoint_model_cache`.

The atmosphere cache is saved as a vector of objects instead of being saved as
the entire cache, as is done for the other models. This reduces the file size of
the saved cache. See `checkpoint_model_cache` and
`get_model_cache_to_checkpoint` for more information on how the atmosphere cache
is saved.
"""
function Checkpointer.restore_cache!(sim::ClimaAtmosSimulation, cache_vec)
    comms_ctx = ClimaComms.context(sim.integrator.u.c)
    atmos_cache_itr = Checkpointer.CacheIterator(sim)
    # This assumes that the traversals over the two atmos caches are the same
    for (i, obj_restored) in enumerate(atmos_cache_itr)
        obj_saved = cache_vec[i]
        Checkpointer.restore!(obj_restored, obj_saved, comms_ctx; name = "", ignore = Set())
    end
    return nothing
end


"""
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux_toa})

Extension of Interfacer.get_field to get the net TOA radiation, which is a sum of the
upward and downward longwave and shortwave radiation.
"""
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux_toa})
    FT = eltype(sim.integrator.u)

    if hasradiation(sim.integrator)
        face_space = axes(sim.integrator.u.f)
        nz_faces = length(CC.Spaces.vertical_topology(face_space).mesh.faces)

        (; face_lw_flux_dn, face_lw_flux_up, face_sw_flux_dn, face_sw_flux_up) =
            sim.integrator.p.radiation.rrtmgp_model

        LWd_TOA = CC.Fields.level(
            CC.Fields.array2field(FT.(face_lw_flux_dn), face_space),
            nz_faces - CC.Utilities.half,
        )
        LWu_TOA = CC.Fields.level(
            CC.Fields.array2field(FT.(face_lw_flux_up), face_space),
            nz_faces - CC.Utilities.half,
        )
        SWd_TOA = CC.Fields.level(
            CC.Fields.array2field(FT.(face_sw_flux_dn), face_space),
            nz_faces - CC.Utilities.half,
        )
        SWu_TOA = CC.Fields.level(
            CC.Fields.array2field(FT.(face_sw_flux_up), face_space),
            nz_faces - CC.Utilities.half,
        )

        return @. -(LWd_TOA + SWd_TOA - LWu_TOA - SWu_TOA)
    else
        return FT[0]
    end
end

function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:energy})
    integrator = sim.integrator
    p = integrator.p

    # return total energy and (if EquilibriumMicrophysics0M) the energy lost due to precipitation removal
    microphysics_model = integrator.p.atmos.microphysics_model
    if microphysics_model isa CA.EquilibriumMicrophysics0M
        (; ᶜT, ᶜq_liq_rai, ᶜq_ice_sno, ᶜS_ρq_tot) = p.precomputed
        (; ᶜΦ) = p.core
        thermo_params = get_thermo_params(sim)
        return integrator.u.c.ρe_tot .-
               ᶜS_ρq_tot .*
               CA.e_tot_0M_precipitation_sources_helper.(
            Ref(thermo_params),
            ᶜT,
            ᶜq_liq_rai,
            ᶜq_ice_sno,
            ᶜΦ,
        ) .* float(integrator.dt)
    else
        return integrator.u.c.ρe_tot
    end
end

# helpers for get_field extensions, dipatchable on different moisture model options and radiation modes

surface_rain_flux(::CA.DryModel, integrator) = eltype(integrator.u)(0)
function surface_rain_flux(
    ::Union{
        CA.EquilibriumMicrophysics0M,
        CA.NonEquilibriumMicrophysics1M,
        CA.NonEquilibriumMicrophysics2M,
        CA.NonEquilibriumMicrophysics2MP3,
    },
    integrator,
)
    return integrator.p.precomputed.surface_rain_flux
end

surface_snow_flux(::CA.DryModel, integrator) = eltype(integrator.u)(0)
function surface_snow_flux(
    ::Union{
        CA.EquilibriumMicrophysics0M,
        CA.NonEquilibriumMicrophysics1M,
        CA.NonEquilibriumMicrophysics2M,
        CA.NonEquilibriumMicrophysics2MP3,
    },
    integrator,
)
    return integrator.p.precomputed.surface_snow_flux
end

surface_radiation_flux(::Nothing, integrator) = eltype(integrator.u)(0)
surface_radiation_flux(::CA.RRTMGPI.AbstractRRTMGPMode, integrator) =
    CC.Fields.level(integrator.p.radiation.ᶠradiation_flux, CC.Utilities.half)

moisture_flux(::CA.DryModel, integrator) = eltype(integrator.u)(0)
moisture_flux(
    ::Union{
        CA.EquilibriumMicrophysics0M,
        CA.NonEquilibriumMicrophysics1M,
        CA.NonEquilibriumMicrophysics2M,
        CA.NonEquilibriumMicrophysics2MP3,
    },
    integrator,
) = CC.Geometry.WVector.(integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot)

ρq_tot(::CA.DryModel, integrator) = eltype(integrator.u)(0)
ρq_tot(
    ::Union{
        CA.EquilibriumMicrophysics0M,
        CA.NonEquilibriumMicrophysics1M,
        CA.NonEquilibriumMicrophysics2M,
        CA.NonEquilibriumMicrophysics2MP3,
    },
    integrator,
) = integrator.u.c.ρq_tot

# extensions required by the Interfacer
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:air_pressure}) =
    CC.Fields.level(sim.integrator.p.precomputed.ᶜp, 1)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:air_temperature}) =
    CC.Fields.level(sim.integrator.p.precomputed.ᶜT, 1)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:air_density}) =
    CC.Fields.level(sim.integrator.u.c.ρ, 1)
# When CO2 is stored as a 1D Array in the tracers cache, we access
# it from there. Otherwise, we pull a fixed value from ClimaParams.
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:co2}) =
    :co2 in propertynames(sim.integrator.p.tracers) ? sim.integrator.p.tracers.co2[1] :
    CAP.trace_gas_params(sim.integrator.p.params).CO2_fixed_value

function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:diffuse_fraction})
    # Diffuse fraction doesn't matter when we don't have radiation, so return zero
    FT = eltype(sim.integrator.u)
    hasradiation(sim.integrator) || return zero(FT)

    radiation_model = sim.integrator.p.radiation.rrtmgp_model
    # only take the first level
    total_flux_dn = radiation_model.face_sw_flux_dn[1, :]
    lowest_face_space = CC.Spaces.level(axes(sim.integrator.u.f), CC.Utilities.half)
    if radiation_model.radiation_mode isa CA.RRTMGPInterface.GrayRadiation
        diffuse_fraction = zero(total_flux_dn)
    else
        direct_flux_dn = radiation_model.face_sw_direct_flux_dn[1, :]
        FT = eltype(total_flux_dn)
        diffuse_fraction =
            clamp.(
                (
                    (x, y) -> y > zero(y) ? x / y : zero(y)
                ).(total_flux_dn .- direct_flux_dn, total_flux_dn),
                zero(FT),
                one(FT),
            )
    end
    return CC.Fields.array2field(diffuse_fraction, lowest_face_space)
end
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:liquid_precipitation}) =
    surface_rain_flux(sim.integrator.p.atmos.microphysics_model, sim.integrator)
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:LW_d})
    # If we don't have radiation, downwelling LW is zero
    FT = eltype(sim.integrator.u)
    hasradiation(sim.integrator) || return zero(FT)

    return CC.Fields.level(
        CC.Fields.array2field(
            sim.integrator.p.radiation.rrtmgp_model.face_lw_flux_dn,
            axes(sim.integrator.u.f),
        ),
        CC.Utilities.half,
    )
end
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:total_specific_humidity}) =
    CC.Fields.level(sim.integrator.p.precomputed.ᶜq_tot_safe, 1)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:liquid_specific_humidity}) =
    CC.Fields.level(sim.integrator.p.precomputed.ᶜq_liq_rai, 1)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:ice_specific_humidity}) =
    CC.Fields.level(sim.integrator.p.precomputed.ᶜq_ice_sno, 1)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:snow_precipitation}) =
    surface_snow_flux(sim.integrator.p.atmos.microphysics_model, sim.integrator)
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:SW_d})
    # If we don't have radiation, downwelling SW is zero
    FT = eltype(sim.integrator.u)
    hasradiation(sim.integrator) || return zero(FT)

    return CC.Fields.level(
        CC.Fields.array2field(
            sim.integrator.p.radiation.rrtmgp_model.face_sw_flux_dn,
            axes(sim.integrator.u.f),
        ),
        CC.Utilities.half,
    )
end
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:water}) =
    ρq_tot(sim.integrator.p.atmos.microphysics_model, sim.integrator)

function Interfacer.update_field!(
    sim::ClimaAtmosSimulation,
    ::Val{:surface_temperature},
    field,
)
    Interfacer.remap!(sim.integrator.p.precomputed.sfc_conditions.T_sfc, field)
end
function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_humidity}, csf)
    # NOTE: This update_field! takes as argument the entire coupler fields struct, instead
    # of a single field. This is unlike most of other functions, so we may want to revisit it.

    # Remap coupler fields onto the atmosphere surface space
    atmos_surface_space = get_surface_space(sim)
    temp_field_surface = sim.integrator.p.scratch.ᶠtemp_field_level
    @assert axes(temp_field_surface) == atmos_surface_space

    # Compute surface humidity on the coupler space, then remap to atmosphere surface space
    Interfacer.get_field!(csf.scalar_temp1, sim, Val(:air_temperature))
    T_atmos = csf.scalar_temp1
    Interfacer.get_field!(csf.scalar_temp2, sim, Val(:total_specific_humidity))
    q_tot_atmos = csf.scalar_temp2
    Interfacer.get_field!(csf.scalar_temp3, sim, Val(:air_density))
    ρ_atmos = csf.scalar_temp3

    surface_fluxes_params = FluxCalculator.get_surface_params(sim)

    # TODO: Is csf.height_int .- csf.height_sfc allocating?
    # TODO: Add two scratch fields for q_liq_atmos and q_ice_atmos
    csf.scalar_temp4 .=
        SF.surface_density.(
            surface_fluxes_params,
            T_atmos,
            ρ_atmos,
            csf.T_sfc,
            csf.height_int .- csf.height_sfc,
            q_tot_atmos,
            0, # q_liq
            0, # q_ice
        )
    ρ_sfc = csf.scalar_temp4

    thermo_params = get_thermo_params(sim)
    csf.scalar_temp1 .= TD.q_vap_saturation.(thermo_params, csf.T_sfc, ρ_sfc, 0, 0)

    # Remap surface temperature and humidity to atmosphere surface space
    q_sfc_atmos = Interfacer.remap(atmos_surface_space, csf.scalar_temp1)
    sim.integrator.p.precomputed.sfc_conditions.q_vap_sfc .= q_sfc_atmos
end
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:height_int}) =
    CC.Spaces.level(CC.Fields.coordinate_field(sim.integrator.u.c).z, 1)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:height_sfc}) =
    CC.Spaces.level(CC.Fields.coordinate_field(sim.integrator.u.f).z, CC.Utilities.half)
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:u_int})
    # NOTE: This calculation is copied from ClimaAtmos (and is allocating! Fix me if you can!)
    int_local_geometry_values =
        CC.Fields.level(CC.Fields.local_geometry_field(sim.integrator.u.c), 1)
    int_u_values = CC.Spaces.level(sim.integrator.p.precomputed.ᶜu, 1)
    return CA.projected_vector_data.(CA.CT1, int_u_values, int_local_geometry_values)
end
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:v_int})
    # NOTE: This calculation is copied from ClimaAtmos (and is allocating! Fix me if you can!)
    int_local_geometry_values =
        CC.Fields.level(CC.Fields.local_geometry_field(sim.integrator.u.c), 1)
    int_u_values = CC.Spaces.level(sim.integrator.p.precomputed.ᶜu, 1)
    return CA.projected_vector_data.(CA.CT2, int_u_values, int_local_geometry_values)
end

# extensions required by the Interfacer
function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:emissivity}, field)
    if hasradiation(sim.integrator)
        # Remap field onto the atmosphere surface space in scratch field
        temp_field_surface = sim.integrator.p.scratch.ᶠtemp_field_level
        Interfacer.remap!(temp_field_surface, field)
        # Set each row (band) of the emissivity matrix by transposing the vector returned from `field2array`
        sim.integrator.p.radiation.rrtmgp_model.surface_emissivity .=
            CC.Fields.field2array(temp_field_surface)'
    end
    return nothing
end
function Interfacer.update_field!(
    sim::ClimaAtmosSimulation,
    ::Val{:surface_direct_albedo},
    field,
)
    if hasradiation(sim.integrator)
        # Remap field onto the atmosphere surface space in scratch field
        temp_field_surface = sim.integrator.p.scratch.ᶠtemp_field_level
        Interfacer.remap!(temp_field_surface, field)
        # Set each row (band) of the albedo matrix by transposing the vector returned from `field2array`
        sim.integrator.p.radiation.rrtmgp_model.direct_sw_surface_albedo .=
            CC.Fields.field2array(temp_field_surface)'
    end
    return nothing
end

function Interfacer.update_field!(
    sim::ClimaAtmosSimulation,
    ::Val{:surface_diffuse_albedo},
    field,
)
    if hasradiation(sim.integrator)
        # Remap field onto the atmosphere surface space in scratch field
        temp_field_surface = sim.integrator.p.scratch.ᶠtemp_field_level
        Interfacer.remap!(temp_field_surface, field)
        # Set each row (band) of the albedo matrix by transposing the vector returned from `field2array`
        sim.integrator.p.radiation.rrtmgp_model.diffuse_sw_surface_albedo .=
            CC.Fields.field2array(temp_field_surface)'
    end
    return nothing
end

function FluxCalculator.update_turbulent_fluxes!(sim::ClimaAtmosSimulation, fields)
    (; F_lh, F_sh, F_turb_moisture, F_turb_ρτxz, F_turb_ρτyz) = fields
    atmos_surface_space = get_surface_space(sim)
    temp_field_surface = sim.integrator.p.scratch.ᶠtemp_field_level

    Y = sim.integrator.u
    surface_local_geometry =
        CC.Fields.level(CC.Fields.local_geometry_field(Y.f), CC.Utilities.half)
    surface_normal = @. CA.C3(CA.unit_basis_vector_data(CA.C3, surface_local_geometry))

    # get template objects for the contravariant components of the momentum fluxes (required by Atmos boundary conditions)
    # NOTE: This is allocating! Fix me if you can!
    vec_ct12_ct1 = @. CA.CT12(
        CA.CT1(CA.unit_basis_vector_data(CA.CT1, surface_local_geometry)),
        surface_local_geometry,
    )
    vec_ct12_ct2 = @. CA.CT12(
        CA.CT2(CA.unit_basis_vector_data(CA.CT2, surface_local_geometry)),
        surface_local_geometry,
    )

    # NOTE: This is allocating (F_turb_ρτyz_atmos)! If we had 1 more scratch field, we could avoid this
    Interfacer.remap!(temp_field_surface, F_turb_ρτxz) # F_turb_ρτxz_atmos
    F_turb_ρτyz_atmos = Interfacer.remap(atmos_surface_space, F_turb_ρτyz) # F_turb_ρτyz_atmos
    sim.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ .= (
        surface_normal .⊗
        CA.C12.(
            temp_field_surface .* vec_ct12_ct1 .+ F_turb_ρτyz_atmos .* vec_ct12_ct2,
            surface_local_geometry,
        )
    )

    # Remap summed turbulent energy fluxes into the scratch field
    Interfacer.remap!(temp_field_surface, F_lh .+ F_sh)
    sim.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot .=
        temp_field_surface .* surface_normal # (shf + lhf)

    # Remap turbulent moisture flux into the scratch field
    Interfacer.remap!(temp_field_surface, F_turb_moisture)
    sim.integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot .=
        temp_field_surface .* surface_normal # (evap)

    Interfacer.remap!(
        sim.integrator.p.precomputed.sfc_conditions.obukhov_length,
        fields.L_MO,
    )
    Interfacer.remap!(sim.integrator.p.precomputed.sfc_conditions.ustar, fields.ustar)
    Interfacer.remap!(
        sim.integrator.p.precomputed.sfc_conditions.buoyancy_flux,
        fields.buoyancy_flux,
    )

    return nothing
end

# extensions required by FieldExchanger
Interfacer.step!(sim::ClimaAtmosSimulation, t::Real) =
    Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
function Interfacer.step!(sim::ClimaAtmosSimulation, t::ITime)
    while sim.integrator.t < t
        Interfacer.step!(sim.integrator)
    end
    return nothing
end

"""
Extend Interfacer.add_coupler_fields! to add the fields required for ClimaAtmosSimulation.

The fields added are:
- `:surface_direct_albedo` (for radiation)
- `:surface_diffuse_albedo` (for radiation)
- `:ustar`, `:L_MO`, `:buoyancy_flux` (for EDMF boundary conditions)
"""
function Interfacer.add_coupler_fields!(
    coupler_field_names,
    atmos_sim::ClimaAtmosSimulation,
)
    atmos_coupler_fields =
        [:surface_direct_albedo, :surface_diffuse_albedo, :ustar, :L_MO, :buoyancy_flux]
    push!(coupler_field_names, atmos_coupler_fields...)
end

"""
    Interfacer.close_output_writers(sim::ClimaAtmosSimulation)

Close all output writers used by the atmos simulation.
"""
Interfacer.close_output_writers(sim::ClimaAtmosSimulation) =
    isnothing(sim.output_writers) || foreach(close, sim.output_writers)


# FluxCalculator.get_surface_params required by FluxCalculator (partitioned fluxes)
FluxCalculator.get_surface_params(sim::ClimaAtmosSimulation) =
    CAP.surface_fluxes_params(sim.integrator.p.params)

"""
    Interfacer.set_cache!(sim::ClimaAtmosSimulation, csf)

Set cache variables that cannot be initialized before the initial exchange.
Radiation must be called here because it requires the surface model initial
temperatures, which are set in the atmosphere during the initial exchange.
Any other callbacks which modify the cache should be called here as well.

This function does not set all the cache variables, because many are computed
as part of the tendendencies.
"""
function Interfacer.set_cache!(sim::ClimaAtmosSimulation, csf)
    if hasradiation(sim.integrator)
        CA.rrtmgp_model_callback!(sim.integrator)
        CA.nogw_model_callback!(sim.integrator)
    end
    return nothing
end

###
### ClimaAtmos.jl model-specific functions (not explicitly required by ClimaCoupler.jl)
###
"""
    get_atmos_config_dict(coupler_config::Dict, atmos_output_dir)

Returns the specified atmospheric configuration (`atmos_config`) overwitten by arguments
in the coupler dictionary (`coupler_config`).
The returned `atmos_config` dictionary will then be used to set up the atmosphere simulation.

The `atmos_config` is mostly a copy of the `coupler_config`, which has already been created
by combining the default and provided atmosphere and coupler configuration files.
Please see the documentation of `get_coupler_config_dict` for more details about how that is done.

This function makes the following changes to `coupler_config` when creating `atmos_config`:
- Renaming for consistency with ClimaAtmos.jl
- TOML parameter file selection (see details below)
- Setting the atmos output directory to the provided `atmos_output_dir`
- Adding extra atmosphere diagnostics
- Manually setting the restart file suffix

The TOML parameter file to use is chosen using the following priority:
If a coupler TOML file is provided, it is used.
Otherwise we use an atmos TOML file if it's provided.
If neither file is provided, the default ClimaAtmos parameters are used.
This is analogous to the logic in `get_coupler_config_dict` for selecting
the coupler configuration TOML file.
"""
function get_atmos_config_dict(
    coupler_config::Dict,
    atmos_output_dir,
    coupled_param_dict::CP.ParamDict,
    comms_ctx,
)
    atmos_config = deepcopy(coupler_config)

    # Ensure Atmos's own checkpoints are synced up with ClimaCoupler, so that we
    # can pick up from where we have left. NOTE: This should not be needed, but
    # there is no easy way to initialize ClimaAtmos with a different t_start
    atmos_config["dt_save_state_to_disk"] = coupler_config["checkpoint_dt"]
    atmos_config["log_progress"] = coupler_config["atmos_log_progress"]

    # The Atmos `get_simulation` function expects the atmos config to contains its timestep size
    # in the `dt` field. If there is a `dt_atmos` field in coupler_config, we add it to the atmos config as `dt`
    dt_atmos =
        haskey(coupler_config, "dt_atmos") ? coupler_config["dt_atmos"] :
        coupler_config["dt"]
    atmos_config["dt"] = dt_atmos

    # Set up the atmosphere parameter dictionary (`toml_dict`)
    # If we can't find the atmos TOML file at the relative path, prepend pkgdir(ClimaAtmos)
    atmos_toml = map(atmos_config["toml"]) do file
        isfile(file) ? file : joinpath(pkgdir(CA), file)
    end
    # Use atmos toml only if coupler toml is not defined
    coupler_toml = atmos_config["coupler_toml"]
    toml = isempty(coupler_toml) ? atmos_toml : coupler_toml
    if !isempty(toml)
        @info "Overwriting Atmos parameters from input TOML file(s): $toml"
        atmos_config["toml"] = toml
    end

    # Override atmos parameters with coupled parameters
    FT = atmos_config["FLOAT_TYPE"] == "Float64" ? Float64 : Float32
    override_file = CP.merge_toml_files(atmos_config["toml"], override = true)
    atmos_toml = CP.create_toml_dict(FT; override_file)
    toml_dict = CP.merge_override_default_values(coupled_param_dict, atmos_toml)

    # Specify atmos output directory to be inside the coupler output directory
    atmos_config["output_dir_style"] = "RemovePreexisting"
    atmos_config["output_dir"] = atmos_output_dir

    # Add all extra atmos diagnostic entries into the vector of atmos diagnostics
    atmos_config["diagnostics"] =
        haskey(atmos_config, "diagnostics") ?
        vcat(atmos_config["diagnostics"], coupler_config["extra_atmos_diagnostics"]) :
        coupler_config["extra_atmos_diagnostics"]

    # set restart file to the initial file saved in this location if it is not nothing
    # TODO this is hardcoded and should be fixed once we have a better restart system
    if !isnothing(atmos_config["restart_file"])
        atmos_config["restart_file"] =
            replace(atmos_config["restart_file"], "active" => "0000")
    end

    # Construct the AtmosConfig object
    config_files = atmos_config["toml"]
    job_id = atmos_config["job_id"]
    config = CA.override_default_config(atmos_config)

    TD = typeof(toml_dict)
    PA = typeof(config)
    C = typeof(comms_ctx)
    CF = typeof(config_files)
    return CA.AtmosConfig{FT, TD, PA, C, CF}(
        toml_dict,
        config,
        comms_ctx,
        config_files,
        job_id,
    )
end

"""
    get_thermo_params(sim::ClimaAtmosSimulation)

Returns the thermodynamic parameters from the atmospheric model simulation object.
"""
get_thermo_params(sim::ClimaAtmosSimulation) =
    CAP.thermodynamics_params(sim.integrator.p.params)

"""
    dss_state!(sim::ClimaAtmosSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on atmosphere simulations.
These sims don't store a dss buffer in their cache, so we must allocate
one here.
"""
function dss_state!(sim::ClimaAtmosSimulation)
    Y = sim.integrator.u
    for key in propertynames(Y)
        field = getproperty(Y, key)
        buffer = CC.Spaces.create_dss_buffer(field)
        CC.Spaces.weighted_dss!(field, buffer)
    end
end

"""
    climaatmos_restart_path(output_dir_root, t)

Look at the most recent output in `output_dir_root` and find a checkpoint for time `t`.
"""
function climaatmos_restart_path(output_dir_root, t)
    isdir(output_dir_root) || error("$(output_dir_root) does not exist")
    name_rx = r"output_(\d\d\d\d)"
    existing_outputs = filter(x -> !isnothing(match(name_rx, x)), readdir(output_dir_root))

    day = floor(Int, t / (60 * 60 * 24))
    sec = floor(Int, t % (60 * 60 * 24))

    # Walk back the folders and tyr to find a checkpoint
    for output in sort(existing_outputs, rev = true)
        previous_folder = joinpath(output_dir_root, output)
        restart_file = joinpath(previous_folder, "clima_atmos", "day$day.$sec.hdf5")
        ispath(restart_file) && return restart_file
    end
    error("Restart file for time $t not found")
end

# Additional ClimaAtmos getter methods for plotting debug fields
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:w})
    w_c = ones(CC.Spaces.horizontal_space(sim.domain.face_space))
    parent(w_c) .= parent(
        CC.Fields.level(
            CC.Geometry.WVector.(sim.integrator.u.f.u₃),
            5 .+ CC.Utilities.half,
        ),
    )
    return w_c
end
specific_humidity(::CA.DryModel, integrator) = [eltype(integrator.u)(0)]
specific_humidity(
    ::Union{
        CA.EquilibriumMicrophysics0M,
        CA.NonEquilibriumMicrophysics1M,
        CA.NonEquilibriumMicrophysics2M,
        CA.NonEquilibriumMicrophysics2MP3,
    },
    integrator,
) = integrator.u.c.ρq_tot
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:ρq_tot}) =
    specific_humidity(sim.integrator.p.atmos.microphysics_model, sim.integrator)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:ρe_tot}) = sim.integrator.u.c.ρe_tot

Plotting.debug_plot_fields(sim::ClimaAtmosSimulation) =
    (:w, :ρq_tot, :ρe_tot, :liquid_precipitation, :snow_precipitation)


end
