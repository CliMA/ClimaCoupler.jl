# atmos_init: for ClimaAtmos interface
import StaticArrays
import Statistics
import LinearAlgebra
import ClimaAtmos as CA
import ClimaAtmos: set_surface_albedo!
import ClimaAtmos.Parameters as CAP
import ClimaCore as CC
import ClimaCore.Geometry: ⊗
import SurfaceFluxes as SF
import Thermodynamics as TD
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaUtilities.TimeManager: ITime

include("climaatmos_extra_diags.jl")

if pkgversion(CA) < v"0.28.6"
    # Allow cache to be moved to CPU (this is a little bit of type piracy, but we
    # allow it in this particular file)
    CC.Adapt.@adapt_structure CA.AtmosCache
    CC.Adapt.@adapt_structure CA.RRTMGPInterface.RRTMGPModel
end

include("../shared/restore.jl")

###
### Functions required by ClimaCoupler.jl for an AtmosModelSimulation
###
struct ClimaAtmosSimulation{P, D, I, OW} <: Interfacer.AtmosModelSimulation
    params::P
    domain::D
    integrator::I
    output_writers::OW
end

function hasradiation(integrator)
    return !isnothing(integrator.p.atmos.radiation_mode)
end

function hasmoisture(integrator)
    return !(integrator.p.atmos.moisture_model isa CA.DryModel)
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
        if pkgversion(CA) >= v"0.29.0"
            surface_rain_flux = integrator.p.precomputed.surface_rain_flux
            surface_snow_flux = integrator.p.precomputed.surface_snow_flux
        else
            surface_rain_flux = integrator.p.precipitation.surface_rain_flux
            surface_snow_flux = integrator.p.precipitation.surface_snow_flux
        end
        surface_rain_flux .= FT(0)
        surface_snow_flux .= FT(0)
    end

    microphysics_model =
        pkgversion(CA) < v"0.31.0" ? integrator.p.atmos.precip_model : integrator.p.atmos.microphysics_model
    if microphysics_model isa CA.Microphysics0Moment
        if pkgversion(CA) >= v"0.29.0"
            ᶜS_ρq_tot = integrator.p.precomputed.ᶜS_ρq_tot
            ᶜS_ρq_tot .= FT(0)
        else
            ᶜS_ρq_tot = integrator.p.precipitation.ᶜS_ρq_tot
            ᶜ3d_rain = integrator.p.precipitation.ᶜ3d_rain
            ᶜ3d_snow = integrator.p.precipitation.ᶜ3d_snow
            ᶜS_ρq_tot .= FT(0)
            ᶜ3d_rain .= FT(0)
            ᶜ3d_snow .= FT(0)
        end
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

# Extension of CA.set_surface_albedo! to set the surface albedo to 0.38 at the beginning of the simulation,
# so the initial callback initialization doesn't lead to NaNs in the radiation model.
function CA.set_surface_albedo!(Y, p, t, ::CA.CouplerAlbedo)
    if float(t) == 0
        FT = eltype(Y)
        # set initial insolation initial conditions
        !(p.atmos.insolation isa CA.IdealizedInsolation) && CA.set_insolation_variables!(Y, p, t, p.atmos.insolation)
        # set surface albedo to 0.38
        @warn "Setting surface albedo to 0.38 at the beginning of the simulation"
        p.radiation.rrtmgp_model.direct_sw_surface_albedo .= FT(0.38)
        p.radiation.rrtmgp_model.diffuse_sw_surface_albedo .= FT(0.38)
    else
        nothing
    end
end

"""
    get_surface_space(sim::ClimaAtmosSimulation)

Get the surface space of the atmosphere simulation, which is the bottom
level of the face space.
"""
get_surface_space(sim::ClimaAtmosSimulation) =
    axes(CC.Spaces.level(CC.Fields.coordinate_field(sim.domain.face_space).z, CC.Utilities.half))

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

function Checkpointer.restore_cache!(sim::ClimaAtmosSimulation, new_cache)
    comms_ctx = ClimaComms.context(sim.integrator.u.c)
    restore!(
        Checkpointer.get_model_cache(sim),
        new_cache,
        comms_ctx;
        ignore = Set([:rc, :params, :ghost_buffer, :hyperdiffusion_ghost_buffer, :data_handler, :graph_context]),
    )
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

        (; face_lw_flux_dn, face_lw_flux_up, face_sw_flux_dn, face_sw_flux_up) = sim.integrator.p.radiation.rrtmgp_model

        LWd_TOA = CC.Fields.level(CC.Fields.array2field(FT.(face_lw_flux_dn), face_space), nz_faces - CC.Utilities.half)
        LWu_TOA = CC.Fields.level(CC.Fields.array2field(FT.(face_lw_flux_up), face_space), nz_faces - CC.Utilities.half)
        SWd_TOA = CC.Fields.level(CC.Fields.array2field(FT.(face_sw_flux_dn), face_space), nz_faces - CC.Utilities.half)
        SWu_TOA = CC.Fields.level(CC.Fields.array2field(FT.(face_sw_flux_up), face_space), nz_faces - CC.Utilities.half)

        return @. -(LWd_TOA + SWd_TOA - LWu_TOA - SWu_TOA)
    else
        return FT[0]
    end
end

function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:energy})
    integrator = sim.integrator
    p = integrator.p


    # return total energy and (if Microphysics0Moment) the energy lost due to precipitation removal
    microphysics_model =
        pkgversion(CA) < v"0.31.0" ? integrator.p.atmos.precip_model : integrator.p.atmos.microphysics_model
    if microphysics_model isa CA.Microphysics0Moment
        ᶜts = p.precomputed.ᶜts
        ᶜΦ = p.core.ᶜΦ
        if pkgversion(CA) >= v"0.29.0"
            ᶜS_ρq_tot = p.precomputed.ᶜS_ρq_tot
        else
            ᶜS_ρq_tot = p.precipitation.ᶜS_ρq_tot
        end
        thermo_params = get_thermo_params(sim)
        return integrator.u.c.ρe_tot .-
               ᶜS_ρq_tot .* CA.e_tot_0M_precipitation_sources_helper.(Ref(thermo_params), ᶜts, ᶜΦ) .*
               float(integrator.dt)
    else
        return integrator.u.c.ρe_tot
    end
end

# helpers for get_field extensions, dipatchable on different moisture model options and radiation modes

surface_rain_flux(::CA.DryModel, integrator) = eltype(integrator.u)(0)
function surface_rain_flux(::Union{CA.EquilMoistModel, CA.NonEquilMoistModel}, integrator)
    if pkgversion(CA) >= v"0.29.0"
        return integrator.p.precomputed.surface_rain_flux
    else
        return integrator.p.precipitation.surface_rain_flux
    end
end

surface_snow_flux(::CA.DryModel, integrator) = eltype(integrator.u)(0)
function surface_snow_flux(::Union{CA.EquilMoistModel, CA.NonEquilMoistModel}, integrator)
    if pkgversion(CA) >= v"0.29.0"
        return integrator.p.precomputed.surface_snow_flux
    else
        return integrator.p.precipitation.surface_snow_flux
    end
end

surface_radiation_flux(::Nothing, integrator) = eltype(integrator.u)(0)
surface_radiation_flux(::CA.RRTMGPI.AbstractRRTMGPMode, integrator) =
    CC.Fields.level(integrator.p.radiation.ᶠradiation_flux, CC.Utilities.half)

moisture_flux(::CA.DryModel, integrator) = eltype(integrator.u)(0)
moisture_flux(::Union{CA.EquilMoistModel, CA.NonEquilMoistModel}, integrator) =
    CC.Geometry.WVector.(integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot)

ρq_tot(::CA.DryModel, integrator) = eltype(integrator.u)(0)
ρq_tot(::Union{CA.EquilMoistModel, CA.NonEquilMoistModel}, integrator) = integrator.u.c.ρq_tot

# extensions required by the Interfacer
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:air_pressure}) =
    CC.Fields.level(sim.integrator.p.precomputed.ᶜp, 1)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:air_temperature}) =
    TD.air_temperature.(
        sim.integrator.p.params.thermodynamics_params,
        CC.Fields.level(sim.integrator.p.precomputed.ᶜts, 1),
    )
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:air_density}) =
    TD.air_density.(sim.integrator.p.params.thermodynamics_params, CC.Fields.level(sim.integrator.p.precomputed.ᶜts, 1))
# When using the `FixedCO2` option, CO2 is stored as that struct type,
# so we access its value (a scalar) from the struct.
# When using the `MaunaLoa` option, CO2 is stored as a 1D Array in the tracers
# cache, so we access it from there.
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:co2}) =
    sim.integrator.p.atmos.co2 isa CA.FixedCO2 ? sim.integrator.p.atmos.co2.value : sim.integrator.p.tracers.co2[1]

function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:diffuse_fraction})
    hasradiation(sim.integrator) || return nothing
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
                ((x, y) -> y > zero(y) ? x / y : zero(y)).(total_flux_dn .- direct_flux_dn, total_flux_dn),
                zero(FT),
                one(FT),
            )
    end
    return CC.Fields.array2field(diffuse_fraction, lowest_face_space)
end
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:liquid_precipitation}) =
    surface_rain_flux(sim.integrator.p.atmos.moisture_model, sim.integrator)
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:LW_d})
    hasradiation(sim.integrator) || return nothing
    return CC.Fields.level(
        CC.Fields.array2field(sim.integrator.p.radiation.rrtmgp_model.face_lw_flux_dn, axes(sim.integrator.u.f)),
        CC.Utilities.half,
    )
end
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:specific_humidity}) =
    CC.Fields.level(ρq_tot(sim.integrator.p.atmos.moisture_model, sim.integrator) ./ sim.integrator.u.c.ρ, 1)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux_sfc}) =
    surface_radiation_flux(sim.integrator.p.atmos.radiation_mode, sim.integrator)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:snow_precipitation}) =
    surface_snow_flux(sim.integrator.p.atmos.moisture_model, sim.integrator)
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:SW_d})
    hasradiation(sim.integrator) || return nothing
    return CC.Fields.level(
        CC.Fields.array2field(sim.integrator.p.radiation.rrtmgp_model.face_sw_flux_dn, axes(sim.integrator.u.f)),
        CC.Utilities.half,
    )
end
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:water}) =
    ρq_tot(sim.integrator.p.atmos.moisture_model, sim.integrator)
function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_temperature}, csf)
    # NOTE: This update_field! takes as argument the entire coupler fields struct, instead
    # of a single field. This is unlike most of other functions, so we may want to revisit
    # it.
    # The rrtmgp_model.surface_temperature field is updated by the RRTMGP callback using the
    # sfc_conditions.ts, so all we need to do is update sfc_conditions.ts
    # Note also that this method updates the entire thermo state, not just the temperature

    # Remap coupler fields onto the atmosphere surface space
    atmos_surface_space = get_surface_space(sim)
    temp_field_surface = sim.integrator.p.scratch.ᶠtemp_field_level
    @assert axes(temp_field_surface) == atmos_surface_space

    # NOTE: This is allocating! If we had 2 more scratch fields, we could avoid this
    T_sfc_atmos = Interfacer.remap(csf.T_sfc, atmos_surface_space)
    q_sfc_atmos = Interfacer.remap(csf.q_sfc, atmos_surface_space)
    # Store `ρ_sfc_atmos` in an atmosphere scratch field on the surface space
    temp_field_surface =
        FluxCalculator.extrapolate_ρ_to_sfc.(
            get_thermo_params(sim),
            sim.integrator.p.precomputed.sfc_conditions.ts,
            T_sfc_atmos,
        ) # ρ_sfc_atmos

    if sim.integrator.p.atmos.moisture_model isa CA.DryModel
        sim.integrator.p.precomputed.sfc_conditions.ts .=
            TD.PhaseDry_ρT.(get_thermo_params(sim), temp_field_surface, T_sfc_atmos) # temp_field_surface = ρ_sfc_atmos
    else
        sim.integrator.p.precomputed.sfc_conditions.ts .=
            TD.PhaseNonEquil_ρTq.(
                get_thermo_params(sim),
                temp_field_surface,
                T_sfc_atmos,
                TD.PhasePartition.(q_sfc_atmos),
            ) # temp_field_surface = ρ_sfc_atmos
    end
end
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:height_int}) =
    CC.Spaces.level(CC.Fields.coordinate_field(sim.integrator.u.c).z, 1)
Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:height_sfc}) =
    CC.Spaces.level(CC.Fields.coordinate_field(sim.integrator.u.f).z, CC.Utilities.half)
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:u_int})
    # NOTE: This calculation is copied from ClimaAtmos (and is allocating! Fix me if you can!)
    int_local_geometry_values = Fields.level(Fields.local_geometry_field(sim.integrator.u.c), 1)
    int_u_values = CC.Spaces.level(sim.integrator.p.precomputed.ᶜu, 1)
    return CA.projected_vector_data.(CA.CT1, int_u_values, int_local_geometry_values)
end
function Interfacer.get_field(sim::ClimaAtmosSimulation, ::Val{:v_int})
    # NOTE: This calculation is copied from ClimaAtmos (and is allocating! Fix me if you can!)
    int_local_geometry_values = Fields.level(Fields.local_geometry_field(sim.integrator.u.c), 1)
    int_u_values = CC.Spaces.level(sim.integrator.p.precomputed.ᶜu, 1)
    return CA.projected_vector_data.(CA.CT2, int_u_values, int_local_geometry_values)
end

# extensions required by the Interfacer
function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:emissivity}, field)
    hasradiation(sim.integrator) || return nothing
    # Remap field onto the atmosphere surface space in scratch field
    temp_field_surface = sim.integrator.p.scratch.ᶠtemp_field_level
    Interfacer.remap!(temp_field_surface, field)
    # sets all 16 bands (rows) to the same values
    sim.integrator.p.radiation.rrtmgp_model.surface_emissivity .=
        reshape(CC.Fields.field2array(temp_field_surface), 1, length(parent(temp_field_surface)))
end
function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_direct_albedo}, field)
    hasradiation(sim.integrator) || return nothing
    # Remap field onto the atmosphere surface space in scratch field
    temp_field_surface = sim.integrator.p.scratch.ᶠtemp_field_level
    Interfacer.remap!(temp_field_surface, field)
    sim.integrator.p.radiation.rrtmgp_model.direct_sw_surface_albedo .=
        reshape(CC.Fields.field2array(temp_field_surface), 1, length(parent(temp_field_surface)))
end

function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_diffuse_albedo}, field)
    hasradiation(sim.integrator) || return nothing
    # Remap field onto the atmosphere surface space in scratch field
    temp_field_surface = sim.integrator.p.scratch.ᶠtemp_field_level
    Interfacer.remap!(temp_field_surface, field)
    sim.integrator.p.radiation.rrtmgp_model.diffuse_sw_surface_albedo .=
        reshape(CC.Fields.field2array(temp_field_surface), 1, length(parent(temp_field_surface)))
end

function Interfacer.update_field!(sim::ClimaAtmosSimulation, ::Val{:turbulent_fluxes}, fields)
    (; F_lh, F_sh, F_turb_moisture, F_turb_ρτxz, F_turb_ρτyz) = fields
    atmos_surface_space = get_surface_space(sim)
    temp_field_surface = sim.integrator.p.scratch.ᶠtemp_field_level

    Y = sim.integrator.u
    surface_local_geometry = CC.Fields.level(CC.Fields.local_geometry_field(Y.f), CC.Utilities.half)
    surface_normal = @. CA.C3(CA.unit_basis_vector_data(CA.C3, surface_local_geometry))

    # get template objects for the contravariant components of the momentum fluxes (required by Atmos boundary conditions)
    # NOTE: This is allocating! Fix me if you can!
    vec_ct12_ct1 = @. CA.CT12(CA.CT1(CA.unit_basis_vector_data(CA.CT1, surface_local_geometry)), surface_local_geometry)
    vec_ct12_ct2 = @. CA.CT12(CA.CT2(CA.unit_basis_vector_data(CA.CT2, surface_local_geometry)), surface_local_geometry)

    # NOTE: This is allocating (F_turb_ρτyz_atmos)! If we had 1 more scratch field, we could avoid this
    Interfacer.remap!(temp_field_surface, F_turb_ρτxz) # F_turb_ρτxz_atmos
    F_turb_ρτyz_atmos = Interfacer.remap(F_turb_ρτyz, atmos_surface_space) # F_turb_ρτyz_atmos
    sim.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ .= (
        surface_normal .⊗
        CA.C12.(temp_field_surface .* vec_ct12_ct1 .+ F_turb_ρτyz_atmos .* vec_ct12_ct2, surface_local_geometry)
    )

    # Remap summed turbulent energy fluxes into the scratch field
    Interfacer.remap!(temp_field_surface, F_lh .+ F_sh)
    sim.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot .= temp_field_surface .* surface_normal # (shf + lhf)

    # Remap turbulent moisture flux into the scratch field
    Interfacer.remap!(temp_field_surface, F_turb_moisture)
    sim.integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot .= temp_field_surface .* surface_normal # (evap)

    Interfacer.remap!(sim.integrator.p.precomputed.sfc_conditions.obukhov_length, fields.L_MO)
    Interfacer.remap!(sim.integrator.p.precomputed.sfc_conditions.ustar, fields.ustar)
    Interfacer.remap!(sim.integrator.p.precomputed.sfc_conditions.buoyancy_flux, fields.buoyancy_flux)

    return nothing
end

# extensions required by FieldExchanger
Interfacer.step!(sim::ClimaAtmosSimulation, t::Real) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
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
function Interfacer.add_coupler_fields!(coupler_field_names, atmos_sim::ClimaAtmosSimulation)
    atmos_coupler_fields = [:surface_direct_albedo, :surface_diffuse_albedo, :ustar, :L_MO, :buoyancy_flux]
    push!(coupler_field_names, atmos_coupler_fields...)
end

"""
    Interfacer.close_output_writers(sim::ClimaAtmosSimulation)

Close all output writers used by the atmos simulation.
"""
Interfacer.close_output_writers(sim::ClimaAtmosSimulation) =
    isnothing(sim.output_writers) || foreach(close, sim.output_writers)


# FluxCalculator.get_surface_params required by FluxCalculator (partitioned fluxes)
FluxCalculator.get_surface_params(sim::ClimaAtmosSimulation) = CAP.surface_fluxes_params(sim.integrator.p.params)

"""
    Interfacer.set_cache!(sim::ClimaAtmosSimulation)

Set cache variables that cannot be initialized before the initial exchange.
Radiation must be called here because it requires the surface model initial
temperatures, which are set in the atmosphere during the initial exchange.
Any other callbacks which modify the cache should be called here as well.

This function does not set all the cache variables, because many are computed
as part of the tendendencies.
"""
function Interfacer.set_cache!(sim::ClimaAtmosSimulation)
    if hasradiation(sim.integrator)
        CA.rrtmgp_model_callback!(sim.integrator)
        if pkgversion(CA) == v"0.30" && !isnothing(sim.integrator.p.atmos.non_orographic_gravity_wave)
            # In version 0.30, nogw_model_callback crashes when there are no gravity waves,
            # see CA #3792
            CA.nogw_model_callback!(sim.integrator)
        else
            pkgversion(CA) > v"0.30" && CA.nogw_model_callback!(sim.integrator)
        end
    end
    return nothing
end

###
### ClimaAtmos.jl model-specific functions (not explicitly required by ClimaCoupler.jl)
###
"""
    get_atmos_config_dict(coupler_dict::Dict, atmos_output_dir)

Returns the specified atmospheric configuration (`atmos_config`) overwitten by arguments
in the coupler dictionary (`config_dict`).
The returned `atmos_config` dictionary will then be used to set up the atmosphere simulation.

The `atmos_config_repo` flag allows us to use a configuration specified within
either the ClimaCoupler or ClimaAtmos repos, which is useful for direct
coupled/atmos-only comparisons.

In this function, parameters are overwritten in a specific order, from lowest to highest priority:
    1. Default atmos config
    2. Provided atmos config file (if any)
    3. TOML parameter file
    5. Output directory and timestep are set explicitly based on the coupler config

The TOML parameter file to use is chosen using the following priority:
If a coupler TOML file is provided, it is used. Otherwise we use an atmos TOML
file if it's provided. If neither is provided, we use a default coupler TOML file.
"""
function get_atmos_config_dict(coupler_dict::Dict, atmos_output_dir)
    atmos_config_file = coupler_dict["atmos_config_file"]
    atmos_config_repo = coupler_dict["atmos_config_repo"]
    # override default or specified configs with coupler arguments, and set the correct atmos config_file
    if atmos_config_repo == "ClimaCoupler"
        @assert !isnothing(atmos_config_file) "Must specify `atmos_config_file` within ClimaCoupler."

        @info "Using Atmos configuration from ClimaCoupler in $atmos_config_file"
        atmos_config = merge(
            CA.override_default_config(joinpath(pkgdir(ClimaCoupler), atmos_config_file)),
            coupler_dict,
            Dict("config_file" => atmos_config_file),
        )
    elseif atmos_config_repo == "ClimaAtmos"
        if isnothing(atmos_config_file)
            @info "Using Atmos default configuration"
            atmos_config = merge(CA.default_config_dict(), coupler_dict, Dict("config_file" => atmos_config_file))
        else
            @info "Using Atmos configuration from $atmos_config_file"
            atmos_config = merge(
                CA.override_default_config(joinpath(pkgdir(CA), atmos_config_file)),
                coupler_dict,
                Dict("config_file" => atmos_config_file),
            )
        end
    else
        error("Invalid `atmos_config_repo`; please use \"ClimaCoupler\" or \"ClimaAtmos\"")
    end

    # use atmos toml if coupler toml is not defined
    # If we can't find the file at the relative path, prepend pkgdir(ClimaAtmos)
    atmos_toml = map(atmos_config["toml"]) do file
        isfile(file) ? file : joinpath(pkgdir(CA), file)
    end
    coupler_toml = atmos_config["coupler_toml"]
    toml = isempty(coupler_toml) ? atmos_toml : coupler_toml
    if !isempty(toml)
        @info "Overwriting Atmos parameters from input TOML file(s): $toml"
        atmos_config["toml"] = toml
    end

    # Specify atmos output directory to be inside the coupler output directory
    atmos_config["output_dir_style"] = "RemovePreexisting"
    atmos_config["output_dir"] = atmos_output_dir

    # Ensure Atmos's own checkpoints are synced up with ClimaCoupler, so that we
    # can pick up from where we have left. NOTE: This should not be needed, but
    # there is no easy way to initialize ClimaAtmos with a different t_start
    atmos_config["dt_save_state_to_disk"] = coupler_dict["checkpoint_dt"]

    # Add all extra atmos diagnostic entries into the vector of atmos diagnostics
    atmos_config["diagnostics"] =
        haskey(atmos_config, "diagnostics") ?
        vcat(atmos_config["diagnostics"], coupler_dict["extra_atmos_diagnostics"]) :
        coupler_dict["extra_atmos_diagnostics"]

    # The Atmos `get_simulation` function expects the atmos config to contains its timestep size
    # in the `dt` field. If there is a `dt_atmos` field in coupler_dict, we add it to the atmos config as `dt`
    dt_atmos = haskey(coupler_dict, "dt_atmos") ? coupler_dict["dt_atmos"] : coupler_dict["dt"]
    atmos_config["dt"] = dt_atmos

    # set restart file to the initial file saved in this location if it is not nothing
    # TODO this is hardcoded and should be fixed once we have a better restart system
    if !isnothing(atmos_config["restart_file"])
        atmos_config["restart_file"] = replace(atmos_config["restart_file"], "active" => "0000")
    end
    return atmos_config
end

"""
    get_thermo_params(sim::ClimaAtmosSimulation)

Returns the thermodynamic parameters from the atmospheric model simulation object.
"""
get_thermo_params(sim::ClimaAtmosSimulation) = CAP.thermodynamics_params(sim.integrator.p.params)

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
