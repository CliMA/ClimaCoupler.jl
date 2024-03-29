# atmos_init: for ClimaAtmos pre-AMIP interface
using StaticArrays
using Statistics: mean
using LinearAlgebra: norm

import ClimaAtmos as CA
import ClimaAtmos: CT1, CT2, CT12, CT3, C3, C12, unit_basis_vector_data, ⊗
import SurfaceFluxes as SF
using ClimaCore
using ClimaCore.Utilities: half

import ClimaCoupler.Interfacer: AtmosModelSimulation
import ClimaCoupler.FluxCalculator:
    atmos_turbulent_fluxes!,
    calculate_surface_air_density,
    PartitionedStateFluxes,
    extrapolate_ρ_to_sfc,
    get_surface_params,
    water_albedo_from_atmosphere!
import ClimaCoupler.Interfacer: get_field, update_field!, name
import ClimaCoupler.Checkpointer: get_model_prog_state
import ClimaCoupler.FieldExchanger: update_sim!, step!, reinit!
import ClimaCoupler.Utilities: swap_space!

include("climaatmos_extra_diags.jl")

###
### Functions required by ClimaCoupler.jl for an AtmosModelSimulation
###
struct ClimaAtmosSimulation{P, Y, D, I} <: AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
name(::ClimaAtmosSimulation) = "ClimaAtmosSimulation"

function atmos_init(::Type{FT}, atmos_config_dict::Dict) where {FT}
    # By passing `parsed_args` to `AtmosConfig`, `parsed_args` overwrites the default atmos config
    atmos_config_dict["surface_albedo"] = "CouplerAlbedo"
    atmos_config = CA.AtmosConfig(atmos_config_dict)
    simulation = CA.get_simulation(atmos_config)
    (; integrator) = simulation
    Y = integrator.u
    center_space = axes(Y.c.ρe_tot)
    face_space = axes(Y.f.u₃)
    spaces = (; center_space = center_space, face_space = face_space)
    if :ρe_int in propertynames(Y.c)
        @warn("Running with ρe_int in coupled mode is not tested yet.", maxlog = 1)
    end

    # define shorter references for long variable names to increase readability
    ρ_flux_h_tot = integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot
    ρ_flux_q_tot = integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot
    ᶠradiation_flux = integrator.p.radiation.ᶠradiation_flux
    ρ_flux_uₕ = integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ
    col_integrated_rain = integrator.p.precipitation.col_integrated_rain
    col_integrated_snow = integrator.p.precipitation.col_integrated_snow
    ᶜS_ρq_tot = integrator.p.precipitation.ᶜS_ρq_tot
    ᶜ3d_rain = integrator.p.precipitation.ᶜ3d_rain
    ᶜ3d_snow = integrator.p.precipitation.ᶜ3d_snow

    # set initial fluxes to zero
    @. ρ_flux_h_tot = ClimaCore.Geometry.Covariant3Vector(FT(0.0))
    @. ρ_flux_q_tot = ClimaCore.Geometry.Covariant3Vector(FT(0.0))
    @. ᶠradiation_flux = ClimaCore.Geometry.WVector(FT(0))
    ρ_flux_uₕ.components .= Ref(SMatrix{1, 2}([FT(0), FT(0)]))
    col_integrated_rain .= FT(0)
    col_integrated_snow .= FT(0)
    ᶜS_ρq_tot .= FT(0)
    ᶜ3d_rain .= FT(0)
    ᶜ3d_snow .= FT(0)

    sim = ClimaAtmosSimulation(integrator.p.params, Y, spaces, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

"""
    get_model_prog_state(sim::ClimaAtmosSimulation)

Extension of Checkpointer.get_model_prog_state to get the model state.
"""
function get_model_prog_state(sim::ClimaAtmosSimulation)
    return sim.integrator.u
end

"""
    get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux_toa})

Extension of Interfacer.get_field to get the net TOA radiation, which is a sum of the
upward and downward longwave and shortwave radiation.
"""
function get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux_toa})
    FT = eltype(atmos_sim.integrator.u)

    if atmos_sim.integrator.p.radiation.radiation_model != nothing
        face_space = axes(atmos_sim.integrator.u.f)
        nz_faces = length(ClimaCore.Spaces.vertical_topology(face_space).mesh.faces)

        (; face_lw_flux_dn, face_lw_flux_up, face_sw_flux_dn, face_sw_flux_up) =
            atmos_sim.integrator.p.radiation.radiation_model

        LWd_TOA = ClimaCore.Fields.level(CA.RRTMGPI.array2field(FT.(face_lw_flux_dn), face_space), nz_faces - half)
        LWu_TOA = ClimaCore.Fields.level(CA.RRTMGPI.array2field(FT.(face_lw_flux_up), face_space), nz_faces - half)
        SWd_TOA = ClimaCore.Fields.level(CA.RRTMGPI.array2field(FT.(face_sw_flux_dn), face_space), nz_faces - half)
        SWu_TOA = ClimaCore.Fields.level(CA.RRTMGPI.array2field(FT.(face_sw_flux_up), face_space), nz_faces - half)

        return @. -(LWd_TOA + SWd_TOA - LWu_TOA - SWu_TOA)
    else
        return FT(0)
    end
end

function get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:energy})
    thermo_params = get_thermo_params(atmos_sim)

    ᶜS_ρq_tot = atmos_sim.integrator.p.precipitation.ᶜS_ρq_tot
    ᶜts = atmos_sim.integrator.p.precomputed.ᶜts
    ᶜΦ = atmos_sim.integrator.p.core.ᶜΦ

    # return total energy and (if Microphysics0Moment) the energy lost due to precipitation removal
    if atmos_sim.integrator.p.atmos.precip_model isa CA.Microphysics0Moment
        ᶜS_ρq_tot = atmos_sim.integrator.p.precipitation.ᶜS_ρq_tot
        ᶜts = atmos_sim.integrator.p.precomputed.ᶜts
        ᶜΦ = atmos_sim.integrator.p.core.ᶜΦ
        return atmos_sim.integrator.u.c.ρe_tot .-
               ᶜS_ρq_tot .* CA.e_tot_0M_precipitation_sources_helper.(Ref(thermo_params), ᶜts, ᶜΦ) .*
               atmos_sim.integrator.dt
    else
        return atmos_sim.integrator.u.c.ρe_tot
    end
end

# extensions required by the Interfacer
get_field(sim::ClimaAtmosSimulation, ::Val{:air_density}) =
    TD.air_density.(thermo_params, sim.integrator.p.precomputed.ᶜts)
get_field(sim::ClimaAtmosSimulation, ::Val{:air_temperature}) =
    TD.air_temperature.(thermo_params, sim.integrator.p.precomputed.ᶜts)
get_field(sim::ClimaAtmosSimulation, ::Val{:liquid_precipitation}) = sim.integrator.p.precipitation.col_integrated_rain
get_field(sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux_sfc}) =
    ClimaCore.Fields.level(sim.integrator.p.radiation.ᶠradiation_flux, half)
get_field(sim::ClimaAtmosSimulation, ::Val{:snow_precipitation}) = sim.integrator.p.precipitation.col_integrated_snow
get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_energy_flux}) =
    ClimaCore.Geometry.WVector.(sim.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot)
get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_moisture_flux}) =
    ClimaCore.Geometry.WVector.(sim.integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot)
get_field(sim::ClimaAtmosSimulation, ::Val{:thermo_state_int}) =
    ClimaCore.Spaces.level(sim.integrator.p.precomputed.ᶜts, 1)
get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:water}) = atmos_sim.integrator.u.c.ρq_tot

# extensions required by FluxCalculator (partitioned fluxes)
get_field(sim::ClimaAtmosSimulation, ::Val{:height_int}) =
    ClimaCore.Spaces.level(ClimaCore.Fields.coordinate_field(sim.integrator.u.c).z, 1)
get_field(sim::ClimaAtmosSimulation, ::Val{:height_sfc}) =
    ClimaCore.Spaces.level(ClimaCore.Fields.coordinate_field(sim.integrator.u.f).z, half)
function get_field(sim::ClimaAtmosSimulation, ::Val{:uv_int})
    uₕ_int = ClimaCore.Geometry.UVVector.(ClimaCore.Spaces.level(sim.integrator.u.c.uₕ, 1))
    return @. StaticArrays.SVector(uₕ_int.components.data.:1, uₕ_int.components.data.:2)
end

function update_field!(atmos_sim::ClimaAtmosSimulation, ::Val{:co2}, field)
    if atmos_sim.integrator.p.atmos.radiation_mode isa CA.RRTMGPI.GrayRadiation
        @warn("Gray radiation model initialized, skipping CO2 update", maxlog = 1)
        return
    else
        atmos_sim.integrator.p.radiation.radiation_model.volume_mixing_ratio_co2 .= mean(parent(field))
    end
end
# extensions required by the Interfacer
function update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_temperature}, csf)
    sim.integrator.p.radiation.radiation_model.surface_temperature .= CA.RRTMGPI.field2array(csf.T_S)
end

function update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_direct_albedo}, field)
    sim.integrator.p.radiation.radiation_model.direct_sw_surface_albedo .=
        reshape(CA.RRTMGPI.field2array(field), 1, length(parent(field)))
end

function update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_diffuse_albedo}, field)
    sim.integrator.p.radiation.radiation_model.diffuse_sw_surface_albedo .=
        reshape(CA.RRTMGPI.field2array(field), 1, length(parent(field)))
end

function update_field!(sim::ClimaAtmosSimulation, ::Val{:turbulent_fluxes}, fields)
    (; F_turb_energy, F_turb_moisture, F_turb_ρτxz, F_turb_ρτyz) = fields

    Y = sim.integrator.u
    surface_local_geometry = ClimaCore.Fields.level(ClimaCore.Fields.local_geometry_field(Y.f), ClimaCore.Fields.half)
    surface_normal = @. C3(unit_basis_vector_data(C3, surface_local_geometry))

    # get template objects for the contravariant components of the momentum fluxes (required by Atmos boundary conditions)
    vec_ct12_ct1 = @. CT12(CT2(unit_basis_vector_data(CT1, surface_local_geometry)), surface_local_geometry)
    vec_ct12_ct2 = @. CT12(CT2(unit_basis_vector_data(CT2, surface_local_geometry)), surface_local_geometry)

    sim.integrator.p.precomputed.sfc_conditions.ρ_flux_uₕ .= (
        surface_normal .⊗
        C12.(
            swap_space!(ones(axes(vec_ct12_ct1)), F_turb_ρτxz) .* vec_ct12_ct1 .+
            swap_space!(ones(axes(vec_ct12_ct2)), F_turb_ρτyz) .* vec_ct12_ct2,
            surface_local_geometry,
        )
    )

    parent(sim.integrator.p.precomputed.sfc_conditions.ρ_flux_h_tot) .= parent(F_turb_energy) .* parent(surface_normal) # (shf + lhf)
    parent(sim.integrator.p.precomputed.sfc_conditions.ρ_flux_q_tot) .=
        parent(F_turb_moisture) .* parent(surface_normal) # (evap)

    # TODO: see if Atmos can rever to a simpler solution
end

# extensions required by FieldExchanger
step!(sim::ClimaAtmosSimulation, t) = step!(sim.integrator, t - sim.integrator.t, true)
reinit!(sim::ClimaAtmosSimulation) = reinit!(sim.integrator)

function update_sim!(atmos_sim::ClimaAtmosSimulation, csf, turbulent_fluxes)
    update_field!(atmos_sim, Val(:surface_direct_albedo), csf.surface_direct_albedo)
    update_field!(atmos_sim, Val(:surface_diffuse_albedo), csf.surface_diffuse_albedo)
    update_field!(atmos_sim, Val(:surface_temperature), csf)

    if turbulent_fluxes isa PartitionedStateFluxes
        update_field!(atmos_sim, Val(:turbulent_fluxes), csf)
    end
end

"""
    calculate_surface_air_density(atmos_sim::ClimaAtmosSimulation, T_S::ClimaCore.Fields.Field)

Extension for this function to calculate surface density.
"""
function calculate_surface_air_density(atmos_sim::ClimaAtmosSimulation, T_S::ClimaCore.Fields.Field)
    thermo_params = get_thermo_params(atmos_sim)
    ts_int = get_field(atmos_sim, Val(:thermo_state_int))
    extrapolate_ρ_to_sfc.(Ref(thermo_params), ts_int, swap_space!(ones(axes(ts_int.ρ)), T_S))
end

# get_surface_params required by FluxCalculator (partitioned fluxes)
get_surface_params(sim::ClimaAtmosSimulation) = CAP.surface_fluxes_params(sim.integrator.p.params)

###
### ClimaAtmos.jl model-specific functions (not explicitly required by ClimaCoupler.jl)
###
"""
    get_atmos_config(coupler_dict::Dict)

Returns the specified atmospheric configuration (`atmos_config_dict`) overwitten by arguments
in the coupler dictionary (`config_dict`).
"""
function get_atmos_config(coupler_dict)
    atmos_config_file = coupler_dict["atmos_config_file"]
    # override default or specified configs with coupler arguments, and set the correct atmos config_file
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

    # use coupler toml if atmos is not defined
    atmos_toml_file = atmos_config["toml"]
    coupler_toml_file = coupler_dict["coupler_toml_file"]
    default_toml_file = "toml/default_coarse.toml"

    toml_file = !isempty(atmos_toml_file) ? joinpath(pkgdir(CA), atmos_toml_file[1]) : nothing
    toml_file = !isnothing(coupler_toml_file) ? joinpath(pkgdir(ClimaCoupler), coupler_toml_file) : toml_file
    toml_file = isnothing(toml_file) ? joinpath(pkgdir(ClimaCoupler), default_toml_file) : toml_file

    if !isnothing(toml_file)
        @info "Overwriting Atmos parameters from $toml_file"
        atmos_config = merge(atmos_config, Dict("toml" => [toml_file]))
    end

    # specify atmos output directory to be inside the coupler output directory
    atmos_output_dir = joinpath(
        coupler_dict["coupler_output_dir"],
        joinpath(coupler_dict["mode_name"], coupler_dict["run_name"]),
        "clima_atmos",
    )
    atmos_config = merge(atmos_config, Dict("output_dir" => atmos_output_dir))

    return atmos_config
end


# flux calculation borrowed from atmos
"""
    CoupledMoninObukhov()
A modified version of a Monin-Obukhov surface for the Coupler, see the link below for more information
https://clima.github.io/SurfaceFluxes.jl/dev/SurfaceFluxes/#Monin-Obukhov-Similarity-Theory-(MOST)
"""
struct CoupledMoninObukhov end
"""
    coupler_surface_setup(::CoupledMoninObukhov, p, csf_sfc = (; T = nothing, z0m = nothing, z0b = nothing, beta = nothing, q_vap = nothing))

Sets up `surface_setup` as a `ClimaCore.Fields.Field` of `SurfaceState`s.
"""
function coupler_surface_setup(
    ::CoupledMoninObukhov,
    p;
    csf_sfc = (; T = nothing, z0m = nothing, z0b = nothing, beta = nothing, q_vap = nothing),
)

    surface_state(z0m, z0b, T, beta, q_vap) = CA.SurfaceConditions.SurfaceState(;
        parameterization = CA.SurfaceConditions.MoninObukhov(; z0m, z0b),
        T,
        beta,
        q_vap,
    )
    surface_state_field = @. surface_state(csf_sfc.z0m, csf_sfc.z0b, csf_sfc.T, csf_sfc.beta, csf_sfc.q_vap)
    return surface_state_field
end

"""
    get_new_cache(atmos_sim::ClimaAtmosSimulation, csf)

Returns a new `p` with the updated surface conditions.
"""
function get_new_cache(atmos_sim::ClimaAtmosSimulation, csf)
    csf_sfc = (; T = csf.T_S, z0m = csf.z0m_S, z0b = csf.z0b_S, beta = csf.beta, q_vap = csf.q_sfc)
    modified_atmos_cache(atmos_sim, csf_sfc)
end

"""
    modified_atmos_cache(atmos_sim, csf_sfc)

Returns a new `p` with the updated surface conditions.
"""
function modified_atmos_cache(atmos_sim, csf_sfc)

    p = atmos_sim.integrator.p

    coupler_sfc_setup = coupler_surface_setup(CoupledMoninObukhov(), p; csf_sfc = csf_sfc)

    p_names = propertynames(p)
    p_values = map(x -> x == :sfc_setup ? coupler_sfc_setup : getproperty(p, x), p_names) # TODO: use merge here

    (; zip(p_names, p_values)...)
end

"""
    atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)

Computes turbulent surface fluxes using ClimaAtmos's `update_surface_conditions!`. This
requires that we define a new temporary parameter Tuple, `new_p`, and save the new surface state
in it. We do not want `new_p` to live in the atmospheric model permanently, because that would also
trigger flux calculation during Atmos `step!`. We only want to trigger this once per coupling
timestep from ClimaCoupler.

For debigging atmos, we can set the following atmos defaults:
 csf.z0m_S .= 1.0e-5
 csf.z0b_S .= 1.0e-5
 csf.beta .= 1
 csf = merge(csf, (;q_sfc = nothing))
"""
function atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)

    if isnothing(atmos_sim.integrator.p.sfc_setup) # trigger flux calculation if not done in Atmos internally
        new_p = get_new_cache(atmos_sim, csf)
        CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t)
        atmos_sim.integrator.p.precomputed.sfc_conditions .= new_p.precomputed.sfc_conditions
    end
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
        buffer = ClimaCore.Spaces.create_dss_buffer(field)
        ClimaCore.Spaces.weighted_dss!(field, buffer)
    end
end

"""
    water_albedo_from_atmosphere!(atmos_sim::ClimaAtmosSimulation, direct_albedo::ClimaCore.Fields.Field, diffuse_albedo::ClimaCore.Fields.Field)

Extension to calculate the water surface albedo from wind speed and insolation. It can be used for prescribed ocean and lakes.
"""
function water_albedo_from_atmosphere!(
    atmos_sim::ClimaAtmosSimulation,
    direct_albedo::ClimaCore.Fields.Field,
    diffuse_albedo::ClimaCore.Fields.Field,
)

    Y = atmos_sim.integrator.u
    p = atmos_sim.integrator.p
    t = atmos_sim.integrator.t

    radiation_model = atmos_sim.integrator.p.radiation.radiation_model
    FT = eltype(Y)
    λ = FT(0) # spectral wavelength (not used for now)

    # update for current zenith angle
    CA.set_insolation_variables!(Y, p, t)

    bottom_coords = ClimaCore.Fields.coordinate_field(ClimaCore.Spaces.level(Y.c, 1))
    μ = CA.RRTMGPI.array2field(radiation_model.cos_zenith, axes(bottom_coords))
    FT = eltype(atmos_sim.integrator.u)
    α_model = CA.RegressionFunctionAlbedo{FT}()


    # set the direct and diffuse surface albedos
    direct_albedo .= CA.surface_albedo_direct(α_model).(λ, μ, norm.(ClimaCore.Fields.level(Y.c.uₕ, 1)))
    diffuse_albedo .= CA.surface_albedo_diffuse(α_model).(λ, μ, norm.(ClimaCore.Fields.level(Y.c.uₕ, 1)))

end
