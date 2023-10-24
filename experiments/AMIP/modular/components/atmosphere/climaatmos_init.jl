# atmos_init: for ClimaAtmos pre-AMIP interface
import ClimaAtmos as CA
import ClimaAtmos: CT1, CT2, CT12, CT3, C3, C12, unit_basis_vector_data, ⊗
import ClimaCoupler.FluxCalculator:
    atmos_turbulent_fluxes!,
    calculate_surface_air_density,
    PartitionedStateFluxes,
    extrapolate_ρ_to_sfc,
    get_surface_params
using ClimaCore: Fields.level, Geometry
import ClimaCoupler.FieldExchanger: get_thermo_params
import ClimaCoupler.Interfacer: get_field, update_field!, name, get_model_state_vector
using StaticArrays

import SurfaceFluxes as SF



# the clima atmos `integrator` is now defined
struct ClimaAtmosSimulation{P, Y, D, I} <: AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
name(::ClimaAtmosSimulation) = "ClimaAtmosSimulation"

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
    atmos_toml_file = coupler_dict["atmos_toml_file"]
    coupler_toml_file = coupler_dict["coupler_toml_file"]
    toml_file = isnothing(coupler_toml_file) ? joinpath(pkgdir(CA), atmos_toml_file) : joinpath(pkgdir(ClimaCoupler), coupler_toml_file)

    if !isnothing(toml_file)
        @info "Overwriting Atmos parameters from $toml_file"
        atmos_config = merge(atmos_config, Dict("toml" => [toml_file]))
    end
    return atmos_config
end

function atmos_init(::Type{FT}, atmos_config_dict::Dict) where {FT}

    # By passing `parsed_args` to `AtmosConfig`, `parsed_args` overwrites the default atmos config
    atmos_config = CA.AtmosConfig(atmos_config_dict)
    integrator = CA.get_integrator(atmos_config)
    Y = integrator.u
    center_space = axes(Y.c.ρe_tot)
    face_space = axes(Y.f.u₃)
    spaces = (; center_space = center_space, face_space = face_space)
    if :ρe_int in propertynames(Y.c)
        @warn("Running with ρe_int in coupled mode is not tested yet.", maxlog = 1)
    end

    # set initial fluxes to zero
    @. integrator.p.sfc_conditions.ρ_flux_h_tot = Geometry.Covariant3Vector(FT(0.0))
    @. integrator.p.sfc_conditions.ρ_flux_q_tot = Geometry.Covariant3Vector(FT(0.0))
    @. integrator.p.sfc_conditions.ρ_flux_uₕ.components = zeros(axes(integrator.p.sfc_conditions.ρ_flux_uₕ.components))
    parent(integrator.p.ᶠradiation_flux) .= parent(zeros(axes(integrator.p.ᶠradiation_flux)))
    integrator.p.col_integrated_rain .= FT(0)
    integrator.p.col_integrated_snow .= FT(0)

    sim = ClimaAtmosSimulation(integrator.p.params, Y, spaces, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

# extensions required by the Interfacer
get_field(sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux}) =
    Fields.level(sim.integrator.p.ᶠradiation_flux, half)
function get_field(sim::ClimaAtmosSimulation, ::Val{:liquid_precipitation})
    ρ_liq = CAP.ρ_cloud_liq(sim.integrator.p.params)
    sim.integrator.p.col_integrated_rain .* ρ_liq .+ sim.integrator.p.col_integrated_snow .* ρ_liq # kg/m^2/s
end
function get_field(sim::ClimaAtmosSimulation, ::Val{:snow_precipitation})
    ρ_liq = CAP.ρ_cloud_liq(sim.integrator.p.params)
    sim.integrator.p.col_integrated_snow .* ρ_liq .* 0  # kg/m^2/s
end

get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_energy_flux}) =
    Geometry.WVector.(sim.integrator.p.sfc_conditions.ρ_flux_h_tot)
get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_moisture_flux}) =
    Geometry.WVector.(sim.integrator.p.sfc_conditions.ρ_flux_q_tot)

get_field(sim::ClimaAtmosSimulation, ::Val{:thermo_state_int}) = Spaces.level(sim.integrator.p.ᶜts, 1)

# extensions required by FluxCalculator (partitioned fluxes)
get_field(sim::ClimaAtmosSimulation, ::Val{:height_int}) =
    Spaces.level(Fields.coordinate_field(sim.integrator.u.c).z, 1)
get_field(sim::ClimaAtmosSimulation, ::Val{:height_sfc}) =
    Spaces.level(Fields.coordinate_field(sim.integrator.u.f).z, half)
function get_field(sim::ClimaAtmosSimulation, ::Val{:uv_int})
    uₕ_int = Geometry.UVVector.(Spaces.level(sim.integrator.u.c.uₕ, 1))
    return @. StaticArrays.SVector(uₕ_int.components.data.:1, uₕ_int.components.data.:2)
end
get_field(sim::ClimaAtmosSimulation, ::Val{:air_density}) = TD.air_density.(thermo_params, sim.integrator.p.ᶜts)
get_field(sim::ClimaAtmosSimulation, ::Val{:air_temperature}) = TD.air_temperature.(thermo_params, sim.integrator.p.ᶜts)
get_field(sim::ClimaAtmosSimulation, ::Val{:cv_m}) = TD.cv_m.(thermo_params, sim.integrator.p.ᶜts)
get_field(sim::ClimaAtmosSimulation, ::Val{:gas_constant_air}) =
    TD.gas_constant_air.(thermo_params, sim.integrator.p.ᶜts)

get_surface_params(sim::ClimaAtmosSimulation) = CAP.surface_fluxes_params(sim.integrator.p.params)

function update_field!(atmos_sim::ClimaAtmosSimulation, ::Val{:co2_gm}, field)
    if atmos_sim.integrator.p.radiation_model.radiation_mode isa CA.RRTMGPI.GrayRadiation
        @warn("Gray radiation model initialized, skipping CO2 update", maxlog = 1)
        return
    else
        atmos_sim.integrator.p.radiation_model.volume_mixing_ratio_co2 .= parent(field)[1]
    end
end
# extensions required by the Interfacer
function update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_temperature}, csf)
    sim.integrator.p.radiation_model.surface_temperature .= CA.RRTMGPI.field2array(csf.T_S)
end

function update_field!(sim::ClimaAtmosSimulation, ::Val{:albedo}, field)
    sim.integrator.p.radiation_model.diffuse_sw_surface_albedo .=
        reshape(CA.RRTMGPI.field2array(field), 1, length(parent(field)))
    sim.integrator.p.radiation_model.direct_sw_surface_albedo .=
        reshape(CA.RRTMGPI.field2array(field), 1, length(parent(field)))
end

# get_surface_params required by FluxCalculator (partitioned fluxes)
function update_field!(sim::ClimaAtmosSimulation, ::Val{:turbulent_fluxes}, fields)
    (; F_turb_energy, F_turb_moisture, F_turb_ρτxz, F_turb_ρτyz) = fields

    Y = sim.integrator.u
    surface_local_geometry = Fields.level(Fields.local_geometry_field(Y.f), Fields.half)
    surface_normal = @. C3(unit_basis_vector_data(C3, surface_local_geometry))

    # get template objects for the contravariant components of the momentum fluxes (required by Atmos boundary conditions)
    vec_ct12_ct1 = @. CT12(CT2(unit_basis_vector_data(CT1, surface_local_geometry)), surface_local_geometry)
    vec_ct12_ct2 = @. CT12(CT2(unit_basis_vector_data(CT2, surface_local_geometry)), surface_local_geometry)

    sim.integrator.p.sfc_conditions.ρ_flux_uₕ .= (
        surface_normal .⊗
        C12.(
            swap_space!(ones(axes(vec_ct12_ct1)), F_turb_ρτxz) .* vec_ct12_ct1 .+
            swap_space!(ones(axes(vec_ct12_ct2)), F_turb_ρτyz) .* vec_ct12_ct2,
            surface_local_geometry,
        )
    )

    parent(sim.integrator.p.sfc_conditions.ρ_flux_h_tot) .= parent(F_turb_energy) .* parent(surface_normal) # (shf + lhf)
    parent(sim.integrator.p.sfc_conditions.ρ_flux_q_tot) .= parent(F_turb_moisture) .* parent(surface_normal) # (evap)

    # TODO: see if Atmos can rever to a simpler solution
end

# extensions required by FieldExchanger
step!(sim::ClimaAtmosSimulation, t) = step!(sim.integrator, t - sim.integrator.t, true)

reinit!(sim::ClimaAtmosSimulation) = reinit!(sim.integrator)

function update_sim!(atmos_sim::ClimaAtmosSimulation, csf, turbulent_fluxes)
    update_field!(atmos_sim, Val(:albedo), csf.albedo)
    update_field!(atmos_sim, Val(:surface_temperature), csf)

    if turbulent_fluxes isa PartitionedStateFluxes
        update_field!(atmos_sim, Val(:turbulent_fluxes), csf)
    end
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

Sets up `surface_setup` as a `Fields.Field` of `SurfaceState`s.
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

    p = atmos_sim.integrator.p

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
    get_thermo_params(sim::ClimaAtmosSimulation)

Returns the thermodynamic parameters from the atmospheric model simulation object.
"""
get_thermo_params(sim::ClimaAtmosSimulation) = CAP.thermodynamics_params(sim.integrator.p.params)

"""
    calculate_surface_air_density(atmos_sim::ClimaAtmosSimulation, T_S::Fields.Field)

Extension for this  to to calculate surface density.
"""
function calculate_surface_air_density(atmos_sim::ClimaAtmosSimulation, T_S::Fields.Field)
    thermo_params = get_thermo_params(atmos_sim)
    ts_int = get_field(atmos_sim, Val(:thermo_state_int))
    extrapolate_ρ_to_sfc.(Ref(thermo_params), ts_int, swap_space!(ones(axes(ts_int.ρ)), T_S))
end

"""
    get_model_state_vector(sim::ClimaAtmosSimulation)

Extension of Checkpointer.get_model_state_vector to get the model state.
"""
function get_model_state_vector(sim::ClimaAtmosSimulation)
    return sim.integrator.u
end

"""
    get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:F_radiative_TOA})

Extension of Interfacer.get_field to get the net TOA radiation, which is a sum of the
upward and downward longwave and shortwave radiation.
"""
function get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:F_radiative_TOA})
    radiation = atmos_sim.integrator.p.radiation_model
    FT = eltype(atmos_sim.integrator.u)
    # save radiation source
    if radiation != nothing
        face_space = axes(atmos_sim.integrator.u.f)
        z = parent(Fields.coordinate_field(face_space).z)
        Δz_top = round(FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1]))
        n_faces = length(z[:, 1, 1, 1, 1])

        LWd_TOA = Fields.level(
            CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_dn), face_space),
            n_faces - half,
        )
        LWu_TOA = Fields.level(
            CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_up), face_space),
            n_faces - half,
        )
        SWd_TOA = Fields.level(
            CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_dn), face_space),
            n_faces - half,
        )
        SWu_TOA = Fields.level(
            CA.RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_up), face_space),
            n_faces - half,
        )

        return @. -(LWd_TOA + SWd_TOA - LWu_TOA - SWu_TOA) # [W/m^2]
    else
        return FT(0)
    end
end

get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:energy}) = atmos_sim.integrator.u.c.ρe_tot

get_field(atmos_sim::ClimaAtmosSimulation, ::Val{:water}) = atmos_sim.integrator.u.c.ρq_tot

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
        buffer = Spaces.create_dss_buffer(field)
        Spaces.weighted_dss!(field, buffer)
    end
end



# DEBUG below




"""
    atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)

Computes turbulent surface fluxes using ClimaAtmos's `update_surface_conditions!`. This
requires that we define a new temporary parameter Tuple, `new_p`, and save the new surface state
in it. We do not want `new_p` to live in the atmospheric model permanently, because that would also
trigger flux calculation during Atmos `step!`. We only want to trigger this once per coupling
timestep from ClimaCoupler.
"""
function atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)

    csf.z0m_S .= 1.0e-5
    csf.z0b_S .= 1.0e-5
    csf.beta .= 1
    csf = merge(csf, (;q_sfc = nothing))

    if isnothing(atmos_sim.integrator.p.sfc_setup) # trigger flux calculation if not done in Atmos internally
        new_p = get_new_cache(atmos_sim, csf)
        CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t)
        atmos_sim.integrator.p.sfc_conditions .= new_p.sfc_conditions
    end
end

# function atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)
#     integrator = atmos_sim.integrator
#     params = integrator.p.params

#     sfc_ts = integrator.p.sfc_conditions.ts

#     thermo_params = CAP.thermodynamics_params(params)

#     # set surface ts (T, rho and q)
#     # new_ts = @. surface_ts(csf.ρ_sfc, csf.T_S, csf.q_sfc, thermo_params)
#     # parent(sfc_ts) .= parent(new_ts)

#     # update_surface_conditions_coupler!(integrator.u, integrator.p, integrator.t)
# end

# function surface_ts(ρ_sfc, T_sfc, q_sfc, thermo_params)
#     # Assume that the surface is water with saturated air directly
#     # above it.
#     phase = TD.Liquid()
#     q_vap = TD.q_vap_saturation_generic(thermo_params, T_sfc, ρ_sfc, phase)
#     q = TD.PhasePartition(q_vap)
#     return TD.PhaseNonEquil_ρTq(thermo_params, ρ_sfc, T_sfc, q) # nonequil surface :(
# end


# # # function below allows changes in beta and q_sfc
# function atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)
#     # new_p = get_new_cache(atmos_sim, csf)
#     # CA.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t)
#     # atmos_sim.integrator.p.sfc_conditions .= new_p.sfc_conditions
#     integrator = atmos_sim.integrator
#     params = integrator.p.params

#     ᶜts = integrator.p.ᶜts

#     sfc_ts = integrator.p.sfc_conditions.ts
#     int_ts = Fields.level(integrator.p.ᶜts, 1)

#     sfc_z_values = Fields.field_values(Fields.coordinate_field(sfc_ts).z)
#     int_z_values = Fields.field_values(Fields.coordinate_field(int_ts).z)
#     int_ts_values = Fields.field_values(int_ts)

#     int_u_values = Fields.field_values(Fields.level(integrator.u.c.uₕ, 1))

#     int_local_geometry = Fields.local_geometry_field(Fields.level(ᶜts, 1))
#     int_local_geometry_values = Fields.field_values(int_local_geometry)

#     thermo_params = CAP.thermodynamics_params(params)
#     surface_params = CAP.surface_fluxes_params(params)

#     # set surface ts (T, rho and q)
#     new_ts = @. surface_ts(csf.ρ_sfc, csf.T_S, csf.q_sfc, thermo_params)
#     parent(sfc_ts) .= parent(new_ts)

#     # set surface beta
#     beta = ones(axes(sfc_ts))
#     parent(beta) .= parent(csf.beta) # we can;t pass beta to atmos
#     beta_values = Fields.field_values(beta)

#     sfc_ts_values = Fields.field_values(sfc_ts)

#     parameterization = integrator.p.sfc_setup.parameterization
#     sc_val = @. surface_conditions_coupler(
#         sfc_ts_values,
#         sfc_z_values,
#         int_ts_values,
#         ClimaAtmos.SurfaceConditions.projected_vector_data(CT1, int_u_values, int_local_geometry_values),
#         ClimaAtmos.SurfaceConditions.projected_vector_data(CT2, int_u_values, int_local_geometry_values),
#         int_z_values,
#         surface_params,
#         parameterization.z0m,
#         parameterization.z0b,
#         beta_values,
#         # Fields.field_values(beta),
#     )
#     sc = similar(sfc_ts, SF.SurfaceFluxConditions{FT})
#     parent(sc) .= parent(sc_val)

#     surface_local_geometry = Fields.level(Fields.local_geometry_field(integrator.u.f), Fields.half)
#     @. integrator.p.sfc_conditions =  ClimaAtmos.SurfaceConditions.atmos_surface_conditions(
#         sc,
#         sfc_ts,
#         surface_local_geometry,
#         integrator.p.atmos,
#         params,
#     )

#     # update_surface_conditions_coupler!(integrator.u, integrator.p, integrator.t) # ideally we'd use this to calculate the flues, but this is not possible because we cannot set beta or q_sfc
# end

# function surface_conditions_coupler(
#     sfc_ts_values,
#     sfc_z_values,
#     int_ts_values,
#     interior_u,
#     interior_v,
#     interior_z,
#     surface_params,
#     z0m,
#     z0b,
#     beta,
# )
#     FT = eltype(sfc_ts_values)
#     surface_values =
#         SF.SurfaceValues(sfc_z_values, StaticArrays.SVector(FT(0), FT(0)), sfc_ts_values)
#     interior_values = SF.InteriorValues(
#         interior_z,
#         StaticArrays.SVector(interior_u, interior_v),
#         int_ts_values,
#     )


#     surface_inputs = SF.ValuesOnly(
#         interior_values,
#         surface_values,
#         z0m,#(parameterization.z0m),
#         z0b,# parameterization.z0b),
#         FT(1), # gustiness
#         beta, # beta
#     )
#     return SF.surface_conditions(surface_params, surface_inputs)
# end


# # adapted atmos init cond setter
# using ClimaAtmos
# import ClimaAtmos.SurfaceConditions: update_surface_conditions!

# function update_surface_conditions_coupler!(Y, p, t)
#     # Need to extract the field values so that we can do
#     # a DataLayout broadcast rather than a Field broadcast
#     # because we are mixing surface and interior fields
#     # if isnothing(p.sfc_setup)
#     #     p.is_init[] && set_dummy_surface_conditions!(p)
#     #     return
#     # end

#     set_precomputed_quantities_coupler!(Y, p, t)

#     sfc_local_geometry_values = Fields.field_values(
#         Fields.level(Fields.local_geometry_field(Y.f), Fields.half),
#     )
#     int_local_geometry_values =
#         Fields.field_values(Fields.level(Fields.local_geometry_field(Y.c), 1))
#     (; ᶜts, ᶜu, sfc_conditions, params, sfc_setup, atmos) = p
#     int_ts_values = Fields.field_values(Fields.level(ᶜts, 1))
#     int_u_values = Fields.field_values(Fields.level(ᶜu, 1))
#     int_z_values =
#         Fields.field_values(Fields.level(Fields.coordinate_field(Y.c).z, 1))
#     sfc_conditions_values = Fields.field_values(sfc_conditions)
#     wrapped_sfc_setup = ClimaAtmos.SurfaceConditions.sfc_setup_wrapper(sfc_setup)
#     sfc_temp_var =
#         p.atmos.surface_model isa ClimaAtmos.PrognosticSurfaceTemperature ?
#         (; sfc_prognostic_temp = Fields.field_values(Y.sfc.T)) : (;)
#     @. sfc_conditions_values = ClimaAtmos.SurfaceConditions.surface_state_to_conditions(
#         ClimaAtmos.SurfaceConditions.surface_state(
#             wrapped_sfc_setup,
#             sfc_local_geometry_values,
#             int_z_values,
#             t,
#         ),
#         sfc_local_geometry_values,
#         int_ts_values,
#         ClimaAtmos.SurfaceConditions.projected_vector_data(CT1, int_u_values, int_local_geometry_values),
#         ClimaAtmos.SurfaceConditions.projected_vector_data(CT2, int_u_values, int_local_geometry_values),
#         int_z_values,
#         params,
#         atmos,
#         sfc_temp_var...,
#         )
#     # @info "passed check!! :)"
#     return nothing
# end

# # overwrite this function
# function ClimaAtmos.SurfaceConditions.update_surface_conditions!(Y, p, t)
#     if t ≈ 0.0

#         set_precomputed_quantities_coupler!(Y, p, t)
#         update_surface_conditions_coupler!(Y, p, t)
#         @info "init sfc cond :)"
#         @info t
#     # else
#     #     @info "ignoring atmos"
#     #     @info t
#     end

#     return nothing
# end


# function set_precomputed_quantities_coupler!(Y, p, t)

#     (; energy_form, moisture_model, turbconv_model) = p.atmos
#     thermo_params = CAP.thermodynamics_params(p.params)
#     n = CA.n_mass_flux_subdomains(turbconv_model)
#     thermo_args = (thermo_params, energy_form, moisture_model)
#     (; ᶜspecific, ᶜu, ᶠu³, ᶜK, ᶜts, ᶜp, ᶜΦ) = p
#     ᶠuₕ³ = p.ᶠtemp_CT3

#     @. ᶜspecific = CA.specific_gs(Y.c)
#     CA.set_ᶠuₕ³!(ᶠuₕ³, Y)

#     # TODO: We might want to move this to dss! (and rename dss! to something
#     # like enforce_constraints!).
#     CA.set_velocity_at_surface!(Y, ᶠuₕ³, turbconv_model)

#     CA.set_velocity_quantities!(ᶜu, ᶠu³, ᶜK, Y.f.u₃, Y.c.uₕ, ᶠuₕ³)
#     if n > 0
#         # TODO: In the following increments to ᶜK, we actually need to add
#         # quantities of the form ᶜρaχ⁰ / ᶜρ⁰ and ᶜρaχʲ / ᶜρʲ to ᶜK, rather than
#         # quantities of the form ᶜρaχ⁰ / ᶜρ and ᶜρaχʲ / ᶜρ. However, we cannot
#         # compute ᶜρ⁰ and ᶜρʲ without first computing ᶜts⁰ and ᶜtsʲ, both of
#         # which depend on the value of ᶜp, which in turn depends on ᶜts. Since
#         # ᶜts depends on ᶜK (at least when the energy_form is TotalEnergy), this
#         # means that the amount by which ᶜK needs to be incremented is a
#         # function of ᶜK itself. So, unless we run a nonlinear solver here, this
#         # circular dependency will prevent us from computing the exact value of
#         # ᶜK. For now, we will make the anelastic approximation ᶜρ⁰ ≈ ᶜρʲ ≈ ᶜρ.
#         # add_sgs_ᶜK!(ᶜK, Y, ᶜρa⁰, ᶠu₃⁰, turbconv_model)
#         # @. ᶜK += Y.c.sgs⁰.ρatke / Y.c.ρ
#         # TODO: We should think more about these increments before we use them.
#     end
#     @. ᶜts = CA.ts_gs(thermo_args..., ᶜspecific, ᶜK, ᶜΦ, Y.c.ρ)
#     @. ᶜp = TD.air_pressure(thermo_params, ᶜts)

#     if energy_form isa CA.TotalEnergy
#         (; ᶜh_tot) = p
#         @. ᶜh_tot =
#             TD.total_specific_enthalpy(thermo_params, ᶜts, ᶜspecific.e_tot)
#     end

# end