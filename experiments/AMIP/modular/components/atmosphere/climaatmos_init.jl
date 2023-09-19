# atmos_init: for ClimaAtmos pre-AMIP interface
import ClimaAtmos
using ClimaAtmos: RRTMGPI
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

# the clima atmos `integrator` is now defined
struct ClimaAtmosSimulation{P, Y, D, I} <: AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
name(::ClimaAtmosSimulation) = "ClimaAtmosSimulation"

function atmos_init(::Type{FT}, parsed_args::Dict) where {FT}

    atmos_config = ClimaAtmos.AtmosConfig(argparse_settings(); parsed_args)
    integrator = ClimaAtmos.get_integrator(atmos_config)
    Y = integrator.u
    center_space = axes(Y.c.ρe_tot)
    face_space = axes(Y.f.u₃)
    spaces = (; center_space = center_space, face_space = face_space)
    if :ρe_int in propertynames(Y.c)
        @warn("Running with ρe_int in coupled mode is not tested yet.")
    end

    # set initial fluxes to zero
    @. integrator.p.sfc_conditions.ρ_flux_h_tot = Geometry.Covariant3Vector(FT(0.0))
    @. integrator.p.sfc_conditions.ρ_flux_q_tot = Geometry.Covariant3Vector(FT(0.0))
    @. integrator.p.sfc_conditions.ρ_flux_uₕ.components = zeros(axes(integrator.p.sfc_conditions.ρ_flux_uₕ.components))
    parent(integrator.p.ᶠradiation_flux) .= parent(zeros(axes(integrator.p.ᶠradiation_flux)))
    integrator.p.col_integrated_rain .= FT(0)
    integrator.p.col_integrated_snow .= FT(0)


    ClimaAtmosSimulation(integrator.p.params, Y, spaces, integrator)
end

# extensions required by the Interfacer
get_field(sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux}) =
    Fields.level(sim.integrator.p.ᶠradiation_flux, half)
function get_field(sim::ClimaAtmosSimulation, ::Val{:liquid_precipitation})
    ρ_liq = CAP.ρ_cloud_liq(sim.integrator.p.params)
    (sim.integrator.p.col_integrated_rain .+ sim.integrator.p.col_integrated_snow) .* ρ_liq # kg/m^2/s ; all fallen snow melts for now
end
function get_field(sim::ClimaAtmosSimulation, ::Val{:snow_precipitation}) # kg/m^2/s
    ρ_liq = CAP.ρ_cloud_liq(sim.integrator.p.params)
    sim.integrator.p.col_integrated_snow .* ρ_liq .* FT(0)
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

# extensions required by the Interfacer
function update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_temperature}, csf)
    sim.integrator.p.radiation_model.surface_temperature .= RRTMGPI.field2array(csf.T_S)
end

function update_field!(sim::ClimaAtmosSimulation, ::Val{:albedo}, field)
    sim.integrator.p.radiation_model.diffuse_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(field), 1, length(parent(field)))
    sim.integrator.p.radiation_model.direct_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(field), 1, length(parent(field)))
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

    surface_state(z0m, z0b, T, beta, q_vap) = ClimaAtmos.SurfaceConditions.SurfaceState(;
        parameterization = ClimaAtmos.SurfaceConditions.MoninObukhov(; z0m, z0b),
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
    atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)

Computes turbulent surface fluxes using ClimaAtmos's `update_surface_conditions!`. This
requires that we define a new temporary parameter Tuple, `new_p`, and save the new surface state
in it. We do not want `new_p` to live in the atmospheric model permanently, because that would also
trigger flux calculation during Atmos `step!`. We only want to trigger this once per coupling
timestep from ClimaCoupler.
"""
function atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)
    new_p = get_new_cache(atmos_sim, csf)
    ClimaAtmos.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t)
    atmos_sim.integrator.p.sfc_conditions .= new_p.sfc_conditions
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
