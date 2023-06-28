# atmos_init: for ClimaAtmos pre-AMIP interface
import ClimaAtmos
using ClimaAtmos: RRTMGPI
import ClimaCoupler.FluxCalculator: compute_atmos_turbulent_fluxes!, calculate_surface_air_density
using ClimaCore: Fields.level, Geometry
import ClimaCoupler.FieldExchanger: get_thermo_params

# the clima atmos `integrator` is now defined
struct ClimaAtmosSimulation{P, Y, D, I} <: AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
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

    ClimaAtmosSimulation(integrator.p.params, Y, spaces, integrator)
end

# required by the Interfacer
get_field(sim::ClimaAtmosSimulation, ::Val{:radiative_energy_flux}) =
    Fields.level(sim.integrator.p.ᶠradiation_flux, half)
get_field(sim::ClimaAtmosSimulation, ::Val{:liquid_precipitation}) =
    sim.integrator.p.col_integrated_rain .+ sim.integrator.p.col_integrated_snow # all fallen snow melts for now
get_field(sim::ClimaAtmosSimulation, ::Val{:snow_precipitation}) = sim.integrator.p.col_integrated_snow .* FT(0)

get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_energy_flux}) =
    Geometry.WVector.(sim.integrator.p.sfc_conditions.ρ_flux_h_tot)
get_field(sim::ClimaAtmosSimulation, ::Val{:turbulent_moisture_flux}) =
    Geometry.WVector.(sim.integrator.p.sfc_conditions.ρ_flux_q_tot)

function update_field!(sim::ClimaAtmosSimulation, ::Val{:surface_temperature}, field)
    sim.integrator.p.radiation_model.surface_temperature .= RRTMGPI.field2array(field)
end
function update_field!(sim::ClimaAtmosSimulation, ::Val{:albedo}, field)
    sim.integrator.p.radiation_model.diffuse_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(field), 1, length(parent(field)))
    sim.integrator.p.radiation_model.direct_sw_surface_albedo .=
        reshape(RRTMGPI.field2array(field), 1, length(parent(field)))
end

step!(sim::ClimaAtmosSimulation, t) = step!(sim.integrator, t - sim.integrator.t, true)

reinit!(sim::ClimaAtmosSimulation) = reinit!(sim.integrator)

"""
    update_sim!(atmos_sim::ClimaAtmosSimulation, csf)

Updates the surface fields for temperature and albedo.

# Arguments
- `atmos_sim`: [ClimaAtmosSimulation] containing an atmospheric model simulation object.
- `csf`: [NamedTuple] containing coupler fields.
"""
function update_sim!(atmos_sim::ClimaAtmosSimulation, csf, turbulent_fluxes)
    update_field!(atmos_sim, Val(:albedo), csf.albedo)
    update_field!(atmos_sim, Val(:surface_temperature), csf.T_S)

    if turbulent_fluxes !== FluxCalculator.CombinedAtmosGrid()
        Fields.bycolumn(axes(csf.T_S)) do colidx
            coupler_fields = (; F_ρτxz = csf.F_ρτxz[colidx], F_ρτyz = csf.F_ρτyz[colidx], F_shf = csf.F_shf[colidx], F_lhf = csf.F_lhf[colidx], F_evap = csf.F_evap[colidx])
            update_turbulent_fluxes_point!(atmos_sim, coupler_fields, colidx)
        end
    end
end

function update_turbulent_fluxes_point!(sim::ClimaAtmosSimulation, fields, colidx)
    (; F_ρτxz , F_ρτyz , F_shf , F_lhf , F_evap ) = fields

    ρ_int = Spaces.level(sim.integrator.u.c.ρ , 1)

    @. sim.integrator.p.dif_flux_energy_bc[colidx] = - Geometry.WVector(F_shf[colidx] + F_lhf[colidx]) # Geometry.WVector(outputs[colidx].F_shf + outputs[colidx].F_lhf,)
    @. sim.integrator.p.dif_flux_ρq_tot_bc[colidx] = - Geometry.WVector(F_evap[colidx])
    @. sim.integrator.p.dif_flux_uₕ_bc[colidx] = - Geometry.Contravariant3Vector(sim.integrator.p.surface_normal[colidx]) ⊗ Geometry.Covariant12Vector(Geometry.UVVector(F_ρτxz / ρ_int, F_ρτyz / ρ_int)[colidx],)

end


# flux calculation borrowed from atmos
"""
    CoupledMoninObukhov()
A modified version of a Monin-Obukhov surface for the Coupler, see the link below for more information
https://clima.github.io/SurfaceFluxes.jl/dev/SurfaceFluxes/#Monin-Obukhov-Similarity-Theory-(MOST)
"""
struct CoupledMoninObukhov end
"""
    coupler_surface_setup(::CoupledMoninObukhov, p, csf_sfc = (; T = nothing, z0m = nothing, z0b = nothing, beta = nothing))

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
    compute_atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)

Computes turbulent surface fluxes using ClimaAtmos's `update_surface_conditions!`. This
requires that we define a new temporary parameter Tuple, `new_p`, and save the new surface state
in it. We do not want `new_p` to live in the atmospheric model permanently, because that would also
trigger flux calculation during Atmos `step!`. We only want to trigger this once per coupling
timestep from ClimaCoupler.
"""
function compute_atmos_turbulent_fluxes!(atmos_sim::ClimaAtmosSimulation, csf)

    p = atmos_sim.integrator.p

    thermo_params = get_thermo_params(atmos_sim)
    ts_int = get_field(atmos_sim, Val(:thermo_state_int))
    FT = eltype(csf.T_S)

    # surface density is needed for q_sat and requires atmos and sfc states, so it is calculated and saved in the coupler
    # parent(csf.ρ_sfc) .=
    #     parent(extrapolate_ρ_to_sfc.(Ref(thermo_params), ts_int, swap_space!(ones(axes(ts_int.ρ)), csf.T_S)))

    # # surface specific q_vap is calculated as a bulk quantity here, assuming a saturated surface.
    # # TODO - use land's q_sfc, but ClimaLSM need to be modified to ingest model type
    # # NB: q_sfc (q_sfc is actually saturated!) is a passive variable in Land and nonexistent in the slabs, so this is ok for testing,
    # # though this q_sfc does not account for `q_vap_saturation_generic(..., TD.Ice())`. To approximately account for undersaturation of the surface,
    # # we can use beta to adjust evaporation to appoximate undersaturation.
    # q_vap_sfc = TD.q_vap_saturation_generic.(thermo_params, csf.T_S, csf.ρ_sfc, TD.Liquid())

    coupler_sfc_setup = coupler_surface_setup(
        CoupledMoninObukhov(),
        p;
        csf_sfc = (; T = csf.T_S, z0m = csf.z0m_S, z0b = csf.z0b_S, beta = csf.beta, q_vap = csf.q_sfc .* FT(1.0)),
    )

    p_names = propertynames(p)
    p_values = map(x -> x == :sfc_setup ? coupler_sfc_setup : getproperty(p, x), p_names)

    new_p = (; zip(p_names, p_values)...)

    ClimaAtmos.SurfaceConditions.update_surface_conditions!(atmos_sim.integrator.u, new_p, atmos_sim.integrator.t)

    p.sfc_conditions .= new_p.sfc_conditions

end

get_thermo_params(sim::ClimaAtmosSimulation) = CAP.thermodynamics_params(sim.integrator.p.params)
get_field(sim::ClimaAtmosSimulation, ::Val{:thermo_state_int}) = Spaces.level(sim.integrator.p.ᶜts, 1)
get_field(sim::ClimaAtmosSimulation, ::Val{:height_int})  = Spaces.level(Fields.coordinate_field(sim.integrator.u.c).z, 1)
get_field(sim::ClimaAtmosSimulation, ::Val{:height_sfc})  = Spaces.level(Fields.coordinate_field(sim.integrator.u.f).z, half)
function get_field(sim::ClimaAtmosSimulation, ::Val{:uv_int})
    uₕ_int = Geometry.UVVector.(Spaces.level(sim.integrator.u.c.uₕ, 1))
    return @. StaticArrays.SVector(uₕ_int.components.data.:1, uₕ_int.components.data.:2)
end


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
    extrapolate_ρ_to_sfc(thermo_params, ts_int, T_sfc)

Uses the ideal gas law and hydrostatic balance to extrapolate for surface density pointwise.
"""
function extrapolate_ρ_to_sfc(thermo_params, ts_in, T_sfc)
    T_int = TD.air_temperature(thermo_params, ts_in)
    Rm_int = TD.gas_constant_air(thermo_params, ts_in)
    ρ_air = TD.air_density(thermo_params, ts_in)
    ρ_air * (T_sfc / T_int)^(TD.cv_m(thermo_params, ts_in) / Rm_int)
end

"""
    update_turbulent_fluxes_point!(sim::ClimaAtmosSimulation, fields, colidx)

Updates the turbulent fluxes in the `integrator` of `sim` using the fields in `fields` at the column index `colidx`.
"""
function update_turbulent_fluxes_point!(sim::ClimaAtmosSimulation, fields, colidx)
    (; F_ρτxz , F_ρτyz , F_shf , F_lhf , F_evap ) = fields

    ρ_int = Spaces.level(sim.integrator.u.c.ρ , 1)

    @. sim.integrator.p.dif_flux_energy_bc[colidx] = - Geometry.WVector(F_shf[colidx] + F_lhf[colidx]) # Geometry.WVector(outputs[colidx].F_shf + outputs[colidx].F_lhf,)
    @. sim.integrator.p.dif_flux_ρq_tot_bc[colidx] = - Geometry.WVector(F_evap[colidx])
    @. sim.integrator.p.dif_flux_uₕ_bc[colidx] = - Geometry.Contravariant3Vector(sim.integrator.p.surface_normal[colidx]) ⊗ Geometry.Covariant12Vector(Geometry.UVVector(F_ρτxz / ρ_int, F_ρτyz / ρ_int)[colidx],)

end