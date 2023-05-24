# test calculate_and_send_turbulent_fluxes!(cs)
import CLIMAParameters as CP
import Thermodynamics as TD
import Thermodynamics.Parameters as TP

using ClimaCore: Fields
using ClimaCoupler: Utilities, TestHelper
using Test
using SurfaceFluxes

FT = Float64
include("components/flux_calculator.jl")

# Parameter set generated in ClimaAtmos GCM run (copied from SurfaceFluxes.jl)
const sf_params = SurfaceFluxes.Parameters.SurfaceFluxesParameters{
    FT,
    SurfaceFluxes.UniversalFunctions.BusingerParams{FT},
    TP.ThermodynamicsParameters{FT},
}(
    0.4f0,
    SurfaceFluxes.UniversalFunctions.BusingerParams{FT}(0.74f0, 4.7f0, 4.7f0, 2.5f0, 4.45f0),
    TP.ThermodynamicsParameters{FT}(
        273.16f0,
        100000.0f0,
        1859.0f0,
        4181.0f0,
        2100.0f0,
        2.5008f6,
        2.8344f6,
        611.657f0,
        273.16f0,
        273.15f0,
        150.0f0,
        1000.0f0,
        298.15f0,
        6864.8f0,
        10513.6f0,
        0.2857143f0,
        8.31446f0,
        0.02897f0,
        0.01801528f0,
        290.0f0,
        220.0f0,
        9.80616f0,
        233.0f0,
        1.0f0,
    ),
)

struct TestAtmos{P, Y, D, I} <: AtmosModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
get_surface_scheme(sim::TestAtmos) = sim.params.surface_scheme
get_surface_params(::TestAtmos) = sf_params
function get_thermo_params(sim::TestAtmos)
    FT = sim.params.FT
    aliases = string.(fieldnames(TP.ThermodynamicsParameters))
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    TP.ThermodynamicsParameters{FT}(; pairs...)
end
get_height_int_point(sim::TestAtmos, colidx) = 1
get_uv_int_point(sim::TestAtmos, colidx) = 1
get_thermo_state_point(sim::TestAtmos, colidx) = TD.PhaseEquil_ρTq.(get_thermo_params(sim), sim.integrator.ρ, sim.integrator.T, sim.integrator.q)[colidx]
get_height_sfc_point(sim::TestAtmos, colidx) = 0
get_air_density(::TestAtmos, thermo_params, thermo_state) = TD.air_density.(thermo_params, thermo_state)
get_air_temperature(::TestAtmos, thermo_params, thermo_state_int) = TD.air_temperature.(thermo_params, thermo_state_int)
get_cv_m(::TestAtmos, thermo_params, thermo_state_int) = TD.cv_m.(thermo_params, thermo_state_int)
get_gas_constant_air(::TestAtmos, thermo_params, thermo_state_int)  = TD.gas_constant_air.(thermo_params, thermo_state_int)
get_q_vap_saturation_generic(::TestAtmos, thermo_params, thermo_state_int) = TD.q_vap_saturation_generic.(thermo_params, T_sfc, ρ_sfc, TD.Liquid())


function update_calculated_fluxes_point!(sim::TestAtmos, fields, colidx)
    (; F_ρτxz , F_ρτyz , F_shf , F_lhf , F_evap ) = fields

    ρ_int = Spaces.level(sim.integrator.u.c.ρ , 1)

    @. sim.integrator.p.energy_bc[colidx] = - F_shf[colidx] + F_lhf[colidx] # Geometry.WVector(outputs[colidx].F_shf + outputs[colidx].F_lhf,)
    @. sim.integrator.p.ρq_tot_bc[colidx] = - F_evap[colidx]
    @. sim.integrator.p.uₕ_bc[colidx] = - (F_ρτxz / ρ_int)[colidx] # x-compoennt only

end

struct TestOcean{M, Y, D, I} <: SurfaceModelSimulation
    model::M
    Y_init::Y
    domain::D
    integrator::I
end
get_temperature_point(sim::TestOcean, colidx) = 300.0
get_humidity_point(sim::TestOcean, colidx) = 0.01
get_z0m_point(sim::TestOcean, colidx) = 0.01
get_z0b_point(sim::TestOcean, colidx)= 0.01
get_beta_point(sim::TestOcean, colidx) = 1.0
get_area_fraction(sim::TestOcean) = ones(axes(sim.integrator.p.F_aero))
function update_calculated_fluxes_point!(sim::TestOcean, fields, colidx)
    (; F_shf, F_lhf) = fields
    @. sim.integrator.p.F_aero[colidx] = F_shf + F_lhf
end
issaturated(::TestOcean, q) = isnothing(q)

boundary_space = TestHelper.create_space(FT)

params = (; surface_scheme = MoninObukhovScheme(), FT = FT)
p = (; energy_bc = ones(boundary_space), ρq_tot_bc = ones(boundary_space), uₕ_bc = ones(boundary_space))
Y_init = (; ρ = ones(boundary_space), T = ones(boundary_space), q = ones(boundary_space))
integrator = (; Y_init... , p = p, )
atmos_sim = TestAtmos(params, Y_init, nothing, integrator)

p = (; F_aero = ones(boundary_space),)
Y_init = (; ρ = ones(boundary_space), T = ones(boundary_space), q = ones(boundary_space))
integrator = (; Y_init... , p = p, )
ocean_sim = TestOcean(nothing, Y_init, nothing, integrator)

sims = (; atmos_sim, ocean_sim)

coupler_cache_names = (:T_S, :albedo, :F_R_sfc, :F_R_toa, :P_liq, :P_snow, :P_net, :F_lhf, :F_shf, :F_ρτxz, :F_ρτyz,  :F_evap)
coupler_fields = NamedTuple{coupler_cache_names}(ntuple(i -> Fields.zeros(boundary_space), length(coupler_cache_names)))

# @test
calculate_and_send_turbulent_fluxes!(sims, coupler_fields, boundary_space)






