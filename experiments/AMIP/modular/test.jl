# test calculate_and_send_turbulent_fluxes!(cs)
import CLIMAParameters as CP
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import SurfaceFluxes as SF

using ClimaCore: Fields
using ClimaCoupler: Utilities, TestHelper
using Test
using StaticArrays

FT = Float64
include("components/flux_calculator.jl")

# Parameter set generated in ClimaAtmos GCM run (copied from SurfaceFluxes.jl)
const sf_params = SF.Parameters.SurfaceFluxesParameters{
    FT,
    SF.UniversalFunctions.BusingerParams{FT},
    TP.ThermodynamicsParameters{FT},
}(
    0.4f0,
    SF.UniversalFunctions.BusingerParams{FT}(0.74f0, 4.7f0, 4.7f0, 2.5f0, 4.45f0),
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

# atmos sim
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
get_height_int(sim::TestAtmos) = sim.integrator.p.z
get_height_int_point(sim::TestAtmos, colidx) = get_height_int(sim)[colidx]
get_uv_int_point(sim::TestAtmos, colidx) = @. StaticArrays.SVector(sim.integrator.p.u, sim.integrator.p.v)[colidx]
get_thermo_state(sim::TestAtmos) = TD.PhaseEquil_ρTq.(get_thermo_params(sim), sim.integrator.ρ, sim.integrator.T, sim.integrator.q)
get_thermo_state_point(sim::TestAtmos, colidx) = get_thermo_state(sim)[colidx]
get_height_sfc(sim::TestAtmos) = sim.integrator.p.z_sfc
get_height_sfc_point(sim::TestAtmos, colidx) = get_height_sfc(sim)[colidx]
get_air_density(::TestAtmos, thermo_params, thermo_state) = TD.air_density.(thermo_params, thermo_state)
get_air_temperature(::TestAtmos, thermo_params, thermo_state_int) = TD.air_temperature.(thermo_params, thermo_state_int)
get_cv_m(::TestAtmos, thermo_params, thermo_state_int) = TD.cv_m.(thermo_params, thermo_state_int)
get_gas_constant_air(::TestAtmos, thermo_params, thermo_state_int)  = TD.gas_constant_air.(thermo_params, thermo_state_int)
get_q_vap_saturation_generic(::TestAtmos, thermo_params, thermo_state_int) = TD.q_vap_saturation_generic.(thermo_params, T_sfc, ρ_sfc, TD.Liquid())
function update_calculated_fluxes_point!(sim::TestAtmos, fields, colidx)
    (; F_ρτxz, F_shf, F_lhf, F_evap ) = fields

    ρ_int = sim.integrator.ρ
    @. sim.integrator.p.energy_bc[colidx] = - F_shf[colidx] + F_lhf[colidx] # Geometry.WVector(outputs[colidx].F_shf + outputs[colidx].F_lhf,)
    @. sim.integrator.p.ρq_tot_bc[colidx] = - F_evap[colidx]
    @. sim.integrator.p.uₕ_bc[colidx] = - (F_ρτxz / ρ_int)[colidx] # x-compoennt only

end

# ocean sim
struct TestOcean{M, Y, D, I} <: SurfaceModelSimulation
    model::M
    Y_init::Y
    domain::D
    integrator::I
end
get_temperature(sim::TestOcean) = sim.integrator.T
get_temperature_point(sim::TestOcean, colidx) = get_temperature(sim)[colidx]
get_humidity(sim::TestOcean) = sim.integrator.p.q
get_humidity_point(sim::TestOcean, colidx) = get_humidity(sim)[colidx]
get_z0m_point(sim::TestOcean, colidx) = sim.integrator.p.z0m
get_z0b_point(sim::TestOcean, colidx)= sim.integrator.p.z0b
get_beta_point(sim::TestOcean, colidx) = sim.integrator.p.beta[colidx]
get_area_fraction(sim::TestOcean) = ones(axes(sim.integrator.p.F_aero))
get_heat_transfer_coefficient_point(sim::TestOcean, colidx) = sim.integrator.p.Ch
get_drag_transfer_coefficient_point(sim::TestOcean, colidx) = sim.integrator.p.Cd
function update_calculated_fluxes_point!(sim::TestOcean, fields, colidx)
    (; F_shf, F_lhf) = fields
    @. sim.integrator.p.F_aero[colidx] = F_shf + F_lhf
end
issaturated(::TestOcean, q) = isnothing(q)


@testset "calculate correct fluxes" begin
    surface_scheme_list = (MoninObukhovScheme(), BulkScheme())
    for scheme in surface_scheme_list
        boundary_space = TestHelper.create_space(FT)

        params = (; surface_scheme = scheme, FT = FT)

        p = (; energy_bc = zeros(boundary_space), ρq_tot_bc = zeros(boundary_space), uₕ_bc = ones(boundary_space), z = ones(boundary_space), z_sfc = zeros(boundary_space), u = ones(boundary_space), v = ones(boundary_space))
        Y_init = (; ρ = ones(boundary_space) .* 1.2, T = ones(boundary_space) .* 310, q = zeros(boundary_space))
        integrator = (; Y_init... , p = p, )
        atmos_sim = TestAtmos(params, Y_init, nothing, integrator)

        p = (; F_aero = zeros(boundary_space), z0m = FT(0.01), z0b = FT(0.01), beta = ones(boundary_space), q = zeros(boundary_space), Cd = FT(0.01), Ch = FT(0.01))
        Y_init = (; T = ones(boundary_space) .* 300.0,)
        integrator = (; Y_init... , p = p, )
        ocean_sim = TestOcean(nothing, Y_init, nothing, integrator)

        model_sims = (; atmos_sim, ocean_sim);

        coupler_cache_names = (:T_S, :albedo, :F_R_sfc, :F_R_toa, :P_liq, :P_snow, :P_net, :F_lhf, :F_shf, :F_ρτxz, :F_ρτyz,  :F_evap)
        fields = NamedTuple{coupler_cache_names}(ntuple(i -> Fields.zeros(boundary_space), length(coupler_cache_names)))

        # @test
        calculate_and_send_turbulent_fluxes!(model_sims, fields, boundary_space)

        windspeed =  @. hypot(atmos_sim.integrator.p.u, atmos_sim.integrator.p.v)


        thermo_params = get_thermo_params(atmos_sim)
        thermo_state_int = get_thermo_state(atmos_sim)
        thermo_state_sfc = surface_thermo_state(
                            atmos_sim,
                            ocean_sim,
                            thermo_params,
                            get_temperature(ocean_sim),#ocean_sim.integrator.T,
                            get_humidity(ocean_sim),
                            thermo_state_int,
                        )
        ρ_sfc = get_air_density(atmos_sim, thermo_params, thermo_state_sfc)

        cpm = get_cv_m(atmos_sim, thermo_params, thermo_state_int) .+ get_gas_constant_air(atmos_sim, thermo_params, thermo_state_int) # cp = R + cv
        gz = (get_height_int(atmos_sim) .- get_height_sfc(atmos_sim)) .* 9.81
        shf_analytical = @. (cpm * (ocean_sim.integrator.T - atmos_sim.integrator.T) - gz) * ocean_sim.integrator.p.Ch * ρ_sfc * windspeed #-ρ_sfc * Ch * windspeed(sc) * (cp_m * ΔT + ΔΦ)

        if scheme == BulkScheme()
            @test isapprox(parent(shf_analytical)[1], parent(fields.F_shf)[1], rtol = 1e-6)
        end
        @test parent(fields.F_evap)[1] ≈ FT(0)
        @test parent(fields.F_lhf)[1] ≈ FT(0)
    end

    # TODO: add test for the moist case

end

