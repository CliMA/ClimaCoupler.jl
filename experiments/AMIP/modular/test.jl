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
get_surface_params(::TestAtmos) = sf_params

function get_thermo_params(sim::TestAtmos)
    FT = sim.params.FT
    aliases = string.(fieldnames(TP.ThermodynamicsParameters))
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    TP.ThermodynamicsParameters{FT}(; pairs...)
end

get_field(sim::TestAtmos, ::Val{:height_int}) = sim.integrator.p.z
get_field(sim::TestAtmos, ::Val{:height_sfc}) = sim.integrator.p.z_sfc
get_field(sim::TestAtmos, ::Val{:uv_int_point}) = @. StaticArrays.SVector(sim.integrator.p.u, sim.integrator.p.v)
get_field(sim::TestAtmos, ::Val{:thermo_state_int}) = TD.PhaseEquil_ρTq.(get_thermo_params(sim), sim.integrator.ρ, sim.integrator.T, sim.integrator.q)
get_field(sim::TestAtmos, ::Val{:air_density}) = sim.integrator.ρ
get_field(sim::TestAtmos, ::Val{:air_temperature}) = sim.integrator.T

# get_q_vap_saturation_generic(::TestAtmos, thermo_params, thermo_state_int) = TD.q_vap_saturation_generic.(thermo_params, T_sfc, ρ_sfc, TD.Liquid())
function update_turbulent_fluxes_point!(sim::TestAtmos, fields, colidx)
    (; F_ρτxz, F_shf, F_lhf, F_evap ) = fields
    ρ_int = sim.integrator.ρ[colidx]
    @. sim.integrator.p.energy_bc[colidx] = - (F_shf + F_lhf)
    @. sim.integrator.p.ρq_tot_bc[colidx] = - F_evap
    @. sim.integrator.p.uₕ_bc[colidx] = - (F_ρτxz / ρ_int) # x-compoennt only for this test
end

# ocean sim
struct TestOcean{M, Y, D, I} <: SurfaceModelSimulation
    model::M
    Y_init::Y
    domain::D
    integrator::I
end
get_field(sim::TestOcean, ::Val{:air_temperature}) = sim.integrator.T
get_field(sim::TestOcean, ::Val{:air_humidity}) = sim.integrator.p.q
get_field(sim::TestOcean, ::Val{:z0m}) = sim.integrator.p.z0m
get_field(sim::TestOcean, ::Val{:z0b}) = sim.integrator.p.z0b
get_field(sim::TestOcean, ::Val{:beta}) = sim.integrator.p.beta
get_field(sim::TestOcean, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
get_field(sim::TestOcean, ::Val{:area_mask}, colidx) = Regridder.binary_mask.(get_field(sim, Val{:area_fraction})[colidx], threshold = eps())

get_field(sim::TestOcean, ::Val{:heat_transfer_coefficient}) = sim.integrator.p.Ch
get_field(sim::TestOcean, ::Val{:drag_coefficient}) = sim.integrator.p.Cd

# overriding the default here
function surface_thermo_state(sim::TestOcean, thermo_params, thermo_state_int, colidx)
    T_sfc = get_field(sim, Val(:air_temperature), colidx)
    ρ_sfc = thermo_state_int.ρ # arbitrary
    q_sfc = get_field(sim, Val(:air_humidity), colidx) # read from cache
    @. TD.PhaseEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, q_sfc)
end


function update_turbulent_fluxes_point!(sim::TestOcean, fields, colidx)
    (; F_shf, F_lhf) = fields
    @. sim.integrator.p.F_aero[colidx] = F_shf + F_lhf
end

# issaturated(::TestOcean, q) = isnothing(q)

@testset "calculate correct fluxes: dry" begin
    surface_scheme_list = (MoninObukhovScheme(), BulkScheme())
    for scheme in surface_scheme_list
        boundary_space = TestHelper.create_space(FT);

        params = (; surface_scheme = scheme, FT = FT);

        # atmos
        p = (; energy_bc = zeros(boundary_space), ρq_tot_bc = zeros(boundary_space), uₕ_bc = ones(boundary_space), z = ones(boundary_space), z_sfc = zeros(boundary_space), u = ones(boundary_space), v = ones(boundary_space));
        Y_init = (; ρ = ones(boundary_space) .* 1.2, T = ones(boundary_space) .* 310, q = zeros(boundary_space));
        integrator = (; Y_init... , p = p, );
        atmos_sim = TestAtmos(params, Y_init, nothing, integrator);

        # ocean
        p = (; F_aero = zeros(boundary_space), z0m = FT(0.01), z0b = FT(0.01), beta = ones(boundary_space), q = zeros(boundary_space), Cd = FT(0.01), Ch = FT(0.01), area_fraction = ones(boundary_space));
        Y_init = (; T = ones(boundary_space) .* 300.0,);
        integrator = (; Y_init... , p = p, );
        ocean_sim = TestOcean(nothing, Y_init, nothing, integrator);

        model_sims = (; atmos_sim, ocean_sim);

        coupler_cache_names = (:T_S, :albedo, :F_R_sfc, :F_R_toa, :P_liq, :P_snow, :P_net, :F_lhf, :F_shf, :F_ρτxz, :F_ρτyz,  :F_evap)
        fields = NamedTuple{coupler_cache_names}(ntuple(i -> Fields.zeros(boundary_space), length(coupler_cache_names)));

        # calculate turbulent fluxes
        thermo_params = get_thermo_params(atmos_sim)
        calculate_and_send_turbulent_fluxes!(model_sims, fields, boundary_space, scheme, thermo_params)

        # calculating the fluxes twice ensures that no accumulation occurred
        calculate_and_send_turbulent_fluxes!(model_sims, fields, boundary_space, scheme, thermo_params)

        windspeed =  @. hypot(atmos_sim.integrator.p.u, atmos_sim.integrator.p.v)

        thermo_params = get_thermo_params(atmos_sim)
        thermo_state_int = get_field(atmos_sim, Val(:thermo_state_int))

        surface_thermo_states = similar(thermo_state_int)
        Fields.bycolumn(boundary_space) do colidx
            surface_thermo_states[colidx] .= surface_thermo_state(
                                ocean_sim,
                                thermo_params,
                                thermo_state_int[colidx],
                                colidx)
        end

        if scheme == BulkScheme()
            ρ_sfc = get_field(atmos_sim, Val(:air_density))
            cpm = TD.cv_m.(thermo_params, thermo_state_int) .+ TD.gas_constant_air.(thermo_params, thermo_state_int) # cp = R + cv
            gz = (get_field(atmos_sim, Val(:height_int))  .- get_field(atmos_sim, Val(:height_sfc)) ) .* 9.81
            shf_analytical = @. (cpm * (ocean_sim.integrator.T - atmos_sim.integrator.T) - gz) * ocean_sim.integrator.p.Ch * ρ_sfc * windspeed #-ρ_sfc * Ch * windspeed(sc) * (cp_m * ΔT + ΔΦ)

            colidx = Fields.ColumnIndex{2}((1, 1), 73)
            @test isapprox(parent(shf_analytical[colidx]), parent(fields.F_shf[colidx]), rtol = 1e-6)
        end
        @test parent(fields.F_evap)[1] ≈ FT(0)
        @test parent(fields.F_lhf)[1] ≈ FT(0)
    end

    # TODO: add test for the moist case

end

# more granular unit tests
struct DummySimulation <: SurfaceModelSimulation end
@testset "extra_aux_update" begin
    sim = DummySimulation()
    @test extra_aux_update(sim) == nothing
end

@testset "get_field" begin
    sim = DummySimulation()
    boundary_space = TestHelper.create_space(FT);
    colidx = Fields.ColumnIndex{2}((1, 1), 73)
    get_field(::DummySimulation, ::Val{:var}) = ones(boundary_space)
    @test parent(get_field(sim, Val(:var0), colidx))[1] == FT(1)

end

@testset "reinit!" begin
    @test reinit!(:sim) isa Warning
end

@testset "update!" begin
    @test update!(:sim, Val(:sim)) isa Warning
end

@testset "name" begin
    @test name(:sim) == "stub  simulation"
end

struct DummySurfaceSimulation <: SurfaceModelSimulation end
dummy_sim = DummySurfaceSimulation()
get_field(::DummySurfaceSimulation, ::Val{:air_temperature}) = ones(boundary_space) .* FT(300)
@testset "surface_thermo_state" begin

    options = (; T = [FT(300), FT(300), FT(310)], q = [FT(0.001), FT(0.01), FT(0.01)], )
    results = (; q_sfc =[FT(0.02550593313501068), FT(0.02550593313501068), FT(0.02769252063578527)], ρ_sfc = [FT(1), FT(1), FT(0.9210405029743302)] )

    # dummy surface

    for i in [1,2]

        T = options.T[i]
        q = options.q[i]

        boundary_space = TestHelper.create_space(FT);

        params = (; surface_scheme = scheme, FT = FT);
        colidx = Fields.ColumnIndex{2}((1, 1), 73)

        # atmos


        p = (; energy_bc = zeros(boundary_space), ρq_tot_bc = zeros(boundary_space), uₕ_bc = ones(boundary_space), z = ones(boundary_space), z_sfc = zeros(boundary_space), u = ones(boundary_space), v = ones(boundary_space));
        Y_init = (; ρ = ones(boundary_space) , T = ones(boundary_space) .* T, q = ones(boundary_space) .* q );
        integrator = (; Y_init... , p = p, );
        atmos_sim = TestAtmos(params, Y_init, nothing, integrator);
        thermo_state_int = get_field(atmos_sim, Val(:thermo_state_int))[colidx]
        thermo_params = get_thermo_params(atmos_sim)

        # default thermo state from atmos
        surface_state = surface_thermo_state(dummy_sim, thermo_params, thermo_state_int, colidx)
        @test parent(surface_state.q_tot)[1] ≈ results.q_sfc[i]
        @test parent(surface_state.ρ)[1] ≈ results.ρ_sfc[i]
    end
end



# tests for combined fluxes

#=
collect_atmos_fluxes!
push_atmos_fluxes!
collect_surface_state!
push_surface_state!


extra_aux_update(sim, thermo_params, get_thermo_state(atmos_sim)) == nothing


=#


