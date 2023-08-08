
using ClimaCore: Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaCoupler: Utilities, Regridder, TestHelper
using Test
import ClimaCoupler.FluxCalculator:
    atmos_turbulent_fluxes!,
    combined_turbulent_fluxes!,
    CombinedStateFluxes,
    PartitionedStateFluxes,
    calculate_surface_air_density,
    MoninObukhovScheme,
    BulkScheme,
    partitioned_turbulent_fluxes!,
    get_surface_params,
    update_turbulent_fluxes_point!,
    surface_thermo_state
import ClimaCoupler: Interfacer
FT = Float64

import CLIMAParameters as CP
import Thermodynamics as TD
import Thermodynamics.Parameters as TP
import SurfaceFluxes as SF
import SurfaceFluxes.UniversalFunctions as UF

using StaticArrays

# # test for a simple generic atmos model
struct DummySimulation{S, C} <: Interfacer.AtmosModelSimulation
    state::S
    cache::C
end
function atmos_turbulent_fluxes!(sim::DummySimulation, csf)
    sim.cache.flux .= (csf.T_sfc .- sim.state.T) .* sim.cache.κ ./ sim.cache.dz # Eq. 1
end
struct DummySimulation2{C} <: Interfacer.AtmosModelSimulation
    cache::C
end

@testset "combined_turbulent_fluxes!" begin
    boundary_space = TestHelper.create_space(FT)
    coupler_fields = (; T_sfc = 310 .* ones(boundary_space))
    sim =
        DummySimulation((; T = 300 .* ones(boundary_space)), (; κ = FT(0.01), dz = FT(1), flux = zeros(boundary_space)))
    model_sims = (; atmos_sim = sim)
    flux_types = (CombinedStateFluxes(), PartitionedStateFluxes())
    # the result of Eq 1 above, given these states, is 0.1 W/m2, but under PartitionedStateFluxes() turbulent fluxes are
    # not calculated using this method (using combined surface properties), so the fluxes remain 0.
    results = [FT(0.1), FT(0.0)]
    for (i, t) in enumerate(flux_types)
        sim.cache.flux .= FT(0)
        combined_turbulent_fluxes!(model_sims, coupler_fields, t)
        @test parent(sim.cache.flux)[1] == results[i]
    end
    sim2 = DummySimulation2((; cache = (; flux = zeros(boundary_space))))
    model_sims = (; atmos_sim = sim2)
    @test_throws ErrorException atmos_turbulent_fluxes!(sim2, coupler_fields)

end

@testset "calculate_surface_air_density" begin
    boundary_space = TestHelper.create_space(FT)
    coupler_fields = (; T_sfc = 310 .* ones(boundary_space))
    sim2 = DummySimulation2((; cache = (; flux = zeros(boundary_space))))
    @test_throws ErrorException calculate_surface_air_density(sim2, coupler_fields.T_sfc)
end

# Test for the PartitionedStateFluxes() case
# parameter set generated CLIMAParamerers and helper functions from SurfaceFluxes.jl
function create_uf_parameters(toml_dict, ::UF.BusingerType)
    FT = CP.float_type(toml_dict)
    aliases = ["Pr_0_Businger", "a_m_Businger", "a_h_Businger", "ζ_a_Businger", "γ_Businger"]

    pairs = CP.get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple

    pairs = (;
        Pr_0 = pairs.Pr_0_Businger,
        a_m = pairs.a_m_Businger,
        a_h = pairs.a_h_Businger,
        ζ_a = pairs.ζ_a_Businger,
        γ = pairs.γ_Businger,
    )
    return UF.BusingerParams{FT}(; pairs...)
end

function create_parameters(toml_dict, ufpt)
    FT = CP.float_type(toml_dict)

    ufp = create_uf_parameters(toml_dict, ufpt)
    AUFP = typeof(ufp)

    aliases = string.(fieldnames(TD.Parameters.ThermodynamicsParameters))
    pairs = CP.get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    thermo_params = TD.Parameters.ThermodynamicsParameters{FT}(; pairs...)
    TP = typeof(thermo_params)

    aliases = ["von_karman_const"]
    pairs = CP.get_parameter_values!(toml_dict, aliases, "SurfaceFluxesParameters")
    return SF.Parameters.SurfaceFluxesParameters{FT, AUFP, TP}(; pairs..., ufp, thermo_params)
end

toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
const sf_params = create_parameters(toml_dict, UF.BusingerType())

# atmos sim object and extensions
struct TestAtmos{P, Y, D, I} <: Interfacer.AtmosModelSimulation
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

Interfacer.get_field(sim::TestAtmos, ::Val{:height_int}) = sim.integrator.p.z
Interfacer.get_field(sim::TestAtmos, ::Val{:height_sfc}) = sim.integrator.p.z_sfc
Interfacer.get_field(sim::TestAtmos, ::Val{:uv_int}) = @. StaticArrays.SVector(sim.integrator.p.u, sim.integrator.p.v)
Interfacer.get_field(sim::TestAtmos, ::Val{:thermo_state_int}) =
    TD.PhaseEquil_ρTq.(get_thermo_params(sim), sim.integrator.ρ, sim.integrator.T, sim.integrator.q)
Interfacer.get_field(sim::TestAtmos, ::Val{:air_density}) = sim.integrator.ρ
Interfacer.get_field(sim::TestAtmos, ::Val{:air_temperature}) = sim.integrator.T

function update_sim!(sim::TestAtmos, fields, _)
    (; F_turb_ρτxz, F_turb_energy, F_turb_moisture) = fields
    ρ_int = sim.integrator.ρ
    @. sim.integrator.p.energy_bc = -(F_turb_energy)
    @. sim.integrator.p.ρq_tot_bc = -F_turb_moisture
    @. sim.integrator.p.uₕ_bc = -(F_turb_ρτxz / ρ_int) # x-compoennt only for this test
end

# ocean sim object and extensions
struct TestOcean{M, Y, D, I} <: Interfacer.SurfaceModelSimulation
    model::M
    Y_init::Y
    domain::D
    integrator::I
end
Interfacer.get_field(sim::TestOcean, ::Val{:surface_temperature}) = sim.integrator.T
Interfacer.get_field(sim::TestOcean, ::Val{:air_humidity}) = sim.integrator.p.q
Interfacer.get_field(sim::TestOcean, ::Val{:roughness_momentum}) = sim.integrator.p.z0m
Interfacer.get_field(sim::TestOcean, ::Val{:roughness_buoyancy}) = sim.integrator.p.z0b
Interfacer.get_field(sim::TestOcean, ::Val{:beta}) = sim.integrator.p.beta
Interfacer.get_field(sim::TestOcean, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
Interfacer.get_field(sim::TestOcean, ::Val{:heat_transfer_coefficient}) = sim.integrator.p.Ch
Interfacer.get_field(sim::TestOcean, ::Val{:drag_coefficient}) = sim.integrator.p.Cd
Interfacer.get_field(sim::TestOcean, ::Val{:albedo}) = sim.integrator.p.α

function surface_thermo_state(
    sim::TestOcean,
    thermo_params::TD.Parameters.ThermodynamicsParameters,
    thermo_state_int,
    colidx::Fields.ColumnIndex,
)
    T_sfc = Interfacer.get_field(sim, Val(:surface_temperature), colidx)
    ρ_sfc = thermo_state_int.ρ # arbitrary
    q_sfc = Interfacer.get_field(sim, Val(:air_humidity), colidx) # read from cache
    @. TD.PhaseEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, q_sfc)
end

function update_turbulent_fluxes_point!(sim::TestOcean, fields::NamedTuple, colidx::Fields.ColumnIndex)
    (; F_turb_energy) = fields
    @. sim.integrator.p.F_aero[colidx] = F_turb_energy
end

@testset "calculate correct fluxes: dry" begin
    surface_scheme_list = (MoninObukhovScheme(), BulkScheme())
    for scheme in surface_scheme_list
        boundary_space = TestHelper.create_space(FT)

        params = (; surface_scheme = scheme, FT = FT)

        # atmos
        p = (;
            energy_bc = zeros(boundary_space),
            ρq_tot_bc = zeros(boundary_space),
            uₕ_bc = ones(boundary_space),
            z = ones(boundary_space),
            z_sfc = zeros(boundary_space),
            u = ones(boundary_space),
            v = ones(boundary_space),
        )
        Y_init = (; ρ = ones(boundary_space) .* 1.2, T = ones(boundary_space) .* 310, q = zeros(boundary_space))
        integrator = (; Y_init..., p = p)
        atmos_sim = TestAtmos(params, Y_init, nothing, integrator)

        # ocean
        p = (;
            F_aero = zeros(boundary_space),
            z0m = FT(0.01),
            z0b = FT(0.01),
            beta = ones(boundary_space),
            α = ones(boundary_space) .* 0.5,
            q = zeros(boundary_space),
            Cd = FT(0.01),
            Ch = FT(0.01),
            area_fraction = ones(boundary_space) .* 0.5,
        )
        Y_init = (; T = ones(boundary_space) .* 300.0)
        integrator = (; Y_init..., p = p)
        ocean_sim = TestOcean(nothing, Y_init, nothing, integrator)

        # ocean
        ocean_sim2 = TestOcean(nothing, Y_init, nothing, integrator)

        model_sims = (; atmos_sim, ocean_sim, ocean_sim2)

        coupler_cache_names = (
            :T_S,
            :albedo,
            :F_R_sfc,
            :F_R_toa,
            :P_liq,
            :P_snow,
            :P_net,
            :F_turb_energy,
            :F_turb_ρτxz,
            :F_turb_ρτyz,
            :F_turb_moisture,
        )
        fields = NamedTuple{coupler_cache_names}(ntuple(i -> Fields.zeros(boundary_space), length(coupler_cache_names)))

        # calculate turbulent fluxes
        thermo_params = get_thermo_params(atmos_sim)
        partitioned_turbulent_fluxes!(model_sims, fields, boundary_space, scheme, thermo_params)

        # calculating the fluxes twice ensures that no accumulation occurred (i.e. fluxes are reset to zero each time)
        # TODO: this will need to be extended once flux accumulation is re-enabled
        partitioned_turbulent_fluxes!(model_sims, fields, boundary_space, scheme, thermo_params)

        windspeed = @. hypot(atmos_sim.integrator.p.u, atmos_sim.integrator.p.v)

        thermo_params = get_thermo_params(atmos_sim)
        thermo_state_int = Interfacer.get_field(atmos_sim, Val(:thermo_state_int))

        surface_thermo_states = similar(thermo_state_int)
        Fields.bycolumn(boundary_space) do colidx
            surface_thermo_states[colidx] .=
                surface_thermo_state(ocean_sim, thermo_params, thermo_state_int[colidx], colidx)
        end

        # analytical solution is possible for the BulkScheme() case
        if scheme isa BulkScheme
            ρ_sfc = Interfacer.get_field(atmos_sim, Val(:air_density))
            cpm = TD.cv_m.(thermo_params, thermo_state_int) .+ TD.gas_constant_air.(thermo_params, thermo_state_int) # cp = R + cv
            gz =
                (
                    Interfacer.get_field(atmos_sim, Val(:height_int)) .-
                    Interfacer.get_field(atmos_sim, Val(:height_sfc))
                ) .* 9.81
            shf_analytical = @. (cpm * (ocean_sim.integrator.T - atmos_sim.integrator.T) - gz) *
               ocean_sim.integrator.p.Ch *
               ρ_sfc *
               windspeed #-ρ_sfc * Ch * windspeed(sc) * (cp_m * ΔT + ΔΦ)

            colidx = Fields.ColumnIndex{2}((1, 1), 73) # arbitrary index
            # check the coupler field update
            @test isapprox(parent(shf_analytical[colidx]), parent(fields.F_turb_energy[colidx]), rtol = 1e-6)

            # test the surface field update
            @test parent(fields.F_turb_energy[colidx]) == parent(ocean_sim.integrator.p.F_aero[colidx])

            # test the atmos field update
            update_sim!(atmos_sim, fields, nothing)
            @test parent(fields.F_turb_energy[colidx]) == -parent(atmos_sim.integrator.p.energy_bc[colidx])

        end
        @test parent(fields.F_turb_moisture)[1] ≈ FT(0)

    end
end
@testset "get_surface_params" begin
    @test get_surface_params(TestAtmos([], [], [], [])) == sf_params
    sim = DummySimulation([], [])
    @test_throws ErrorException(
        "get_surface_params is required to be dispatched on" * Interfacer.name(sim) * ", but no method defined",
    ) get_surface_params(DummySimulation([], []))
end

struct DummySurfaceSimulation3{M, Y, D, I} <: Interfacer.SurfaceModelSimulation
    model::M
    Y_init::Y
    domain::D
    integrator::I
end
Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:surface_temperature}) = sim.integrator.T

@testset "update_turbulent_fluxes_point!" begin
    sim = Interfacer.SurfaceStub([])
    sim2 = DummySurfaceSimulation3([], [], [], [])
    colidx = Fields.ColumnIndex{2}((1, 1), 73) # arbitrary index
    @test update_turbulent_fluxes_point!(sim, (;), colidx) == nothing
    @test_throws ErrorException(
        "update_turbulent_fluxes_point! is required to be dispatched on" *
        Interfacer.name(sim2) *
        ", but no method defined",
    ) update_turbulent_fluxes_point!(sim2, (;), colidx) == ErrorException
end

@testset "surface_thermo_state" begin
    boundary_space = TestHelper.create_space(FT)
    _ones = Fields.ones(boundary_space)
    sim = DummySurfaceSimulation3(
        [],
        [],
        [],
        (; T = _ones .* FT(300), ρ = _ones .* FT(1.2), p = (; q = _ones .* FT(0.01))),
    )
    atmos_sim = TestAtmos((; FT = FT), [], [], (; T = _ones .* FT(300), ρ = _ones .* FT(1.2), q = _ones .* FT(0.01)))
    thermo_params = get_thermo_params(atmos_sim)
    thermo_state_int = Interfacer.get_field(atmos_sim, Val(:thermo_state_int))
    colidx = Fields.ColumnIndex{2}((1, 1), 73) # arbitrary index
    @test surface_thermo_state(sim, thermo_params, thermo_state_int[colidx], colidx).ρ == thermo_state_int[colidx].ρ
end
