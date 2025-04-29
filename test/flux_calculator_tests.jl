import Test: @test, @testset, @test_throws
import StaticArrays
import ClimaCore as CC
import ClimaParams
import Thermodynamics as TD
import Thermodynamics.Parameters.ThermodynamicsParameters
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import SurfaceFluxes.UniversalFunctions as UF
import ClimaCoupler: FieldExchanger, FluxCalculator, Interfacer

# simple generic atmos model
struct DummySimulation{S, C} <: Interfacer.AtmosModelSimulation
    state::S
    cache::C
end

struct DummySimulation2{C} <: Interfacer.AtmosModelSimulation
    cache::C
end

# atmos sim object and extensions
struct TestAtmos{P, D, I} <: Interfacer.AtmosModelSimulation
    params::P
    domain::D
    integrator::I
end
Interfacer.name(sim::TestAtmos) = "TestAtmos"
struct TestAtmos2 <: Interfacer.AtmosModelSimulation end
Interfacer.name(sim::TestAtmos2) = "TestAtmos2"

Interfacer.get_field(sim::TestAtmos, ::Val{:height_int}) = sim.integrator.p.z
Interfacer.get_field(sim::TestAtmos, ::Val{:height_sfc}) = sim.integrator.p.z_sfc
Interfacer.get_field(sim::TestAtmos, ::Val{:u_int}) = sim.integrator.p.u
Interfacer.get_field(sim::TestAtmos, ::Val{:v_int}) = sim.integrator.p.v
Interfacer.get_field(sim::TestAtmos, ::Val{:thermo_state_int}) =
    TD.PhaseEquil_ρTq.(get_thermo_params(sim), sim.integrator.ρ, sim.integrator.T, sim.integrator.q)

function FieldExchanger.update_sim!(sim::TestAtmos, fields, _)
    (; F_turb_ρτxz, F_lh, F_sh, F_turb_moisture) = fields
    ρ_int = sim.integrator.ρ
    @. sim.integrator.p.energy_bc = F_lh + F_sh
    @. sim.integrator.p.ρq_tot_bc = F_turb_moisture
    @. sim.integrator.p.uₕ_bc = F_turb_ρτxz / ρ_int # x-component only for this test
end

function get_thermo_params(sim::TestAtmos)
    FT = sim.params.FT
    thermo_params = ThermodynamicsParameters(FT)
    return thermo_params
end

function FluxCalculator.get_surface_params(sim::TestAtmos)
    FT = sim.params.FT
    sf_params = SurfaceFluxesParameters(FT, UF.BusingerParams)
    return sf_params
end


# ocean sim object and extensions
struct TestOcean{M, I} <: Interfacer.SurfaceModelSimulation
    model::M
    integrator::I
end
Interfacer.name(sim::TestOcean) = "TestOcean"

Interfacer.get_field(sim::TestOcean, ::Val{:surface_temperature}) = sim.integrator.T
Interfacer.get_field(sim::TestOcean, ::Val{:air_humidity}) = sim.integrator.p.q
Interfacer.get_field(sim::TestOcean, ::Val{:roughness_momentum}) = sim.integrator.p.z0m
Interfacer.get_field(sim::TestOcean, ::Val{:roughness_buoyancy}) = sim.integrator.p.z0b
Interfacer.get_field(sim::TestOcean, ::Val{:beta}) = sim.integrator.p.beta
Interfacer.get_field(sim::TestOcean, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
Interfacer.get_field(sim::TestOcean, ::Val{:heat_transfer_coefficient}) = sim.integrator.p.Ch
Interfacer.get_field(sim::TestOcean, ::Val{:drag_coefficient}) = sim.integrator.p.Cd
Interfacer.get_field(sim::TestOcean, ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}}) =
    sim.integrator.p.α

function FluxCalculator.surface_thermo_state(sim::TestOcean, thermo_params::ThermodynamicsParameters, thermo_state_int)
    T_sfc = Interfacer.get_field(sim, Val(:surface_temperature))
    ρ_sfc = thermo_state_int.ρ # arbitrary
    q_sfc = Interfacer.get_field(sim, Val(:air_humidity)) # read from cache
    @. TD.PhaseEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, q_sfc)
end

function FluxCalculator.update_turbulent_fluxes!(sim::TestOcean, fields::NamedTuple)
    (; F_lh, F_sh) = fields
    @. sim.integrator.p.F_aero = F_lh + F_sh
end

# simple surface sim object and extensions
struct DummySurfaceSimulation3{M, I} <: Interfacer.SurfaceModelSimulation
    model::M
    integrator::I
end
Interfacer.name(sim::DummySurfaceSimulation3) = "DummySurfaceSimulation3"

Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:surface_temperature}) = sim.integrator.T
Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:heat_transfer_coefficient}) = sim.integrator.p.Ch
Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:drag_coefficient}) = sim.integrator.p.Cd
Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:beta}) = sim.integrator.p.beta

function FluxCalculator.surface_thermo_state(
    sim::DummySurfaceSimulation3,
    thermo_params::ThermodynamicsParameters,
    thermo_state_int,
)
    T_sfc = Interfacer.get_field(sim, Val(:surface_temperature))
    FT = eltype(T_sfc)

    ρ_sfc = @. T_sfc * FT(0) .+ FT(1.2) # arbitrary
    q_sfc = @. T_sfc * FT(0) # dry surface
    @. TD.PhaseEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, q_sfc)
end

function FluxCalculator.water_albedo_from_atmosphere!(::TestAtmos, temp1::CC.Fields.Field, temp2::CC.Fields.Field)
    temp1 .*= 2
    temp2 .*= 3
end

for FT in (Float32, Float64)
    @testset "calculate_surface_air_density for FT=$FT" begin
        boundary_space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4)
        coupler_fields = (; T_sfc = 310 .* ones(boundary_space))
        sim2 = DummySimulation2((; cache = (; flux = zeros(boundary_space))))
        @test_throws ErrorException FluxCalculator.calculate_surface_air_density(sim2, coupler_fields.T_sfc)
    end

    @testset "calculate correct fluxes: dry for FT=$FT" begin
        boundary_space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4)

        params = (; FT = FT)

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
        Y_init = (; ρ = ones(boundary_space) .* FT(1.2), T = ones(boundary_space) .* FT(310), q = zeros(boundary_space))
        integrator = (; Y_init..., p = p)
        atmos_sim = TestAtmos(params, nothing, integrator)

        # ocean
        p = (;
            F_aero = zeros(boundary_space),
            z0m = FT(0.01),
            z0b = FT(0.01),
            beta = ones(boundary_space),
            α = ones(boundary_space) .* FT(0.5),
            q = zeros(boundary_space),
            Cd = FT(0.01),
            Ch = FT(0.01),
            area_fraction = ones(boundary_space) .* FT(0.5),
        )
        Y_init = (; T = ones(boundary_space) .* FT(300))
        integrator = (; Y_init..., p = p)
        ocean_sim = TestOcean(nothing, integrator)

        # ocean
        ocean_sim2 = TestOcean(nothing, integrator)

        model_sims = (; atmos_sim, ocean_sim, ocean_sim2)

        coupler_cache_names = [
            :T_sfc,
            :surface_direct_albedo,
            :surface_diffuse_albedo,
            :P_liq,
            :P_snow,
            :P_net,
            :F_lh,
            :F_sh,
            :F_turb_ρτxz,
            :F_turb_ρτyz,
            :F_turb_moisture,
            :z0m_sfc,
            :z0b_sfc,
            :beta,
            :q_sfc,
            :L_MO,
            :ustar,
            :buoyancy_flux,
        ]
        fields = Interfacer.init_coupler_fields(FT, coupler_cache_names, boundary_space)

        # calculate turbulent fluxes
        thermo_params = get_thermo_params(atmos_sim)
        FluxCalculator.turbulent_fluxes!(model_sims, fields, boundary_space, thermo_params)

        # calculating the fluxes twice ensures that no accumulation occurred (i.e. fluxes are reset to zero each time)
        # TODO: this will need to be extended once flux accumulation is re-enabled
        FluxCalculator.turbulent_fluxes!(model_sims, fields, boundary_space, thermo_params)

        windspeed = @. hypot(atmos_sim.integrator.p.u, atmos_sim.integrator.p.v)

        thermo_params = get_thermo_params(atmos_sim)
        thermo_state_int = Interfacer.get_field(atmos_sim, Val(:thermo_state_int))

        surface_thermo_states = similar(thermo_state_int)
        surface_thermo_states .= FluxCalculator.surface_thermo_state(ocean_sim, thermo_params, thermo_state_int)

        # NOTE: This test is very weak! We should add more stringent tests
        @test Array(parent(fields.F_turb_moisture))[1] ≈ FT(0)
        for p in coupler_cache_names
            @test !any(isnan, getproperty(fields, p))
        end
    end

    @testset "get_surface_params for FT=$FT" begin
        sf_params = SurfaceFluxesParameters(FT, UF.BusingerParams)

        @test FluxCalculator.get_surface_params(TestAtmos((; FT = FT), [], [])) == sf_params
        sim = DummySimulation([], [])
        @test_throws ErrorException(
            "get_surface_params is required to be dispatched on" * Interfacer.name(sim) * ", but no method defined",
        ) FluxCalculator.get_surface_params(DummySimulation([], []))
    end

    @testset "update_turbulent_fluxes! for FT=$FT" begin
        sim = DummySurfaceSimulation3([], [])
        @test_throws ErrorException(
            "update_turbulent_fluxes! is required to be dispatched on" *
            Interfacer.name(sim) *
            ", but no method defined",
        ) FluxCalculator.update_turbulent_fluxes!(sim, (;)) == ErrorException
    end

    @testset "surface_thermo_state for FT=$FT" begin
        boundary_space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4)
        _ones = CC.Fields.ones(boundary_space)
        surface_sim =
            DummySurfaceSimulation3([], (; T = _ones .* FT(300), ρ = _ones .* FT(1.2), p = (; q = _ones .* FT(0.01))))
        atmos_sim = TestAtmos((; FT = FT), [], (; T = _ones .* FT(300), ρ = _ones .* FT(1.2), q = _ones .* FT(0.01)))
        thermo_params = get_thermo_params(atmos_sim)
        thermo_state_int = Interfacer.get_field(atmos_sim, Val(:thermo_state_int))
        @test FluxCalculator.surface_thermo_state(surface_sim, thermo_params, thermo_state_int).ρ == thermo_state_int.ρ
    end

    @testset "water_albedo_from_atmosphere!" begin
        boundary_space = CC.CommonSpaces.CubedSphereSpace(FT; radius = FT(6371e3), n_quad_points = 4, h_elem = 4)
        ocean_sim = Interfacer.SurfaceStub((; α_direct = zeros(boundary_space), α_diffuse = zeros(boundary_space)))
        atmos_sim = TestAtmos(1, 2, 3)
        coupler_fields = (; temp1 = ones(boundary_space), temp2 = ones(boundary_space))
        model_sims = (; atmos_sim, ocean_sim)
        cs = Interfacer.CoupledSimulation{FT}(
            nothing, # comms_ctx
            nothing, # dates
            nothing, # boundary_space
            coupler_fields, # fields
            nothing, # conservation_checks
            (Int(0), Int(1)), # tspan
            0, # Δt_cpl
            Ref(Int(0)), # t
            model_sims, # model_sims
            (;), # callbacks
            (;), # dirs
            nothing, # thermo_params
            nothing, # diags_handler
        )
        FluxCalculator.water_albedo_from_atmosphere!(cs)
        @test sum(parent(cs.model_sims.ocean_sim.cache.α_direct) .- parent(ones(boundary_space)) .* 2) == 0
        @test sum(parent(cs.model_sims.ocean_sim.cache.α_diffuse) .- parent(ones(boundary_space)) .* 3) == 0

        atmos_sim2 = TestAtmos2()
        @test_throws ErrorException(
            "this function is required to be dispatched on" * Interfacer.name(atmos_sim2) * ", but no method defined",
        ) FluxCalculator.water_albedo_from_atmosphere!(atmos_sim2, ones(boundary_space), ones(boundary_space))
    end
end

@testset "SurfaceStub update_turbulent_fluxes!" begin
    FT = Float32
    @test isnothing(FluxCalculator.update_turbulent_fluxes!(Interfacer.SurfaceStub(FT(0)), (;)))
end
