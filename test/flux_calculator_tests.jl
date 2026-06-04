import Test: @test, @testset, @test_throws
import StaticArrays
import ClimaCore as CC
import ClimaParams
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters.SurfaceFluxesParameters
import SurfaceFluxes.UniversalFunctions as UF
import ClimaCoupler: FieldExchanger, FluxCalculator, Interfacer

# simple generic atmos model
struct DummySimulation{S, C} <: Interfacer.AbstractAtmosSimulation
    state::S
    cache::C
end

struct DummySimulation2{C} <: Interfacer.AbstractAtmosSimulation
    cache::C
end

# atmos sim object and extensions
struct TestAtmos{P, D, I} <: Interfacer.AbstractAtmosSimulation
    params::P
    domain::D
    integrator::I
end
struct TestAtmos2 <: Interfacer.AbstractAtmosSimulation end

Interfacer.get_field(sim::TestAtmos, ::Val{:air_temperature}) = sim.integrator.T
Interfacer.get_field(sim::TestAtmos, ::Val{:total_specific_humidity}) = sim.integrator.q
Interfacer.get_field(sim::TestAtmos, ::Val{:liquid_specific_humidity}) =
    CC.Fields.zeros(axes(sim.integrator.q))
Interfacer.get_field(sim::TestAtmos, ::Val{:ice_specific_humidity}) =
    CC.Fields.zeros(axes(sim.integrator.q))
Interfacer.get_field(sim::TestAtmos, ::Val{:air_density}) = sim.integrator.ρ
Interfacer.get_field(sim::TestAtmos, ::Val{:height_int}) = sim.integrator.p.z
Interfacer.get_field(sim::TestAtmos, ::Val{:height_sfc}) = sim.integrator.p.height_sfc
Interfacer.get_field(sim::TestAtmos, ::Val{:u_int}) = sim.integrator.p.u
Interfacer.get_field(sim::TestAtmos, ::Val{:v_int}) = sim.integrator.p.v
Interfacer.get_field(sim::TestAtmos, ::Val{:SW_d}) = CC.Fields.zeros(axes(sim.integrator.T))
Interfacer.get_field(sim::TestAtmos, ::Val{:LW_d}) = CC.Fields.zeros(axes(sim.integrator.T))
Interfacer.get_field(sim::TestAtmos, ::Val{:liquid_precipitation}) =
    CC.Fields.zeros(axes(sim.integrator.T))
Interfacer.get_field(sim::TestAtmos, ::Val{:snow_precipitation}) =
    CC.Fields.zeros(axes(sim.integrator.T))
FieldExchanger.update_sim!(sim::TestAtmos, fields) = nothing

function FluxCalculator.update_turbulent_fluxes!(sim::TestAtmos, fields)
    (; F_turb_ρτxz, F_lh, F_sh, F_turb_moisture) = fields
    ρ_int = sim.integrator.ρ
    @. sim.integrator.p.energy_bc = F_lh + F_sh
    @. sim.integrator.p.ρq_tot_bc = F_turb_moisture
    @. sim.integrator.p.uₕ_bc = F_turb_ρτxz / ρ_int # x-component only for this test
    return nothing
end

function get_thermo_params(sim::TestAtmos)
    FT = sim.params.FT
    thermo_params = TDP.ThermodynamicsParameters(FT)
    return thermo_params
end

function FluxCalculator.get_surface_params(sim::TestAtmos)
    FT = sim.params.FT
    sf_params = SurfaceFluxesParameters(FT, UF.BusingerParams)
    return sf_params
end

# ocean sim object and extensions
struct TestOcean{M, I} <: Interfacer.AbstractSurfaceSimulation
    model::M
    integrator::I
end

Interfacer.get_field(sim::TestOcean, ::Val{:surface_temperature}) = sim.integrator.T
Interfacer.get_field(sim::TestOcean, ::Val{:roughness_momentum}) = sim.integrator.p.z0m
Interfacer.get_field(sim::TestOcean, ::Val{:roughness_buoyancy}) = sim.integrator.p.z0b
Interfacer.get_field(sim::TestOcean, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
Interfacer.get_field(
    sim::TestOcean,
    ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}},
) = sim.integrator.p.α
Interfacer.get_field(sim::TestOcean, ::Val{:emissivity}) = eltype(sim.integrator.T)(1)

function FluxCalculator.update_turbulent_fluxes!(sim::TestOcean, fields::NamedTuple)
    (; F_lh, F_sh) = fields
    @. sim.integrator.p.F_aero = F_lh + F_sh
end

# simple surface sim object and extensions
struct DummySurfaceSimulation3{M, I} <: Interfacer.AbstractSurfaceSimulation
    model::M
    integrator::I
end

Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:surface_temperature}) =
    sim.integrator.T
Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:area_fraction}) =
    sim.integrator.p.area_fraction
Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:emissivity}) =
    eltype(sim.integrator.T)(1)

for FT in (Float32, Float64)
    @testset "calculate correct fluxes: dry for FT=$FT" begin
        boundary_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6.371e6), # in meters
            n_quad_points = 4,
            h_elem = 4,
        )

        params = (; FT = FT)

        # atmos
        p = (;
            energy_bc = zeros(boundary_space),
            ρq_tot_bc = zeros(boundary_space),
            uₕ_bc = ones(boundary_space),
            z = ones(boundary_space),
            height_sfc = zeros(boundary_space),
            u = ones(boundary_space),
            v = ones(boundary_space),
        )
        Y_init = (;
            ρ = ones(boundary_space) .* FT(1.2),
            T = ones(boundary_space) .* FT(310),
            q = zeros(boundary_space),
        )
        integrator = (; Y_init..., p = p)
        atmos_sim = TestAtmos(params, nothing, integrator)

        # ocean
        p = (;
            F_aero = zeros(boundary_space),
            z0m = FT(0.01),
            z0b = FT(0.01),
            α = ones(boundary_space) .* FT(0.5),
            q = zeros(boundary_space),
            area_fraction = ones(boundary_space) .* FT(0.5),
        )
        Y_init = (; T = ones(boundary_space) .* FT(300))
        integrator = (; Y_init..., p = p)
        ocean_sim = TestOcean(nothing, integrator)

        # ocean
        ocean_sim2 = TestOcean(nothing, integrator)

        model_sims = (; atmos_sim, ocean_sim, ocean_sim2)

        coupler_cache_additional = [:surface_direct_albedo, :surface_diffuse_albedo, :P_net]

        coupler_cache_names = Interfacer.default_coupler_fields()
        push!(coupler_cache_names, coupler_cache_additional...)
        fields = Interfacer.init_coupler_fields(FT, coupler_cache_names, boundary_space)

        # import static fields into the coupler fields
        FieldExchanger.import_static_fields!(fields, model_sims)

        # import atmosphere properties into coupler fields
        FieldExchanger.import_atmos_fields!(fields, model_sims)

        # import surface properties into coupler fields
        FieldExchanger.import_combined_surface_fields!(fields, model_sims)

        # calculate turbulent fluxes
        thermo_params = get_thermo_params(atmos_sim)
        FluxCalculator.turbulent_fluxes!(fields, model_sims, thermo_params)

        # Calling turbulent_fluxes! twice without any FluxAccumulators must not
        # double-count: the area-weighted coupler `csf.F_*` fields are zeroed
        # at the start of each call, so the final values should equal a single call.
        FluxCalculator.turbulent_fluxes!(fields, model_sims, thermo_params)

        # Compute expected fluxes
        # Get atmosphere properties
        height_int = Interfacer.get_field(atmos_sim, Val(:height_int))
        uv_int =
            StaticArrays.SVector.(
                Interfacer.get_field(atmos_sim, Val(:u_int)),
                Interfacer.get_field(atmos_sim, Val(:v_int)),
            )
        surface_fluxes_params = FluxCalculator.get_surface_params(atmos_sim)

        # Get surface properties
        height_sfc = Interfacer.get_field(atmos_sim, Val(:height_sfc))
        z0m = Interfacer.get_field(ocean_sim, Val(:roughness_momentum))
        z0b = Interfacer.get_field(ocean_sim, Val(:roughness_buoyancy))
        gustiness = FT(1)

        T_sfc = Interfacer.get_field(ocean_sim, Val(:surface_temperature))
        q_sfc = zeros(boundary_space)
        ρ_sfc =
            SF.surface_density.(
                surface_fluxes_params,
                fields.T_atmos,
                fields.ρ_atmos,
                T_sfc,
                height_int .- height_sfc,
                fields.q_tot_atmos,
                0, # q_liq
                0, # q_ice
            )
        q_sfc .= TD.q_vap_saturation.(thermo_params, T_sfc, ρ_sfc, 0, 0)

        config =
            SF.SurfaceFluxConfig.(
                SF.ConstantRoughnessParams.(z0m, z0b),
                SF.ConstantGustinessSpec.(gustiness),
            )

        fluxes_expected =
            FluxCalculator.get_surface_fluxes.(
                surface_fluxes_params,
                uv_int,
                atmos_sim.integrator.T,
                atmos_sim.integrator.q,
                0, # q_liq
                0, # q_ice
                atmos_sim.integrator.ρ,
                height_int,
                StaticArrays.SVector.(0, 0), # uv_sfc
                T_sfc,
                q_sfc,
                height_sfc,
                0, # d
                Ref(config),
            )

        # Compare expected and computed fluxes
        @test fields.F_turb_ρτxz ≈ fluxes_expected.F_turb_ρτxz
        @test fields.F_turb_ρτyz ≈ fluxes_expected.F_turb_ρτyz
        @test fields.F_lh ≈ fluxes_expected.F_lh
        @test fields.F_sh ≈ fluxes_expected.F_sh
        @test fields.F_turb_moisture ≈ fluxes_expected.F_turb_moisture
    end

    @testset "get_surface_params for FT=$FT" begin
        sf_params = SurfaceFluxesParameters(FT, UF.BusingerParams)

        @test FluxCalculator.get_surface_params(TestAtmos((; FT = FT), [], [])) == sf_params
        sim = DummySimulation([], [])
        @test_throws ErrorException(
            "get_surface_params is required to be dispatched on $(nameof(sim)), but no method defined",
        ) FluxCalculator.get_surface_params(DummySimulation([], []))
    end

    @testset "update_turbulent_fluxes! for FT=$FT" begin
        sim = DummySurfaceSimulation3([], [])
        @test_throws ErrorException(
            "update_turbulent_fluxes! is required to be dispatched on $(nameof(sim)), but no method defined",
        ) FluxCalculator.update_turbulent_fluxes!(sim, (;)) == ErrorException
    end
end

@testset "SurfaceStub update_turbulent_fluxes!" begin
    FT = Float32
    @test isnothing(
        FluxCalculator.update_turbulent_fluxes!(Interfacer.SurfaceStub(FT(0)), (;)),
    )
end

# Minimal recording surface sim used only by the FluxAccumulator tests below.
# `update_turbulent_fluxes!` copies each pushed field into `received` so the
# tests can assert that the pushed values match the expected time-average.
mutable struct RecordingSurface{F} <: Interfacer.AbstractSurfaceSimulation
    received::F
    push_count::Int
end
function FluxCalculator.update_turbulent_fluxes!(sim::RecordingSurface, fields)
    sim.received.F_lh .= fields.F_lh
    sim.received.F_sh .= fields.F_sh
    sim.received.F_turb_moisture .= fields.F_turb_moisture
    sim.received.F_turb_ρτxz .= fields.F_turb_ρτxz
    sim.received.F_turb_ρτyz .= fields.F_turb_ρτyz
    sim.push_count += 1
    return nothing
end

@testset "FluxAccumulator behavior" begin
    for FT in (Float32, Float64)
        boundary_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6.371e6),
            n_quad_points = 4,
            h_elem = 4,
        )

        acc = FluxCalculator.FluxAccumulator(boundary_space)
        # Fresh accumulator: zero fields and zero step counter.
        @test acc.n_steps[] == 0
        @test all(parent(acc.F_lh) .== 0)
        @test all(parent(acc.F_sh) .== 0)
        @test all(parent(acc.F_turb_moisture) .== 0)
        @test all(parent(acc.F_turb_ρτxz) .== 0)
        @test all(parent(acc.F_turb_ρτyz) .== 0)

        # Build four sets of per-surface fluxes with known values; accumulate them.
        contributions = [FT(1), FT(3), FT(5), FT(7)]
        for c in contributions
            fields = (;
                F_lh = ones(boundary_space) .* c,
                F_sh = ones(boundary_space) .* (c + FT(10)),
                F_turb_moisture = ones(boundary_space) .* (c + FT(20)),
                F_turb_ρτxz = ones(boundary_space) .* (c + FT(30)),
                F_turb_ρτyz = ones(boundary_space) .* (c + FT(40)),
            )
            FluxCalculator.accumulate!(acc, fields)
        end
        @test acc.n_steps[] == length(contributions)
        expected_sum_lh = sum(contributions)
        @test all(parent(acc.F_lh) .≈ expected_sum_lh)
        @test all(parent(acc.F_sh) .≈ expected_sum_lh + 4 * FT(10))
        @test all(parent(acc.F_turb_moisture) .≈ expected_sum_lh + 4 * FT(20))
        @test all(parent(acc.F_turb_ρτxz) .≈ expected_sum_lh + 4 * FT(30))
        @test all(parent(acc.F_turb_ρτyz) .≈ expected_sum_lh + 4 * FT(40))

        # push_and_reset! divides by n, pushes the averaged values to the
        # surface, and zeros the accumulator.
        received = (;
            F_lh = zeros(boundary_space),
            F_sh = zeros(boundary_space),
            F_turb_moisture = zeros(boundary_space),
            F_turb_ρτxz = zeros(boundary_space),
            F_turb_ρτyz = zeros(boundary_space),
        )
        sim = RecordingSurface(received, 0)
        FluxCalculator.push_and_reset!(sim, acc)
        @test sim.push_count == 1
        n = length(contributions)
        @test all(parent(sim.received.F_lh) .≈ expected_sum_lh / n)
        @test all(parent(sim.received.F_sh) .≈ (expected_sum_lh + 4 * FT(10)) / n)
        @test all(
            parent(sim.received.F_turb_moisture) .≈ (expected_sum_lh + 4 * FT(20)) / n,
        )
        @test all(parent(sim.received.F_turb_ρτxz) .≈ (expected_sum_lh + 4 * FT(30)) / n)
        @test all(parent(sim.received.F_turb_ρτyz) .≈ (expected_sum_lh + 4 * FT(40)) / n)
        # Reset state
        @test acc.n_steps[] == 0
        @test all(parent(acc.F_lh) .== 0)

        # push_and_reset! with n_steps == 0 is a no-op (doesn't divide by zero,
        # doesn't push).
        FluxCalculator.push_and_reset!(sim, acc)
        @test sim.push_count == 1  # unchanged

        # A second accumulate-then-push cycle behaves the same as the first.
        for c in (FT(2), FT(4))
            fields = (;
                F_lh = ones(boundary_space) .* c,
                F_sh = ones(boundary_space) .* (c + FT(10)),
                F_turb_moisture = ones(boundary_space) .* (c + FT(20)),
                F_turb_ρτxz = ones(boundary_space) .* (c + FT(30)),
                F_turb_ρτyz = ones(boundary_space) .* (c + FT(40)),
            )
            FluxCalculator.accumulate!(acc, fields)
        end
        FluxCalculator.push_and_reset!(sim, acc)
        @test sim.push_count == 2
        @test all(parent(sim.received.F_lh) .≈ FT(3))  # mean of 2 and 4
        @test acc.n_steps[] == 0
    end
end

# Call `turbulent_fluxes!` with a FluxAccumulator, verify the per-surface flux is
# accumulated each call, then `push_and_reset!` averages the contributions and pushes
# to the surface model.
@testset "turbulent_fluxes! with FluxAccumulator" begin
    for FT in (Float32, Float64)
        boundary_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6.371e6),
            n_quad_points = 4,
            h_elem = 4,
        )

        params = (; FT = FT)
        p_atmos = (;
            energy_bc = zeros(boundary_space),
            ρq_tot_bc = zeros(boundary_space),
            uₕ_bc = ones(boundary_space),
            z = ones(boundary_space),
            height_sfc = zeros(boundary_space),
            u = ones(boundary_space),
            v = ones(boundary_space),
        )
        Y_atmos = (;
            ρ = ones(boundary_space) .* FT(1.2),
            T = ones(boundary_space) .* FT(310),
            q = zeros(boundary_space),
        )
        atmos_sim = TestAtmos(params, nothing, (; Y_atmos..., p = p_atmos))

        # Slow ocean: gets accumulator. We will run `turbulent_fluxes!` three
        # times with a slowly increasing T_sfc, then verify push_and_reset!
        # averages the three flux contributions correctly.
        p_ocean = (;
            F_aero = zeros(boundary_space),
            z0m = FT(0.01),
            z0b = FT(0.01),
            α = ones(boundary_space) .* FT(0.5),
            q = zeros(boundary_space),
            area_fraction = ones(boundary_space) .* FT(1),
        )
        ocean_sim = TestOcean(nothing, (; T = ones(boundary_space) .* FT(300), p = p_ocean))

        model_sims = (; atmos_sim, ocean_sim)

        # Build coupler fields and the accumulator for the slow ocean.
        coupler_cache_names = Interfacer.default_coupler_fields()
        push!(coupler_cache_names, :surface_direct_albedo, :surface_diffuse_albedo, :P_net)
        fields = Interfacer.init_coupler_fields(FT, coupler_cache_names, boundary_space)
        FieldExchanger.import_static_fields!(fields, model_sims)
        FieldExchanger.import_atmos_fields!(fields, model_sims)
        FieldExchanger.import_combined_surface_fields!(fields, model_sims)

        ocean_acc = FluxCalculator.FluxAccumulator(boundary_space)
        flux_accumulators = (; ocean_sim = ocean_acc)
        thermo_params = get_thermo_params(atmos_sim)

        # Snapshot F_aero before any flux computation: it should be zero
        # initially. The slow-surface branch should leave it untouched until
        # `push_and_reset!` runs.
        @test all(parent(ocean_sim.integrator.p.F_aero) .== 0)

        # Three coupling-step calls of turbulent_fluxes!. Vary ocean T_sfc so
        # each call produces a different per-surface flux, exercising the
        # averaging path.
        for (i, Tsfc) in enumerate((FT(300), FT(302), FT(304)))
            ocean_sim.integrator.T .= Tsfc
            FluxCalculator.turbulent_fluxes!(
                fields,
                model_sims,
                thermo_params,
                flux_accumulators,
            )
            @test ocean_acc.n_steps[] == i
            # F_aero must still be zero because the slow-surface push has not
            # happened yet — fluxes went into the accumulator only.
            @test all(parent(ocean_sim.integrator.p.F_aero) .== 0)
        end

        # Capture the running-sum values before the divide so we can compute
        # the expected mean and compare to what gets pushed.
        sum_lh = copy(parent(ocean_acc.F_lh))
        sum_sh = copy(parent(ocean_acc.F_sh))
        n = ocean_acc.n_steps[]
        @test n == 3

        FluxCalculator.push_and_reset!(ocean_sim, ocean_acc)

        # After push_and_reset!, the slow surface received the time-averaged
        # F_lh + F_sh in its F_aero BC (per TestOcean's update_turbulent_fluxes!).
        expected_F_aero = sum_lh ./ n .+ sum_sh ./ n
        @test parent(ocean_sim.integrator.p.F_aero) ≈ expected_F_aero
        @test ocean_acc.n_steps[] == 0
        @test all(parent(ocean_acc.F_lh) .== 0)
    end
end
