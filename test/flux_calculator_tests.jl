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
struct TestAtmos2 <: Interfacer.AtmosModelSimulation end

Interfacer.get_field(sim::TestAtmos, ::Val{:air_temperature}) = sim.integrator.T
Interfacer.get_field(sim::TestAtmos, ::Val{:specific_humidity}) = sim.integrator.q
Interfacer.get_field(sim::TestAtmos, ::Val{:air_density}) = sim.integrator.ρ
Interfacer.get_field(sim::TestAtmos, ::Val{:height_int}) = sim.integrator.p.z
Interfacer.get_field(sim::TestAtmos, ::Val{:height_sfc}) = sim.integrator.p.z_sfc
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
struct TestOcean{M, I} <: Interfacer.SurfaceModelSimulation
    model::M
    integrator::I
end

Interfacer.get_field(sim::TestOcean, ::Val{:surface_temperature}) = sim.integrator.T
Interfacer.get_field(sim::TestOcean, ::Val{:roughness_momentum}) = sim.integrator.p.z0m
Interfacer.get_field(sim::TestOcean, ::Val{:roughness_buoyancy}) = sim.integrator.p.z0b
Interfacer.get_field(sim::TestOcean, ::Val{:beta}) = sim.integrator.p.beta
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
struct DummySurfaceSimulation3{M, I} <: Interfacer.SurfaceModelSimulation
    model::M
    integrator::I
end

Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:surface_temperature}) =
    sim.integrator.T
Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:area_fraction}) =
    sim.integrator.p.area_fraction
Interfacer.get_field(sim::DummySurfaceSimulation3, ::Val{:beta}) = sim.integrator.p.beta
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
            z_sfc = zeros(boundary_space),
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
            beta = ones(boundary_space),
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

        coupler_cache_additional = [
            :surface_direct_albedo,
            :surface_diffuse_albedo,
            :P_net,
            :L_MO,
            :ustar,
            :buoyancy_flux,
        ]

        coupler_cache_names = Interfacer.default_coupler_fields()
        push!(coupler_cache_names, coupler_cache_additional...)
        fields = Interfacer.init_coupler_fields(FT, coupler_cache_names, boundary_space)

        # import atmosphere properties into coupler fields
        FieldExchanger.import_atmos_fields!(fields, model_sims)

        # import surface properties into coupler fields
        FieldExchanger.import_combined_surface_fields!(fields, model_sims)

        # calculate turbulent fluxes
        thermo_params = get_thermo_params(atmos_sim)
        FluxCalculator.turbulent_fluxes!(fields, model_sims, thermo_params)

        # calculating the fluxes twice ensures that no accumulation occurred (i.e. fluxes are reset to zero each time)
        # TODO: this will need to be extended once flux accumulation is re-enabled
        FluxCalculator.turbulent_fluxes!(fields, model_sims, thermo_params)

        # Compute expected fluxes
        # Get atmosphere properties
        z_int = Interfacer.get_field(atmos_sim, Val(:height_int))
        uₕ_int =
            StaticArrays.SVector.(
                Interfacer.get_field(atmos_sim, Val(:u_int)),
                Interfacer.get_field(atmos_sim, Val(:v_int)),
            )
        thermo_state_atmos =
            TD.PhaseNonEquil_ρTq.(
                thermo_params,
                atmos_sim.integrator.ρ,
                atmos_sim.integrator.T,
                TD.PhasePartition.(atmos_sim.integrator.q),
            )

        # Get surface properties
        z_sfc = Interfacer.get_field(atmos_sim, Val(:height_sfc))
        z0m = Interfacer.get_field(ocean_sim, Val(:roughness_momentum))
        z0b = Interfacer.get_field(ocean_sim, Val(:roughness_buoyancy))
        gustiness = FT(1)
        beta = Interfacer.get_field(ocean_sim, Val(:beta))

        T_sfc = Interfacer.get_field(ocean_sim, Val(:surface_temperature))
        q_sfc = zeros(boundary_space)
        FluxCalculator.compute_surface_humidity!(
            q_sfc,
            fields.T_atmos,
            fields.q_atmos,
            fields.ρ_atmos,
            T_sfc,
            thermo_params,
        )
        ρ_sfc =
            FluxCalculator.extrapolate_ρ_to_sfc.(thermo_params, thermo_state_atmos, T_sfc)
        thermo_state_sfc =
            TD.PhaseNonEquil_ρTq.(thermo_params, ρ_sfc, T_sfc, TD.PhasePartition.(q_sfc))

        # Use SurfaceFluxes.jl to compute the expected fluxes
        inputs = @. SF.ValuesOnly(
            SF.StateValues(z_int, uₕ_int, thermo_state_atmos), # state_in
            SF.StateValues(                                  # state_sfc
                z_sfc,
                StaticArrays.SVector(FT(0), FT(0)),
                thermo_state_sfc,
            ),
            z0m,
            z0b,
            gustiness,
            beta,
        )
        surface_params = FluxCalculator.get_surface_params(atmos_sim)
        fluxes_expected = FluxCalculator.get_surface_fluxes(inputs, surface_params)

        # Compare expected and computed fluxes
        @test fields.F_turb_ρτxz ≈ fluxes_expected.F_turb_ρτxz
        @test fields.F_turb_ρτyz ≈ fluxes_expected.F_turb_ρτyz
        @test fields.F_lh ≈ fluxes_expected.F_lh
        @test fields.F_sh ≈ fluxes_expected.F_sh
        # The ClimaCore DataLayout underlying the expected moisture flux uses
        # Array instead of SubArray, so we can't compare the fields directly without copy.
        @test copy(fields.F_turb_moisture) ≈ fluxes_expected.F_turb_moisture
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
