import Test: @test, @testset
import ClimaCore as CC
import ClimaCoupler: Interfacer, FieldExchanger, FluxCalculator

include(joinpath("..", "experiments", "ClimaEarth", "test", "TestHelper.jl"))
import .TestHelper

# test for a simple generic atmos model
struct DummySimulation{C} <: Interfacer.AtmosModelSimulation
    cache::C
end
Interfacer.name(::DummySimulation) = "DummySimulation"
Interfacer.get_field(sim::DummySimulation, ::Val{:turbulent_energy_flux}) = sim.cache.turbulent_energy_flux
Interfacer.get_field(sim::DummySimulation, ::Val{:turbulent_moisture_flux}) = sim.cache.turbulent_moisture_flux
Interfacer.get_field(sim::DummySimulation, ::Val{:radiative_energy_flux_sfc}) = sim.cache.radiative_energy_flux_sfc
Interfacer.get_field(sim::DummySimulation, ::Val{:liquid_precipitation}) = sim.cache.liquid_precipitation
Interfacer.get_field(sim::DummySimulation, ::Val{:snow_precipitation}) = sim.cache.snow_precipitation

function FluxCalculator.calculate_surface_air_density(atmos_sim::DummySimulation, T_sfc::CC.Fields.Field)
    FT = CC.Spaces.undertype(axes(T_sfc))
    return T_sfc .* FT(0.0) .+ FT(1.0)
end


# surface field exchange tests
struct TestSurfaceSimulation1{C} <: Interfacer.SurfaceModelSimulation
    cache_field::C
end
Interfacer.name(::TestSurfaceSimulation1) = "TestSurfaceSimulation1"
struct TestSurfaceSimulation2{C} <: Interfacer.SurfaceModelSimulation
    cache_field::C
end
Interfacer.name(::TestSurfaceSimulation2) = "TestSurfaceSimulation2"

Interfacer.get_field(sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2}, ::Val{:surface_temperature}) =
    sim.cache_field
Interfacer.get_field(
    sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2},
    ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}},
) = sim.cache_field
Interfacer.get_field(sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2}, ::Val{:roughness_momentum}) =
    sim.cache_field
Interfacer.get_field(sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2}, ::Val{:roughness_buoyancy}) =
    sim.cache_field
Interfacer.get_field(sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2}, ::Val{:beta}) = sim.cache_field

Interfacer.get_field(sim::TestSurfaceSimulation1, ::Val{:area_fraction}) = sim.cache_field .* 0
Interfacer.get_field(sim::TestSurfaceSimulation2, ::Val{:area_fraction}) =
    sim.cache_field .* CC.Spaces.undertype(axes(sim.cache_field))(0.5)

Interfacer.get_field(sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2}, ::Val{:surface_humidity}) =
    sim.cache_field .* 0
Interfacer.get_field(sim::Union{TestSurfaceSimulation2, TestSurfaceSimulation2}, ::Val{:surface_humidity}) =
    sim.cache_field .* 0

Interfacer.reinit!(::TestSurfaceSimulation1) = nothing
Interfacer.step!(::TestSurfaceSimulation1, _) = nothing

struct TestSurfaceSimulationA <: Interfacer.SurfaceModelSimulation end
struct TestSurfaceSimulationB <: Interfacer.SurfaceModelSimulation end
struct TestSurfaceSimulationC <: Interfacer.SurfaceModelSimulation end
struct TestSurfaceSimulationD <: Interfacer.SurfaceModelSimulation end

# Initialize weights (fractions) and initial values (fields)
Interfacer.get_field(::TestSurfaceSimulationA, ::Val{:random}) = 1.0
Interfacer.get_field(::TestSurfaceSimulationB, ::Val{:random}) = 0.5
Interfacer.get_field(::TestSurfaceSimulationC, ::Val{:random}) = -1.0
Interfacer.get_field(::TestSurfaceSimulationD, ::Val{:random}) = 0.75

Interfacer.get_field(::TestSurfaceSimulationA, ::Val{:area_fraction}) = 0.0
Interfacer.get_field(::TestSurfaceSimulationB, ::Val{:area_fraction}) = 0.5
Interfacer.get_field(::TestSurfaceSimulationC, ::Val{:area_fraction}) = 2.0
Interfacer.get_field(::TestSurfaceSimulationD, ::Val{:area_fraction}) = -10.0

struct DummyStub{C} <: Interfacer.SurfaceModelSimulation
    cache::C
end
Interfacer.get_field(sim::DummyStub, ::Val{:area_fraction}) = sim.cache.area_fraction
function Interfacer.update_field!(sim::DummyStub, ::Val{:area_fraction}, field::CC.Fields.Field)
    sim.cache.area_fraction .= field
end
# atmos sim
struct TestAtmosSimulation{C} <: Interfacer.AtmosModelSimulation
    cache::C
end
Interfacer.name(::TestAtmosSimulation) = "TestAtmosSimulation"
function Interfacer.update_field!(sim::TestAtmosSimulation, ::Val{:surface_direct_albedo}, field)
    parent(sim.cache.albedo_direct) .= parent(field)
end
function Interfacer.update_field!(sim::TestAtmosSimulation, ::Val{:surface_diffuse_albedo}, field)
    parent(sim.cache.albedo_diffuse) .= parent(field)
end
function Interfacer.update_field!(sim::TestAtmosSimulation, ::Val{:roughness_momentum}, field)
    parent(sim.cache.roughness_momentum) .= parent(field)
end
Interfacer.update_field!(sim::TestAtmosSimulation, ::Val{:surface_temperature}, field) = nothing
Interfacer.update_field!(sim::TestAtmosSimulation, ::Val{:roughness_buoyancy}, field) = nothing
Interfacer.update_field!(sim::TestAtmosSimulation, ::Val{:beta}, field) = nothing
Interfacer.update_field!(sim::TestAtmosSimulation, ::Val{:turbulent_fluxes}, field) = nothing

#surface sim
struct TestSurfaceSimulationLand{C} <: Interfacer.SurfaceModelSimulation
    cache::C
end
Interfacer.name(::TestSurfaceSimulationLand) = "TestSurfaceSimulationLand"
function Interfacer.get_field(sim::TestSurfaceSimulationLand, ::Val{:area_fraction})
    FT = CC.Spaces.undertype(axes(sim.cache.turbulent_energy_flux))
    return FT(0.5)
end
function Interfacer.update_field!(sim::TestSurfaceSimulationLand, ::Val{:turbulent_energy_flux}, field)
    parent(sim.cache.turbulent_energy_flux) .= parent(field)
end
function Interfacer.update_field!(sim::TestSurfaceSimulationLand, ::Val{:turbulent_moisture_flux}, field)
    parent(sim.cache.turbulent_moisture_flux) .= parent(field)
end
function Interfacer.update_field!(sim::TestSurfaceSimulationLand, ::Val{:air_density}, field)
    parent(sim.cache.air_density) .= parent(field)
end
Interfacer.update_field!(sim::TestSurfaceSimulationLand, ::Val{:radiative_energy_flux_sfc}, field) = nothing
Interfacer.update_field!(sim::TestSurfaceSimulationLand, ::Val{:liquid_precipitation}, field) = nothing
Interfacer.update_field!(sim::TestSurfaceSimulationLand, ::Val{:snow_precipitation}, field) = nothing

for FT in (Float32, Float64)
    @testset "test dummmy_remap!" begin
        test_space = TestHelper.create_space(FT)
        test_field_ones = CC.Fields.ones(test_space)
        target_field = CC.Fields.zeros(test_space)

        FieldExchanger.dummmy_remap!(target_field, test_field_ones)
        @test parent(target_field) == parent(test_field_ones)
    end

    @testset "test update_surface_fractions!" begin
        test_space = TestHelper.create_space(FT)
        # Construct land fraction of 0s in top half, 1s in bottom half
        land_fraction = CC.Fields.ones(test_space)
        dims = size(parent(land_fraction))
        m = dims[1]
        n = dims[2]
        parent(land_fraction)[1:(m ÷ 2), :, :, :] .= FT(0)

        # Construct ice fraction of 0s on left, 0.5s on right
        ice_d = CC.Fields.zeros(test_space)
        parent(ice_d)[:, (n ÷ 2 + 1):n, :, :] .= FT(0.5)

        # Construct ice fraction of 0s on left, 0.5s on right
        ocean_d = CC.Fields.zeros(test_space)

        # Fill in only the necessary parts of the simulation
        cs = Interfacer.CoupledSimulation{FT}(
            nothing, # comms_ctx
            nothing, # dates
            nothing, # boundary_space
            nothing, # fields
            nothing, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # Δt_cpl
            (;
                ice_sim = DummyStub((; area_fraction = ice_d)),
                ocean_sim = Interfacer.SurfaceStub((; area_fraction = ocean_d)),
                land_sim = DummyStub((; area_fraction = land_fraction)),
            ), # model_sims
            (;), # callbacks
            (;), # dirs
            nothing, # turbulent_fluxes
            nothing, # thermo_params
            nothing, # diags_handler
        )

        FieldExchanger.update_surface_fractions!(cs)

        # Test that sum of fractions is 1 everywhere
        ice_fraction = Interfacer.get_field(cs.model_sims.ice_sim, Val(:area_fraction))
        ocean_fraction = Interfacer.get_field(cs.model_sims.ocean_sim, Val(:area_fraction))
        @test all(parent(ice_fraction .+ ocean_fraction .+ land_fraction) .== FT(1))
    end

    @testset "test combine_surfaces" begin
        test_space = TestHelper.create_space(FT)
        combined_field = CC.Fields.ones(test_space)

        var_name = Val(:random)
        sims = (;
            a = TestSurfaceSimulationA(),
            b = TestSurfaceSimulationB(),
            c = TestSurfaceSimulationC(),
            d = TestSurfaceSimulationD(),
        )

        fractions = (
            Interfacer.get_field(sims.a, Val(:area_fraction)),
            Interfacer.get_field(sims.b, Val(:area_fraction)),
            Interfacer.get_field(sims.c, Val(:area_fraction)),
            Interfacer.get_field(sims.d, Val(:area_fraction)),
        )
        fields = (
            Interfacer.get_field(sims.a, var_name),
            Interfacer.get_field(sims.b, var_name),
            Interfacer.get_field(sims.c, var_name),
            Interfacer.get_field(sims.d, var_name),
        )

        FieldExchanger.combine_surfaces!(combined_field, sims, var_name)
        @test all(parent(combined_field) .== FT(sum(fractions .* fields)))
    end

    @testset "import_atmos_fields! for FT=$FT" begin
        boundary_space = TestHelper.create_space(FT)
        coupler_names = (:F_turb_energy, :F_turb_moisture, :F_radiative, :P_liq, :P_snow, :ρ_sfc, :T_sfc)
        component_names = (
            :turbulent_energy_flux,
            :turbulent_moisture_flux,
            :radiative_energy_flux_sfc,
            :liquid_precipitation,
            :snow_precipitation,
        )
        component_fields =
            NamedTuple{component_names}(ntuple(i -> CC.Fields.ones(boundary_space), length(component_names)))

        model_sims =
            (; atmos_sim = DummySimulation(component_fields), land_sim = TestSurfaceSimulation1(component_fields))

        flux_types = (FluxCalculator.CombinedStateFluxesMOST(), FluxCalculator.PartitionedStateFluxes())
        results = [FT(1), FT(0)]
        for (i, t) in enumerate(flux_types)
            coupler_fields =
                NamedTuple{coupler_names}(ntuple(i -> CC.Fields.zeros(boundary_space), length(coupler_names)))
            FieldExchanger.import_atmos_fields!(coupler_fields, model_sims, boundary_space, t)
            @test Array(parent(coupler_fields.F_turb_energy))[1] == results[i]
            @test Array(parent(coupler_fields.F_turb_moisture))[1] == results[i]
            @test Array(parent(coupler_fields.F_radiative))[1] == results[1]
            @test Array(parent(coupler_fields.P_liq))[1] == results[1]
            @test Array(parent(coupler_fields.P_snow))[1] == results[1]
        end
    end

    @testset "import_combined_surface_fields! for FT=$FT" begin
        # coupler cache setup
        boundary_space = TestHelper.create_space(FT)
        coupler_names =
            (:T_sfc, :z0m_sfc, :z0b_sfc, :surface_direct_albedo, :surface_diffuse_albedo, :beta, :q_sfc, :temp1)

        # coupler cache setup
        exchanged_fields = (
            :surface_temperature,
            :surface_direct_albedo,
            :surface_diffuse_albedo,
            :roughness_momentum,
            :roughness_buoyancy,
            :beta,
        )

        sims = (; a = TestSurfaceSimulation1(ones(boundary_space)), b = TestSurfaceSimulation2(ones(boundary_space)))

        # test the coupler update under CombinedStateFluxesMOST (update all) and PartitionedStateFluxes (update all except
        # surface info needed for the calculation of turbulent fluxes)
        flux_types = (FluxCalculator.CombinedStateFluxesMOST(), FluxCalculator.PartitionedStateFluxes())
        results = [FT(0.5), FT(0)] # 0.5 due to the area fraction weighting
        for (i, t) in enumerate(flux_types)
            coupler_fields =
                NamedTuple{coupler_names}(ntuple(i -> CC.Fields.zeros(boundary_space), length(coupler_names)))
            FieldExchanger.import_combined_surface_fields!(coupler_fields, sims, t)
            @test Array(parent(coupler_fields.T_sfc))[1] == results[1]
            @test Array(parent(coupler_fields.surface_direct_albedo))[1] == results[1]
            @test Array(parent(coupler_fields.surface_diffuse_albedo))[1] == results[1]
            @test Array(parent(coupler_fields.z0m_sfc))[1] == results[i]
            @test Array(parent(coupler_fields.z0b_sfc))[1] == results[i]
            @test Array(parent(coupler_fields.beta))[1] == results[i]
        end
    end

    @testset "update_model_sims! for FT=$FT" begin
        # coupler cache setup
        boundary_space = TestHelper.create_space(FT)
        coupler_field_names = (
            :ρ_sfc,
            :T_sfc,
            :z0m_sfc,
            :z0b_sfc,
            :surface_direct_albedo,
            :surface_diffuse_albedo,
            :beta,
            :F_turb_energy,
            :F_turb_moisture,
            :F_radiative,
            :P_liq,
            :P_snow,
        )
        coupler_fields =
            NamedTuple{coupler_field_names}(ntuple(i -> CC.Fields.ones(boundary_space), length(coupler_field_names)))

        # model cache setup

        atmos_names =
            (:surface_temperature, :albedo_direct, :albedo_diffuse, :roughness_momentum, :roughness_buoyancy, :beta)
        atmos_fields = NamedTuple{atmos_names}(ntuple(i -> CC.Fields.zeros(boundary_space), length(atmos_names)))

        land_names = (
            :turbulent_energy_flux,
            :turbulent_moisture_flux,
            :air_density,
            :radiative_energy_flux_sfc,
            :liquid_precipitation,
            :snow_precipitation,
        )
        land_fields = NamedTuple{land_names}(ntuple(i -> CC.Fields.zeros(boundary_space), length(land_names)))

        model_sims = (;
            atmos_sim = TestAtmosSimulation(atmos_fields),
            land_sim = TestSurfaceSimulationLand(land_fields),
            stub_sim = Interfacer.SurfaceStub((;
                area_fraction = CC.Fields.ones(boundary_space),
                ρ_sfc = CC.Fields.ones(boundary_space),
                albedo_direct = CC.Fields.ones(boundary_space),
                albedo_diffuse = CC.Fields.ones(boundary_space),
            )),
        )
        coupler_fields.surface_diffuse_albedo .= FT(0.5)

        # test the sim update under CombinedStateFluxesMOST (update all) and PartitionedStateFluxes (update all except turbulent fluxes)
        flux_types = (FluxCalculator.CombinedStateFluxesMOST(), FluxCalculator.PartitionedStateFluxes())
        results = [FT(0), FT(1), FT(0.5)]
        for (i, t) in enumerate(flux_types)
            model_sims.atmos_sim.cache.roughness_momentum .= FT(0)
            FieldExchanger.update_model_sims!(model_sims, coupler_fields, t)

            # test atmos
            @test Array(parent(model_sims.atmos_sim.cache.albedo_direct))[1] == results[2]
            @test Array(parent(model_sims.atmos_sim.cache.albedo_diffuse))[1] == results[3]
            if t isa FluxCalculator.CombinedStateFluxesMOST
                @test Array(parent(model_sims.atmos_sim.cache.roughness_momentum))[1] == results[2]
            else
                @test Array(parent(model_sims.atmos_sim.cache.roughness_momentum))[1] == results[1]
            end

            # test variables without updates
            @test Array(parent(model_sims.atmos_sim.cache.surface_temperature))[1] == results[1]
            @test Array(parent(model_sims.atmos_sim.cache.beta))[1] == results[1]
            @test Array(parent(model_sims.atmos_sim.cache.roughness_buoyancy))[1] == results[1]

            # test surface
            @test Array(parent(model_sims.land_sim.cache.turbulent_energy_flux))[1] == results[2] # assuming units / m2
            @test Array(parent(model_sims.land_sim.cache.turbulent_moisture_flux))[1] == results[2]

            # test variables without updates
            @test Array(parent(model_sims.land_sim.cache.radiative_energy_flux_sfc))[1] == results[1]
            @test Array(parent(model_sims.land_sim.cache.liquid_precipitation))[1] == results[1]
            @test Array(parent(model_sims.land_sim.cache.snow_precipitation))[1] == results[1]

            # test stub - albedo should be updated by update_sim!
            @test Array(parent(model_sims.stub_sim.cache.albedo_direct))[1] == results[2]
            @test Array(parent(model_sims.stub_sim.cache.albedo_diffuse))[1] == results[2]

        end
    end
    @testset "reinit_model_sims! for FT=$FT" begin
        @test FieldExchanger.reinit_model_sims!((; stub = TestSurfaceSimulation1(FT(0)))) === nothing
    end

    @testset "step_model_sims! for FT=$FT" begin
        @test FieldExchanger.step_model_sims!((; stub = TestSurfaceSimulation1(FT(0))), 1) === nothing
    end
end
