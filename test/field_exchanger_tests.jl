using Test
using ClimaCore: Meshes, Domains, Topologies, Spaces, Fields, InputOutput
using ClimaCoupler: Utilities, Regridder, TestHelper
import ClimaCoupler.Interfacer:
    update_field!, AtmosModelSimulation, SurfaceModelSimulation, SurfaceStub, get_field, update_field!
import ClimaCoupler.FluxCalculator: CombinedStateFluxes, PartitionedStateFluxes, calculate_surface_air_density
import ClimaCoupler.FieldExchanger:
    import_atmos_fields!,
    import_combined_surface_fields!,
    update_model_sims!,
    reinit_model_sims!,
    reinit!,
    step_model_sims!,
    step!


# test for a simple generic atmos model
struct DummySimulation{C} <: AtmosModelSimulation
    cache::C
end

get_field(sim::DummySimulation, ::Val{:turbulent_energy_flux}) = sim.cache.turbulent_energy_flux
get_field(sim::DummySimulation, ::Val{:turbulent_moisture_flux}) = sim.cache.turbulent_moisture_flux
get_field(sim::DummySimulation, ::Val{:radiative_energy_flux_sfc}) = sim.cache.radiative_energy_flux_sfc
get_field(sim::DummySimulation, ::Val{:liquid_precipitation}) = sim.cache.liquid_precipitation
get_field(sim::DummySimulation, ::Val{:snow_precipitation}) = sim.cache.snow_precipitation

function calculate_surface_air_density(atmos_sim::DummySimulation, T_S::Fields.Field)
    FT = eltype(T_S)
    return T_S .* FT(0.0) .+ FT(1.0)
end


# surface field exchange tests
struct TestSurfaceSimulation1{C} <: SurfaceModelSimulation
    cache_field::C
end
struct TestSurfaceSimulation2{C} <: SurfaceModelSimulation
    cache_field::C
end

get_field(sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2}, ::Val{:surface_temperature}) =
    sim.cache_field .* eltype(sim.cache_field)(1.0)
get_field(
    sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2},
    ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}},
) = sim.cache_field .* eltype(sim.cache_field)(1.0)
get_field(sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2}, ::Val{:roughness_momentum}) =
    sim.cache_field .* eltype(sim.cache_field)(1.0)
get_field(sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2}, ::Val{:roughness_buoyancy}) =
    sim.cache_field .* eltype(sim.cache_field)(1.0)
get_field(sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2}, ::Val{:beta}) =
    sim.cache_field .* eltype(sim.cache_field)(1.0)

get_field(sim::TestSurfaceSimulation1, ::Val{:area_fraction}) = sim.cache_field .* eltype(sim.cache_field)(0)
get_field(sim::TestSurfaceSimulation2, ::Val{:area_fraction}) = sim.cache_field .* eltype(sim.cache_field)(0.5)

get_field(sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2}, ::Val{:surface_humidity}) =
    sim.cache_field .* eltype(sim.cache_field)(0)
get_field(sim::Union{TestSurfaceSimulation2, TestSurfaceSimulation2}, ::Val{:surface_humidity}) =
    sim.cache_field .* eltype(sim.cache_field)(0)


# atmos sim
struct TestAtmosSimulation{C} <: AtmosModelSimulation
    cache::C
end
function update_field!(sim::TestAtmosSimulation, ::Val{:surface_direct_albedo}, field)
    parent(sim.cache.albedo_direct) .= parent(field)
end
function update_field!(sim::TestAtmosSimulation, ::Val{:surface_diffuse_albedo}, field)
    parent(sim.cache.albedo_diffuse) .= parent(field)
end
function update_field!(sim::TestAtmosSimulation, ::Val{:roughness_momentum}, field)
    parent(sim.cache.roughness_momentum) .= parent(field)
end

update_field!(sim::TestAtmosSimulation, ::Val{:roughness_buoyancy}, field) = nothing
update_field!(sim::TestAtmosSimulation, ::Val{:beta}, field) = nothing

#surface sim
struct TestSurfaceSimulationLand{C} <: SurfaceModelSimulation
    cache::C
end
function get_field(sim::TestSurfaceSimulationLand, ::Val{:area_fraction})
    FT = eltype(sim.cache.turbulent_energy_flux)
    return FT(0.5)
end
function update_field!(sim::TestSurfaceSimulationLand, ::Val{:turbulent_energy_flux}, field)
    parent(sim.cache.turbulent_energy_flux) .= parent(field)
end
function update_field!(sim::TestSurfaceSimulationLand, ::Val{:turbulent_moisture_flux}, field)
    parent(sim.cache.turbulent_moisture_flux) .= parent(field)
end

for FT in (Float32, Float64)
    @testset "import_atmos_fields! for FT=$FT" begin
        boundary_space = TestHelper.create_space(FT)
        coupler_names = (:F_turb_energy, :F_turb_moisture, :F_radiative, :P_liq, :P_snow, :ρ_sfc, :T_S)
        atmos_names = (
            :turbulent_energy_flux,
            :turbulent_moisture_flux,
            :radiative_energy_flux_sfc,
            :liquid_precipitation,
            :snow_precipitation,
        )
        atmos_fields = NamedTuple{atmos_names}(ntuple(i -> Fields.ones(boundary_space), length(atmos_names)))

        model_sims = (; atmos_sim = DummySimulation(atmos_fields))

        flux_types = (CombinedStateFluxes(), PartitionedStateFluxes())
        results = [FT(1), FT(0)]
        for (i, t) in enumerate(flux_types)
            coupler_fields = NamedTuple{coupler_names}(ntuple(i -> Fields.zeros(boundary_space), length(coupler_names)))
            import_atmos_fields!(coupler_fields, model_sims, boundary_space, t)
            @test parent(coupler_fields.F_turb_energy)[1] == results[i]
            @test parent(coupler_fields.F_turb_moisture)[1] == results[i]
            @test parent(coupler_fields.F_radiative)[1] == results[1]
            @test parent(coupler_fields.P_liq)[1] == results[1]
            @test parent(coupler_fields.P_snow)[1] == results[1]
        end
    end

    @testset "import_combined_surface_fields! for FT=$FT" begin
        # coupler cache setup
        boundary_space = TestHelper.create_space(FT)
        coupler_names = (:T_S, :z0m_S, :z0b_S, :surface_direct_albedo, :surface_diffuse_albedo, :beta, :q_sfc)

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

        # test the coupler update under CombinedStateFluxes (update all) and PartitionedStateFluxes (update all except
        # surface info needed for the calculation of turbulent fluxes)
        flux_types = (CombinedStateFluxes(), PartitionedStateFluxes())
        results = [FT(0.5), FT(0)] # 0.5 due to the area fraction weighting
        for (i, t) in enumerate(flux_types)
            coupler_fields = NamedTuple{coupler_names}(ntuple(i -> Fields.zeros(boundary_space), length(coupler_names)))
            import_combined_surface_fields!(coupler_fields, sims, boundary_space, t)
            @test parent(coupler_fields.T_S)[1] == results[1]
            @test parent(coupler_fields.surface_direct_albedo)[1] == results[1]
            @test parent(coupler_fields.surface_diffuse_albedo)[1] == results[1]
            @test parent(coupler_fields.z0m_S)[1] == results[i]
            @test parent(coupler_fields.z0b_S)[1] == results[i]
            @test parent(coupler_fields.beta)[1] == results[i]
        end
    end

    @testset "update_model_sims! for FT=$FT" begin
        # coupler cache setup
        boundary_space = TestHelper.create_space(FT)
        coupler_field_names = (
            :ρ_sfc,
            :T_S,
            :z0m_S,
            :z0b_S,
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
            NamedTuple{coupler_field_names}(ntuple(i -> Fields.ones(boundary_space), length(coupler_field_names)))

        # model cache setup

        atmos_names =
            (:surface_temperature, :albedo_direct, :albedo_diffuse, :roughness_momentum, :roughness_buoyancy, :beta)
        atmos_fields = NamedTuple{atmos_names}(ntuple(i -> Fields.zeros(boundary_space), length(atmos_names)))

        land_names = (
            :turbulent_energy_flux,
            :turbulent_moisture_flux,
            :radiative_energy_flux_sfc,
            :liquid_precipitation,
            :snow_precipitation,
            :ρ_sfc,
        )
        land_fields = NamedTuple{land_names}(ntuple(i -> Fields.zeros(boundary_space), length(land_names)))

        model_sims = (;
            atmos_sim = TestAtmosSimulation(atmos_fields),
            land_sim = TestSurfaceSimulationLand(land_fields),
            stub_sim = SurfaceStub((;
                area_fraction = Fields.ones(boundary_space),
                ρ_sfc = Fields.ones(boundary_space),
                albedo_direct = Fields.ones(boundary_space),
                albedo_diffuse = Fields.ones(boundary_space),
            )),
        )
        coupler_fields.surface_diffuse_albedo .= FT(0.5)

        # test the sim update under CombinedStateFluxes (update all) and PartitionedStateFluxes (update all except turbulent fluxes)
        flux_types = (CombinedStateFluxes(), PartitionedStateFluxes())
        results = [FT(0), FT(1), FT(0.5)]
        for (i, t) in enumerate(flux_types)
            model_sims.atmos_sim.cache.roughness_momentum .= FT(0)
            update_model_sims!(model_sims, coupler_fields, t)

            # test atmos
            @test parent(model_sims.atmos_sim.cache.albedo_direct)[1] == results[2]
            @test parent(model_sims.atmos_sim.cache.albedo_diffuse)[1] == results[3]
            if t isa CombinedStateFluxes
                @test parent(model_sims.atmos_sim.cache.roughness_momentum)[1] == results[2]
            else
                @test parent(model_sims.atmos_sim.cache.roughness_momentum)[1] == results[1]
            end

            # unspecified variables
            @test parent(model_sims.atmos_sim.cache.surface_temperature)[1] == results[1]
            @test parent(model_sims.atmos_sim.cache.beta)[1] == results[1]
            @test parent(model_sims.atmos_sim.cache.roughness_buoyancy)[1] == results[1]

            # test surface
            @test parent(model_sims.land_sim.cache.turbulent_energy_flux)[1] == results[2] # assuming units / m2
            @test parent(model_sims.land_sim.cache.turbulent_moisture_flux)[1] == results[2]

            # unspecified variables
            @test parent(model_sims.land_sim.cache.radiative_energy_flux_sfc)[1] == results[1]
            @test parent(model_sims.land_sim.cache.liquid_precipitation)[1] == results[1]
            @test parent(model_sims.land_sim.cache.snow_precipitation)[1] == results[1]

            # test stub - albedo should be updated by update_sim!
            @test parent(model_sims.stub_sim.cache.albedo_direct)[1] == results[2]
            @test parent(model_sims.stub_sim.cache.albedo_diffuse)[1] == results[2]

        end
    end
    @testset "reinit_model_sims! for FT=$FT" begin
        @test reinit_model_sims!((; stub = SurfaceStub(FT(0)))) == nothing
    end

    @testset "step_model_sims! for FT=$FT" begin
        @test step_model_sims!((; stub = SurfaceStub(FT(0))), 1) == nothing
    end
end
