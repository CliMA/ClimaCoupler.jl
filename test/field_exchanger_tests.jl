import Test: @test, @testset
import ClimaCore as CC
import ClimaCoupler: Interfacer, FieldExchanger, FluxCalculator
import Thermodynamics.Parameters as TDP
import ClimaParams # to load TDP extension

# surface field exchange tests
struct TestSurfaceSimulation1{C} <: Interfacer.SurfaceModelSimulation
    cache_field::C
end
struct TestSurfaceSimulation2{C} <: Interfacer.SurfaceModelSimulation
    cache_field::C
end

Interfacer.get_field(
    sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2},
    ::Val{:surface_temperature},
) = sim.cache_field
Interfacer.get_field(
    sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2},
    ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}},
) = sim.cache_field
Interfacer.get_field(
    sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2},
    ::Val{:roughness_momentum},
) = sim.cache_field
Interfacer.get_field(
    sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2},
    ::Val{:roughness_buoyancy},
) = sim.cache_field
Interfacer.get_field(
    sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2},
    ::Val{:beta},
) = sim.cache_field
Interfacer.get_field(
    sim::Union{TestSurfaceSimulation1, TestSurfaceSimulation2},
    ::Val{:emissivity},
) = eltype(sim.cache_field)(1)

Interfacer.get_field(sim::TestSurfaceSimulation1, ::Val{:area_fraction}) =
    sim.cache_field .* CC.Spaces.undertype(axes(sim.cache_field))(0.25)
Interfacer.get_field(sim::TestSurfaceSimulation2, ::Val{:area_fraction}) =
    sim.cache_field .* CC.Spaces.undertype(axes(sim.cache_field))(0.75)

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
Interfacer.get_field(sim::DummyStub, ::Val{:ice_concentration}) = sim.cache.area_fraction
function Interfacer.update_field!(
    sim::DummyStub,
    ::Val{:area_fraction},
    field::CC.Fields.Field,
)
    sim.cache.area_fraction .= field
end
# atmos sim
struct TestAtmosSimulation{C} <: Interfacer.AtmosModelSimulation
    cache::C
end

Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:air_temperature}) =
    sim.cache.air_temperature
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:air_density}) = sim.cache.air_density
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:specific_humidity}) =
    sim.cache.specific_humidity
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:turbulent_energy_flux}) =
    sim.cache.turbulent_energy_flux
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:turbulent_moisture_flux}) =
    sim.cache.turbulent_moisture_flux
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:SW_d}) = sim.cache.SW_d
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:LW_d}) = sim.cache.LW_d
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:liquid_precipitation}) =
    sim.cache.liquid_precipitation
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:snow_precipitation}) =
    sim.cache.snow_precipitation
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:height_int}) =
    CC.Fields.ones(axes(sim.cache.air_temperature))
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:height_sfc}) =
    CC.Fields.zeros(axes(sim.cache.air_temperature))
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:u_int}) =
    CC.Fields.zeros(axes(sim.cache.air_temperature))
Interfacer.get_field(sim::TestAtmosSimulation, ::Val{:v_int}) =
    CC.Fields.zeros(axes(sim.cache.air_temperature))
function Interfacer.update_field!(
    sim::TestAtmosSimulation,
    ::Val{:surface_direct_albedo},
    field,
)
    parent(sim.cache.albedo_direct) .= parent(field)
end
function Interfacer.update_field!(
    sim::TestAtmosSimulation,
    ::Val{:surface_diffuse_albedo},
    field,
)
    parent(sim.cache.albedo_diffuse) .= parent(field)
end
function Interfacer.update_field!(
    sim::TestAtmosSimulation,
    ::Val{:roughness_momentum},
    field,
)
    parent(sim.cache.roughness_momentum) .= parent(field)
end
Interfacer.update_field!(sim::TestAtmosSimulation, ::Val{:surface_temperature}, field) =
    nothing
Interfacer.update_field!(sim::TestAtmosSimulation, ::Val{:roughness_buoyancy}, field) =
    nothing
Interfacer.update_field!(sim::TestAtmosSimulation, ::Val{:beta}, field) = nothing
FluxCalculator.update_turbulent_fluxes!(sim::TestAtmosSimulation, fields) = nothing
Interfacer.step!(sim::TestAtmosSimulation, t) = nothing

#surface sim
struct TestSurfaceSimulationLand{C} <: Interfacer.SurfaceModelSimulation
    cache::C
end
function Interfacer.get_field(sim::TestSurfaceSimulationLand, ::Val{:area_fraction})
    return CC.Fields.ones(axes(sim.cache.liquid_precipitation))
end
Interfacer.get_field(sim::TestSurfaceSimulationLand, ::Val{:surface_direct_albedo}) =
    sim.cache.albedo_direct
Interfacer.get_field(sim::TestSurfaceSimulationLand, ::Val{:surface_diffuse_albedo}) =
    sim.cache.albedo_diffuse
Interfacer.get_field(sim::TestSurfaceSimulationLand, ::Val{:surface_temperature}) =
    sim.cache.surface_temperature
Interfacer.get_field(sim::TestSurfaceSimulationLand, ::Val{:emissivity}) =
    eltype(sim.cache.surface_temperature)(1)
function Interfacer.update_field!(
    sim::TestSurfaceSimulationLand,
    ::Val{:turbulent_energy_flux},
    field,
)
    parent(sim.cache.turbulent_energy_flux) .= parent(field)
end
function Interfacer.update_field!(
    sim::TestSurfaceSimulationLand,
    ::Val{:turbulent_moisture_flux},
    field,
)
    parent(sim.cache.turbulent_moisture_flux) .= parent(field)
end
function Interfacer.update_field!(
    sim::TestSurfaceSimulationLand,
    ::Val{:air_density},
    field,
)
    parent(sim.cache.air_density) .= parent(field)
end
Interfacer.update_field!(sim::TestSurfaceSimulationLand, ::Val{:SW_d}, field) = nothing
Interfacer.update_field!(sim::TestSurfaceSimulationLand, ::Val{:LW_d}, field) = nothing
function Interfacer.update_field!(
    sim::TestSurfaceSimulationLand,
    ::Val{:liquid_precipitation},
    field,
)
    parent(sim.cache.liquid_precipitation) .= parent(field)
end
function Interfacer.update_field!(
    sim::TestSurfaceSimulationLand,
    ::Val{:snow_precipitation},
    field,
)
    parent(sim.cache.snow_precipitation) .= parent(field)
end

for FT in (Float32, Float64)
    @testset "test update_surface_fractions!" begin
        test_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6.371e6), # in meters
            n_quad_points = 4,
            h_elem = 4,
        )
        # Construct land fraction of 0s in first half, 1s in second half
        land_fraction = CC.Fields.ones(test_space)
        nelements = length(CC.Fields.field2array(land_fraction))
        CC.Fields.field2array(land_fraction)[1:(nelements ÷ 2)] .= FT(0)

        # Construct ice fraction of 0.5s on first half, 0s on second half
        ice_d = CC.Fields.zeros(test_space)
        CC.Fields.field2array(ice_d)[1:(nelements ÷ 2)] .= FT(0.5)

        # Construct ice fraction of 0s. During the update, the first half should be set to 0.5
        ocean_d = CC.Fields.zeros(test_space)


        # Fill in only the necessary parts of the simulation
        coupler_fields = Interfacer.init_coupler_fields(
            FT,
            Interfacer.default_coupler_fields(),
            test_space,
        )
        cs = Interfacer.CoupledSimulation{FT}(
            nothing, # dates
            coupler_fields, # fields
            nothing, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # Δt_cpl
            Ref(Int(0)), # t
            Ref(-1), # prev_checkpoint_t
            (;
                ice_sim = DummyStub((; area_fraction = ice_d)),
                ocean_sim = Interfacer.SurfaceStub((; area_fraction = ocean_d)),
                land_sim = DummyStub((; area_fraction = land_fraction)),
            ), # model_sims
            (;), # callbacks
            (;), # dir_paths
            nothing, # thermo_params
            nothing, # diags_handler
        )

        FieldExchanger.update_surface_fractions!(cs)

        # Test that sum of fractions is 1 everywhere
        ice_fraction = Interfacer.get_field(cs.model_sims.ice_sim, Val(:area_fraction))
        ocean_fraction = Interfacer.get_field(cs.model_sims.ocean_sim, Val(:area_fraction))
        @test ice_fraction .+ ocean_fraction .+ land_fraction == ones(test_space)
        # Test that fractions updated correctly (ice and land unchanged, ocean updated)
        @test all(CC.Fields.field2array(ice_fraction)[1:(nelements ÷ 2)] .== FT(0.5))
        @test all(CC.Fields.field2array(ocean_fraction)[1:(nelements ÷ 2)] .== FT(0.5))
    end

    @testset "test combine_surfaces" begin
        test_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6.371e6), # in meters
            n_quad_points = 4,
            h_elem = 4,
        )
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
        temp_field = CC.Fields.zeros(test_space)

        FieldExchanger.combine_surfaces!(combined_field, sims, var_name, temp_field)
        @test combined_field == fill(FT(sum(fractions .* fields)), test_space)
    end

    @testset "import_atmos_fields! for FT=$FT" begin
        boundary_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6.371e6), # in meters
            n_quad_points = 4,
            h_elem = 4,
        )
        # Initialize coupler fields with 0
        coupler_names = Interfacer.default_coupler_fields()
        coupler_fields = Interfacer.init_coupler_fields(FT, coupler_names, boundary_space)
        # Initialize component fields with 1
        component_names = [
            :air_density,
            :air_temperature,
            :specific_humidity,
            :SW_d,
            :LW_d,
            :liquid_precipitation,
            :snow_precipitation,
        ]
        key_types = (component_names...,)
        val_types = Tuple{(FT for _ in 1:length(component_names))...}
        component_fields = ones(NamedTuple{key_types, val_types}, boundary_space)

        model_sims = (;
            atmos_sim = TestAtmosSimulation(component_fields),
            land_sim = TestSurfaceSimulation1(component_fields),
        )

        FieldExchanger.import_atmos_fields!(coupler_fields, model_sims)
        zero_field = zeros(boundary_space)
        @test coupler_fields.F_lh == zero_field
        @test coupler_fields.F_sh == zero_field
        @test coupler_fields.F_turb_moisture == zero_field
        @test coupler_fields.SW_d == model_sims.atmos_sim.cache.SW_d
        @test coupler_fields.LW_d == model_sims.atmos_sim.cache.LW_d
        @test coupler_fields.P_liq == model_sims.atmos_sim.cache.liquid_precipitation
        @test coupler_fields.P_snow == model_sims.atmos_sim.cache.snow_precipitation

    end

    @testset "import_combined_surface_fields! for FT=$FT" begin
        # coupler cache setup
        boundary_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6.371e6), # in meters
            n_quad_points = 4,
            h_elem = 4,
        )
        coupler_names_additional = [:surface_direct_albedo, :surface_diffuse_albedo]
        coupler_names =
            push!(Interfacer.default_coupler_fields(), coupler_names_additional...)
        coupler_fields = Interfacer.init_coupler_fields(FT, coupler_names, boundary_space)

        sims = (;
            atmos_sim = TestAtmosSimulation((; air_temperature = ones(boundary_space),)),
            a = TestSurfaceSimulation1(ones(boundary_space)),
            b = TestSurfaceSimulation2(ones(boundary_space)),
        )

        FieldExchanger.import_static_fields!(coupler_fields, sims)
        FieldExchanger.import_combined_surface_fields!(coupler_fields, sims)

        # Analytically compute expected values and compare
        expected_field_temp =
            (
                Interfacer.get_field(sims.a, Val(:area_fraction)) .*
                Interfacer.get_field(sims.a, Val(:emissivity)) .*
                sims.a.cache_field .^ FT(4) .+
                Interfacer.get_field(sims.b, Val(:area_fraction)) .*
                Interfacer.get_field(sims.b, Val(:emissivity)) .*
                sims.b.cache_field .^ FT(4)
            ) .^ FT(1 / 4)
        @test coupler_fields.T_sfc == expected_field_temp

        expected_field =
            Interfacer.get_field(sims.a, Val(:area_fraction)) .* sims.a.cache_field .+
            Interfacer.get_field(sims.b, Val(:area_fraction)) .* sims.b.cache_field
        @test coupler_fields.surface_direct_albedo == expected_field
        @test coupler_fields.surface_diffuse_albedo == expected_field
    end

    @testset "update_model_sims! for FT=$FT" begin
        # coupler cache setup
        boundary_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6.371e6), # in meters
            n_quad_points = 4,
            h_elem = 4,
        )
        coupler_field_names = [
            Interfacer.default_coupler_fields()
            [:surface_direct_albedo, :surface_diffuse_albedo]
        ]
        coupler_fields =
            Interfacer.init_coupler_fields(FT, coupler_field_names, boundary_space)
        # Initialize with ones
        for p in propertynames(coupler_fields)
            fill!(getproperty(coupler_fields, p), 1)
        end

        # model cache setup

        atmos_names = [
            :surface_temperature,
            :albedo_direct,
            :albedo_diffuse,
            :roughness_momentum,
            :roughness_buoyancy,
            :beta,
        ]
        atmos_fields = Interfacer.init_coupler_fields(FT, atmos_names, boundary_space)

        land_names = [
            :turbulent_energy_flux,
            :turbulent_moisture_flux,
            :air_density,
            :SW_d,
            :LW_d,
            :liquid_precipitation,
            :snow_precipitation,
        ]
        land_fields = Interfacer.init_coupler_fields(FT, land_names, boundary_space)

        model_sims = (;
            atmos_sim = TestAtmosSimulation(atmos_fields),
            land_sim = TestSurfaceSimulationLand(land_fields),
            stub_sim = Interfacer.SurfaceStub((;
                area_fraction = CC.Fields.ones(boundary_space),
                albedo_direct = CC.Fields.ones(boundary_space),
                albedo_diffuse = CC.Fields.ones(boundary_space),
            )),
        )
        coupler_fields.surface_diffuse_albedo .= FT(0.5)

        # test the sim update
        results = [FT(0), FT(1), FT(0.5)]
        model_sims.atmos_sim.cache.roughness_momentum .= FT(0)
        FieldExchanger.update_model_sims!(model_sims, coupler_fields)

        # test atmos
        expected_field = fill(results[2], boundary_space)
        @test model_sims.atmos_sim.cache.albedo_direct == expected_field
        expected_field .= results[3]
        @test model_sims.atmos_sim.cache.albedo_diffuse == expected_field

        # test variables without updates
        expected_field .= results[1]
        @test model_sims.atmos_sim.cache.surface_temperature == expected_field
        @test model_sims.atmos_sim.cache.beta == expected_field
        @test model_sims.atmos_sim.cache.roughness_buoyancy == expected_field

        # test land updates
        expected_field .= results[2]
        @test model_sims.land_sim.cache.liquid_precipitation == expected_field
        @test model_sims.land_sim.cache.snow_precipitation == expected_field

        # test variables without updates
        expected_field .= results[1]
        @test model_sims.land_sim.cache.SW_d == expected_field
        @test model_sims.land_sim.cache.LW_d == expected_field

        # test stub - albedo should be updated by update_sim!
        expected_field .= results[2]
        @test model_sims.stub_sim.cache.albedo_direct == expected_field
        @test model_sims.stub_sim.cache.albedo_diffuse == expected_field
    end

    @testset "step_model_sims! for FT=$FT" begin
        @test FieldExchanger.step_model_sims!(
            (;
                stub = TestSurfaceSimulation1(FT(0)),
                atmos_sim = TestAtmosSimulation(FT(0)),
            ),
            1,
            nothing,
            nothing,
        ) === nothing
    end

    @testset "exchange! for FT=$FT" begin
        # Here we exchange two albedo fields from the surface to the atmos,
        #  and two precipitation fields from the atmos to the surface.
        # coupler cache setup
        boundary_space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6.371e6), # in meters
            n_quad_points = 4,
            h_elem = 4,
        )
        coupler_field_names = [
            Interfacer.default_coupler_fields()
            [:surface_direct_albedo, :surface_diffuse_albedo]
        ]
        # Initialize coupler fields with 0.5
        key_types = (coupler_field_names...,)
        val_types = Tuple{(FT for _ in 1:length(coupler_field_names))...}
        coupler_fields = zeros(NamedTuple{key_types, val_types}, boundary_space) .+ FT(0.5)

        # model cache setup
        atmos_names = [
            :albedo_direct,
            :albedo_diffuse,
            :liquid_precipitation,
            :snow_precipitation,
            :SW_d,
            :LW_d,
            :air_temperature,
            :specific_humidity,
            :air_density,
        ]
        # Initialize atmos fields with 1
        key_types = (atmos_names...,)
        val_types = Tuple{(FT for _ in 1:length(atmos_names))...}
        atmos_fields = ones(NamedTuple{key_types, val_types}, boundary_space)

        land_names = [
            :albedo_direct,
            :albedo_diffuse,
            :liquid_precipitation,
            :snow_precipitation,
            :surface_temperature,
            :air_density,
        ]
        # Here we initialize with `init_coupler_fields` because we want fields of zeros
        land_fields = Interfacer.init_coupler_fields(FT, land_names, boundary_space)

        model_sims = (;
            atmos_sim = TestAtmosSimulation(atmos_fields),
            land_sim = TestSurfaceSimulationLand(land_fields),
        )
        thermo_params = TDP.ThermodynamicsParameters(FT)

        # construct the CoupledSimulation object
        cs = Interfacer.CoupledSimulation{FT}(
            nothing, # start_date
            coupler_fields,
            nothing, # conservation_checks
            nothing, # tspan
            nothing, # dt
            nothing, # t
            Ref(-1), # prev_checkpoint_t
            model_sims,
            (;), # callbacks
            (;), # dir_paths
            thermo_params, # thermo_params
            nothing, # diags_handler
        )

        # perform the exchange
        FieldExchanger.exchange!(cs)

        surface_init_field = CC.Fields.zeros(boundary_space)
        atmos_init_field = CC.Fields.ones(boundary_space)

        # test atmos is updated with surface fields
        @test model_sims.atmos_sim.cache.albedo_direct == surface_init_field
        @test model_sims.atmos_sim.cache.albedo_diffuse == surface_init_field

        # test land is updated with atmos fields
        @test model_sims.land_sim.cache.liquid_precipitation == atmos_init_field
        @test model_sims.land_sim.cache.snow_precipitation == atmos_init_field

        # test coupler fields were updated by atmos and surface
        @test cs.fields.surface_direct_albedo == surface_init_field
        @test cs.fields.surface_diffuse_albedo == surface_init_field
        @test cs.fields.SW_d == atmos_init_field
        @test cs.fields.LW_d == atmos_init_field
        @test cs.fields.P_liq == atmos_init_field
        @test cs.fields.P_snow == atmos_init_field
    end
end
