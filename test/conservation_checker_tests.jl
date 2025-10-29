#=
    Unit tests for ClimaCoupler ConservationChecker, with parsed objects mimicking those in the full coupled system
=#
import Test: @test, @testset
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC
import ClimaCoupler: ConservationChecker, Interfacer
import ClimaCoupler.Utilities: integral

get_slab_energy(slab_sim, T_sfc) =
    slab_sim.integrator.p.params.ρ .* slab_sim.integrator.p.params.c .* T_sfc .*
    slab_sim.integrator.p.params.h

struct TestAtmos{I} <: Interfacer.AtmosModelSimulation
    i::I
end
Interfacer.get_field(s::TestAtmos, ::Val{:radiative_energy_flux_toa}) =
    ones(s.i.space) .* 200
Interfacer.get_field(s::TestAtmos, ::Val{:water}) = s.i.water
Interfacer.get_field(s::TestAtmos, ::Val{:energy}) = s.i.energy

struct TestOcean{I} <: Interfacer.SurfaceModelSimulation
    i::I
end
Interfacer.get_field(s::TestOcean, ::Val{:water}) = zeros(s.i.space)
Interfacer.get_field(s::TestOcean, ::Val{:energy}) =
    ones(s.i.space) .* CC.Spaces.undertype(s.i.space)(1e6)
Interfacer.get_field(s::TestOcean, ::Val{:area_fraction}) =
    ones(s.i.space) .* CC.Spaces.undertype(s.i.space)(0.25)

struct TestLand{I} <: Interfacer.SurfaceModelSimulation
    i::I
end
Interfacer.get_field(s::TestLand, ::Val{:energy}) = zeros(s.i.space)
Interfacer.get_field(s::TestLand, ::Val{:water}) = zeros(s.i.space)
Interfacer.get_field(s::TestLand, ::Val{:area_fraction}) =
    ones(s.i.space) .* CC.Spaces.undertype(s.i.space)(0.25)

for FT in (Float32, Float64)
    @testset "test check_conservation for conservation for FT=$FT" begin
        space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6371e3),
            n_quad_points = 4,
            h_elem = 4,
        )

        # set up model simulations
        initial_energy = ones(space) .* CC.Spaces.undertype(space)(1e6)
        initial_water = ones(space) .* CC.Spaces.undertype(space)(1e5)

        atmos = TestAtmos((; space = space, energy = initial_energy, water = initial_water))
        land = TestOcean((; space = space))
        ocean = TestLand((; space = space))
        ice = Interfacer.SurfaceStub((; area_fraction = CC.Fields.ones(space) .* FT(0.5)))
        model_sims =
            (; atmos_sim = atmos, land_sim = land, ocean_sim = ocean, ice_sim = ice)

        # conservation checkers
        cc = (;
            energy = ConservationChecker.EnergyConservationCheck(model_sims),
            water = ConservationChecker.WaterConservationCheck(model_sims),
        )

        # coupler fields
        cf = Interfacer.init_coupler_fields(
            FT,
            [:P_net, :P_liq, :P_snow, :F_turb_moisture],
            space,
        )
        @. cf.P_liq = -100

        # init
        cs = Interfacer.CoupledSimulation{FT}(
            nothing, # dates
            cf, # fields
            cc, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # Δt_cpl
            Ref(Int(0)), # t
            Ref(-1), # prev_checkpoint_t
            model_sims, # model_sims
            (;), # callbacks
            (;), # dir_paths
            nothing, # thermo_params
            nothing, # diags_handler
        )

        # set non-zero radiation and precipitation
        F_r = 200 .* CC.Fields.ones(space)
        P = cf.P_liq
        Δt = float(cs.Δt_cpl)

        volume = integral(ones(space))

        area_fraction_scaling =
            Interfacer.get_field(land, Val(:area_fraction)) .+
            Interfacer.get_field(ocean, Val(:area_fraction))
        water_from_precipitation = integral(P .* area_fraction_scaling) .* FT(Δt)
        energy_from_radiation = integral(F_r) .* FT(Δt)
        energy_per_unit_cell = CC.Fields.ones(space) .* energy_from_radiation ./ volume
        water_per_unit_cell = CC.Fields.ones(space) .* water_from_precipitation ./ volume

        # analytical solution
        # Only ocean and atmos have energy
        area_fraction_scaling =
            CC.Fields.ones(space) .+ Interfacer.get_field(ocean, Val(:area_fraction))
        tot_energy_an =
            integral(area_fraction_scaling .* FT.(initial_energy)) + energy_from_radiation
        tot_water_an = integral(FT.(initial_water)) - water_from_precipitation

        # run check_conservation!
        ConservationChecker.check_conservation!(cs, runtime_check = true)
        atmos.i.energy .-= energy_per_unit_cell
        atmos.i.water .+= water_per_unit_cell
        ConservationChecker.check_conservation!(cs, runtime_check = true)
        atmos.i.energy .-= energy_per_unit_cell
        atmos.i.water .+= water_per_unit_cell
        ConservationChecker.check_conservation!(cs, runtime_check = true)

        total_energy = cs.conservation_checks.energy.sums.total
        total_water = cs.conservation_checks.water.sums.total

        @test total_energy[end] ==
              (cc.energy.sums.TestAtmos + cc.energy.sums.TestOcean + cc.energy.sums.TestLand + cc.energy.sums.SurfaceStub + cc.energy.sums.toa_net_source)[end]
        @test total_water[end] ==
              (cc.water.sums.TestAtmos + cc.water.sums.TestOcean + cc.water.sums.TestLand + cc.water.sums.SurfaceStub)[end]

        @test abs((total_energy[end] .- tot_energy_an) / total_energy[end]) < 1e-4
        @test abs((total_water[end] .- tot_water_an) / total_water[end]) < 1e-4

        # check that all the totals in each checker have the same length
        for sim in values(model_sims)
            for checker in values(cc)
                ccs = checker.sums
                @test length(ccs.total) == length(getfield(ccs, Symbol(nameof(sim))))
            end
        end

    end

    @testset "test plot_global_conservation with dummy models for FT=$FT" begin
        space = CC.CommonSpaces.CubedSphereSpace(
            FT;
            radius = FT(6371e3),
            n_quad_points = 4,
            h_elem = 4,
        )

        # set up model simulations
        initial_energy = ones(space) .* CC.Spaces.undertype(space)(1e6)
        initial_water = ones(space) .* CC.Spaces.undertype(space)(1e5)

        atmos = TestAtmos((; space = space, energy = initial_energy, water = initial_water))
        land = TestOcean((; space = space))
        ocean = TestLand((; space = space))
        ice = Interfacer.SurfaceStub((; area_fraction = CC.Fields.ones(space) .* FT(0.5)))
        model_sims =
            (; atmos_sim = atmos, land_sim = land, ocean_sim = ocean, ice_sim = ice)

        # conservation checkers
        cc = (;
            energy = ConservationChecker.EnergyConservationCheck(model_sims),
            water = ConservationChecker.WaterConservationCheck(model_sims),
        )

        # coupler fields
        cf = Interfacer.init_coupler_fields(
            FT,
            [:P_net, :P_liq, :P_snow, :F_turb_moisture],
            space,
        )
        @. cf.P_liq = -100

        # init
        cs = Interfacer.CoupledSimulation{FT}(
            nothing, # dates
            cf, # fields
            cc, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # Δt_cpl
            Ref(Int(0)), # t
            Ref(-1), # prev_checkpoint_t
            model_sims, # model_sims
            (;), # callbacks
            (;), # dir_paths
            nothing, # thermo_params
            nothing, # diags_handler
        )

        tot_energy, tot_water = ConservationChecker.check_conservation!(cs)
    end
end
