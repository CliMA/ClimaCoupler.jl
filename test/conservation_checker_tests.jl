#=
    Unit tests for ClimaCoupler ConservationChecker, with parsed objects mimicking those in the full coupled system
=#

using ClimaCoupler: Regridder, TestHelper, Interfacer
using ClimaCoupler.ConservationChecker:
    EnergyConservationCheck, WaterConservationCheck, check_conservation!, plot_global_conservation
using ClimaCore: ClimaCore, Geometry, Meshes, Domains, Topologies, Spaces, Fields, InputOutput
import ClimaCore.InputOutput: read_field
using ClimaLand
using ClimaComms
using Test
using NCDatasets
using Dates
using Downloads

import ClimaCoupler.Interfacer: AtmosModelSimulation, SurfaceModelSimulation, SurfaceStub, get_field, name

REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath("", "regrid_tmp/")

get_slab_energy(slab_sim, T_sfc) =
    slab_sim.integrator.p.params.ρ .* slab_sim.integrator.p.params.c .* T_sfc .* slab_sim.integrator.p.params.h

struct TestAtmos{I} <: Interfacer.AtmosModelSimulation
    i::I
end
name(s::TestAtmos) = "TestAtmos"
get_field(s::TestAtmos, ::Val{:F_radiative_TOA}) = ones(s.i.space) .* 200
get_field(s::TestAtmos, ::Val{:water}) = ones(s.i.space) .* 1
function get_field(s::TestAtmos, ::Val{:energy})
    FT = Domains.float_type(Meshes.domain(s.i.space.grid.topology.mesh))
    ones(s.i.space) .* FT(1e6)
end

struct TestOcean{I} <: Interfacer.SurfaceModelSimulation
    i::I
end
name(s::TestOcean) = "TestOcean"
get_field(s::TestOcean, ::Val{:water}) = ones(s.i.space) .* 0
function get_field(s::TestOcean, ::Val{:energy})
    FT = Domains.float_type(Meshes.domain(s.i.space.grid.topology.mesh))
    ones(s.i.space) .* FT(1e6)
end
function get_field(s::TestOcean, ::Val{:area_fraction})
    FT = Domains.float_type(Meshes.domain(s.i.space.grid.topology.mesh))
    ones(s.i.space) .* FT(0.25)
end

struct TestLand{I} <: Interfacer.SurfaceModelSimulation
    i::I
end
name(s::TestLand) = "TestLand"
get_field(s::TestLand, ::Val{:energy}) = ones(s.i.space) .* 0
get_field(s::TestLand, ::Val{:water}) = ones(s.i.space) .* 0
function get_field(s::TestLand, ::Val{:area_fraction})
    FT = Domains.float_type(Meshes.domain(s.i.space.grid.topology.mesh))
    ones(s.i.space) .* FT(0.25)
end


for FT in (Float32, Float64)
    @testset "test check_conservation for conservation for FT=$FT" begin
        space = TestHelper.create_space(FT)

        # set up model simulations
        atmos = TestAtmos((; space = space))
        land = TestOcean((; space = space))
        ocean = TestLand((; space = space))
        ice = SurfaceStub((; area_fraction = Fields.ones(space) .* FT(0.5)))
        model_sims = (; atmos_sim = atmos, land_sim = land, ocean_sim = ocean, ice_sim = ice)

        # conservation checkers
        cc = (; energy = EnergyConservationCheck(model_sims), water = WaterConservationCheck(model_sims))

        # coupler fields
        cf = (;
            F_radiative_TOA = Fields.ones(space),
            P_net = Fields.zeros(space),
            P_liq = Fields.zeros(space),
            P_snow = Fields.zeros(space),
            F_turb_moisture = Fields.zeros(space),
        )
        @. cf.F_radiative_TOA = 200
        @. cf.P_liq = -100

        # init
        cs = Interfacer.CoupledSimulation{FT}(
            nothing, # comms_ctx
            nothing, # dates
            space, # boundary_space
            cf, # fields
            nothing, # parsed_args
            cc, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # t
            Int(200), # Δt_cpl
            (;), # surface_masks
            model_sims, # model_sims
            (;), # mode
            (), # diagnostics
            (;), # callbacks
            (;), # dirs
        )

        # set non-zero radiation and precipitation
        F_r = cf.F_radiative_TOA
        P = cf.P_liq
        Δt = cs.Δt_cpl

        # analytical solution
        tot_energy_an = sum(FT.(F_r .* 3Δt .+ 1e6 .* 1.25))
        tot_water_an = sum(FT.(.-P .* 3Δt .* 0.5 .+ Fields.ones(space)))

        # run check_conservation!
        check_conservation!(cs, runtime_check = true)
        check_conservation!(cs, runtime_check = true)
        check_conservation!(cs, runtime_check = true)

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
                @test length(ccs.total) == length(getfield(ccs, Symbol(name(sim))))
            end
        end

    end

    @testset "test plot_global_conservation with dummy models for FT=$FT" begin
        space = TestHelper.create_space(FT)

        # set up model simulations
        atmos = TestAtmos((; space = space))
        land = TestOcean((; space = space))
        ocean = TestLand((; space = space))
        ice = SurfaceStub((; area_fraction = Fields.ones(space) .* FT(0.5)))
        model_sims = (; atmos_sim = atmos, land_sim = land, ocean_sim = ocean, ice_sim = ice)

        # conservation checkers
        cc = (; energy = EnergyConservationCheck(model_sims), water = WaterConservationCheck(model_sims))

        # coupler fields
        cf = (;
            F_radiative_TOA = Fields.ones(space),
            P_net = Fields.zeros(space),
            P_liq = Fields.zeros(space),
            P_snow = Fields.zeros(space),
            F_turb_moisture = Fields.zeros(space),
        )
        @. cf.F_radiative_TOA = 200
        @. cf.P_liq = -100

        # init
        cs = Interfacer.CoupledSimulation{FT}(
            nothing, # comms_ctx
            nothing, # dates
            space, # boundary_space
            cf, # fields
            nothing, # parsed_args
            cc, # conservation_checks
            (Int(0), Int(1000)), # tspan
            Int(200), # t
            Int(200), # Δt_cpl
            (;), # surface_masks
            model_sims, # model_sims
            (;), # mode
            (), # diagnostics
            (;), # callbacks
            (;), # dirs
        )

        tot_energy, tot_water = check_conservation!(cs)

        output_plots = "test_cons_plots/"
        mkpath(output_plots)
        plot_global_conservation(
            cs.conservation_checks.energy,
            cs,
            figname1 = output_plots * "energy.png",
            figname2 = output_plots * "energy_log.png",
        )
        plot_global_conservation(
            cs.conservation_checks.water,
            cs,
            figname1 = output_plots * "water.png",
            figname2 = output_plots * "water_log.png",
        )

        # test that files exist
        @test isfile("test_cons_plots/energy.png")
        @test isfile("test_cons_plots/energy_log.png")
        @test isfile("test_cons_plots/water.png")
        @test isfile("test_cons_plots/water_log.png")

        # remove output
        rm(output_plots; recursive = true)
    end
end
