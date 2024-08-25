#=
    Unit tests for ClimaCoupler ConservationChecker, with parsed objects mimicking those in the full coupled system
=#
import Test: @test, @testset
import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends
import ClimaCore as CC
import ClimaCoupler: ConservationChecker, Interfacer

include("TestHelper.jl")
import .TestHelper

REGRID_DIR = @isdefined(REGRID_DIR) ? REGRID_DIR : joinpath("", "regrid_tmp/")

get_slab_energy(slab_sim, T_sfc) =
    slab_sim.integrator.p.params.ρ .* slab_sim.integrator.p.params.c .* T_sfc .* slab_sim.integrator.p.params.h

struct TestAtmos{I} <: Interfacer.AtmosModelSimulation
    i::I
end
Interfacer.name(s::TestAtmos) = "TestAtmos"
Interfacer.get_field(s::TestAtmos, ::Val{:radiative_energy_flux_toa}) = ones(s.i.space) .* 200
Interfacer.get_field(s::TestAtmos, ::Val{:water}) = ones(s.i.space) .* 1
function Interfacer.get_field(s::TestAtmos, ::Val{:energy})
    FT = CC.Domains.float_type(CC.Meshes.domain(s.i.space.grid.topology.mesh))
    ones(s.i.space) .* FT(1e6)
end

struct TestOcean{I} <: Interfacer.SurfaceModelSimulation
    i::I
end
Interfacer.name(s::TestOcean) = "TestOcean"
Interfacer.get_field(s::TestOcean, ::Val{:water}) = ones(s.i.space) .* 0
function Interfacer.get_field(s::TestOcean, ::Val{:energy})
    FT = CC.Domains.float_type(CC.Meshes.domain(s.i.space.grid.topology.mesh))
    ones(s.i.space) .* FT(1e6)
end
function Interfacer.get_field(s::TestOcean, ::Val{:area_fraction})
    FT = CC.Domains.float_type(CC.Meshes.domain(s.i.space.grid.topology.mesh))
    ones(s.i.space) .* FT(0.25)
end

struct TestLand{I} <: Interfacer.SurfaceModelSimulation
    i::I
end
Interfacer.name(s::TestLand) = "TestLand"
Interfacer.get_field(s::TestLand, ::Val{:energy}) = ones(s.i.space) .* 0
Interfacer.get_field(s::TestLand, ::Val{:water}) = ones(s.i.space) .* 0
function Interfacer.get_field(s::TestLand, ::Val{:area_fraction})
    FT = CC.Domains.float_type(CC.Meshes.domain(s.i.space.grid.topology.mesh))
    ones(s.i.space) .* FT(0.25)
end


for FT in (Float32, Float64)
    @testset "test check_conservation for conservation for FT=$FT" begin
        space = TestHelper.create_space(FT)

        # set up model simulations
        atmos = TestAtmos((; space = space))
        land = TestOcean((; space = space))
        ocean = TestLand((; space = space))
        ice = Interfacer.SurfaceStub((; area_fraction = CC.Fields.ones(space) .* FT(0.5)))
        model_sims = (; atmos_sim = atmos, land_sim = land, ocean_sim = ocean, ice_sim = ice)

        # conservation checkers
        cc = (;
            energy = ConservationChecker.EnergyConservationCheck(model_sims),
            water = ConservationChecker.WaterConservationCheck(model_sims),
        )

        # coupler fields
        cf = (;
            radiative_energy_flux_toa = CC.Fields.ones(space),
            P_net = CC.Fields.zeros(space),
            P_liq = CC.Fields.zeros(space),
            P_snow = CC.Fields.zeros(space),
            F_turb_moisture = CC.Fields.zeros(space),
        )
        @. cf.radiative_energy_flux_toa = 200
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
            nothing, # turbulent_fluxes
            nothing, # thermo_params
        )

        # set non-zero radiation and precipitation
        F_r = cf.radiative_energy_flux_toa
        P = cf.P_liq
        Δt = cs.Δt_cpl

        # analytical solution
        tot_energy_an = sum(FT.(F_r .* 3Δt .+ 1e6 .* 1.25))
        tot_water_an = sum(FT.(.-P .* 3Δt .* 0.5 .+ CC.Fields.ones(space)))

        # run check_conservation!
        ConservationChecker.check_conservation!(cs, runtime_check = true)
        ConservationChecker.check_conservation!(cs, runtime_check = true)
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
                @test length(ccs.total) == length(getfield(ccs, Symbol(Interfacer.name(sim))))
            end
        end

    end

    @testset "test plot_global_conservation with dummy models for FT=$FT" begin
        space = TestHelper.create_space(FT)

        # set up model simulations
        atmos = TestAtmos((; space = space))
        land = TestOcean((; space = space))
        ocean = TestLand((; space = space))
        ice = Interfacer.SurfaceStub((; area_fraction = CC.Fields.ones(space) .* FT(0.5)))
        model_sims = (; atmos_sim = atmos, land_sim = land, ocean_sim = ocean, ice_sim = ice)

        # conservation checkers
        cc = (;
            energy = ConservationChecker.EnergyConservationCheck(model_sims),
            water = ConservationChecker.WaterConservationCheck(model_sims),
        )

        # coupler fields
        cf = (;
            radiative_energy_flux_toa = CC.Fields.ones(space),
            P_net = CC.Fields.zeros(space),
            P_liq = CC.Fields.zeros(space),
            P_snow = CC.Fields.zeros(space),
            F_turb_moisture = CC.Fields.zeros(space),
        )
        @. cf.radiative_energy_flux_toa = 200
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
            nothing, # turbulent_fluxes
            nothing, # thermo_params
        )

        tot_energy, tot_water = ConservationChecker.check_conservation!(cs)
    end
end
