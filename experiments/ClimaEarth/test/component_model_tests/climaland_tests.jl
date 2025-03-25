using Test
import ClimaCore as CC
import ClimaCoupler
import Dates

include(joinpath("..", "TestHelper.jl"))
import .TestHelper
include(joinpath("..", "..", "components", "land", "climaland_integrated.jl"))

FT = Float32

@testset "ClimaLandSimulation constructor" begin
    dt = Float64(450)
    tspan = (Float64(0), 3.0dt)
    start_date = Dates.DateTime(2008)
    output_dir = pwd()
    boundary_space = TestHelper.create_space(FT)
    area_fraction = CC.Fields.ones(boundary_space)

    # Construct simulation object
    land_sim = ClimaLandSimulation(FT, dt, tspan, start_date, output_dir, boundary_space, area_fraction)

    # Try taking a timestep
    Interfacer.step!(land_sim, dt)

    # Check that the simulation object is correctly initialized
    @test Interfacer.name(land_sim) == "ClimaLandSimulation"
    @test land_sim.area_fraction == area_fraction

    # Check that the state is correctly initialized
    state_names = propertynames(land_sim.integrator.u)
    @test :canopy in state_names
    @test :soil in state_names
    @test :snow in state_names
    @test :soilco2 in state_names

    # Check that the cache is correctly initialized
    cache_names = propertynames(land_sim.integrator.p)
    @test :canopy in cache_names
    @test :soil in cache_names
    @test :snow in cache_names
    @test :soilco2 in cache_names
    @test :drivers in cache_names

    # Check that the drivers are correctly initialized
    driver_names = propertynames(land_sim.integrator.p.drivers)
    @test driver_names == (:P_liq, :P_snow, :c_co2, :T, :P, :q, :SW_d, :LW_d, :cosθs, :frac_diff, :soc)
end
