using Test
import ClimaCore as CC
import ClimaCoupler
import ClimaCoupler: FluxCalculator, Interfacer
import Dates

include(joinpath("..", "TestHelper.jl"))
import .TestHelper
include(joinpath("..", "..", "components", "land", "climaland_integrated.jl"))
include(joinpath("..", "..", "components", "atmosphere", "climaatmos.jl"))

FT = Float32

# @testset "ClimaLandSimulation constructor" begin
#     dt = Float64(450)
#     tspan = (Float64(0), 3.0dt)
#     start_date = Dates.DateTime(2008)
#     output_dir = pwd()
#     boundary_space = TestHelper.create_space(FT)
#     area_fraction = CC.Fields.ones(boundary_space)

#     # Construct simulation object
#     land_sim = ClimaLandSimulation(FT, dt, tspan, start_date, output_dir, boundary_space, area_fraction)

#     # Try taking a timestep
#     Interfacer.step!(land_sim, dt)

#     # Check that the simulation object is correctly initialized
#     @test Interfacer.name(land_sim) == "ClimaLandSimulation"
#     @test land_sim.area_fraction == area_fraction

#     # Check that the state is correctly initialized
#     state_names = propertynames(land_sim.integrator.u)
#     @test :canopy in state_names
#     @test :soil in state_names
#     @test :snow in state_names
#     @test :soilco2 in state_names

#     # Check that the cache is correctly initialized
#     cache_names = propertynames(land_sim.integrator.p)
#     @test :canopy in cache_names
#     @test :soil in cache_names
#     @test :snow in cache_names
#     @test :soilco2 in cache_names
#     @test :drivers in cache_names

#     # Check that the drivers are correctly initialized
#     driver_names = propertynames(land_sim.integrator.p.drivers)
#     @test driver_names == (:P_liq, :P_snow, :c_co2, :T, :P, :q, :SW_d, :LW_d, :cosθs, :frac_diff, :soc)
# end

@testset "ClimaLandSimulation flux calculations" begin
    dt = Float64(120)
    tspan = (Float64(0), 3.0dt)
    start_date = Dates.DateTime(2008)
    output_dir = pwd()

    # Construct atmos and land simulation objects
    atmos_config_file = joinpath(
        pkgdir(ClimaCoupler),
        "experiments",
        "ClimaEarth",
        "test",
        "component_model_tests",
        "climaatmos_coarse_short.yml",
    )
    atmos_config = CA.AtmosConfig(atmos_config_file; job_id = "atmos_land_flux_test")
    atmos_sim = ClimaAtmosSimulation(atmos_config)

    boundary_space = ClimaCore.Spaces.horizontal_space(atmos_sim.domain.face_space)
    area_fraction = ClimaCore.Fields.ones(boundary_space)
    land_sim = ClimaLandSimulation(FT, dt, tspan, start_date, output_dir, boundary_space, area_fraction)

    # Construct a coupler fields object
    coupler_fluxes = (;
        :F_turb_ρτxz => CC.Fields.zeros(boundary_space),
        :F_turb_ρτyz => CC.Fields.zeros(boundary_space),
        :F_turb_energy => CC.Fields.zeros(boundary_space),
        :F_turb_moisture => CC.Fields.zeros(boundary_space),
    )

    # Compute the surface fluxes
    FluxCalculator.compute_surface_fluxes!(coupler_fluxes, land_sim, atmos_sim, boundary_space, nothing, nothing)

    # Check that the fluxes have been changed
    # TODO only energy flux gets changed
    zero_field = CC.Fields.zeros(boundary_space)
    @test coupler_fluxes.F_turb_ρτxz == zero_field
    @test coupler_fluxes.F_turb_ρτyz == zero_field
    @test coupler_fluxes.F_turb_energy != zero_field
    @test coupler_fluxes.F_turb_moisture == zero_field
end
