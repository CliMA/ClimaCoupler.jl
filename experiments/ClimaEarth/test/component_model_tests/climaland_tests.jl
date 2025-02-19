# import Test: @test, @testset

import ClimaCore as CC
import ClimaCoupler
include(joinpath("..", "TestHelper.jl"))
import .TestHelper
include(joinpath("..", "..", "components", "land", "climaland_integrated.jl"))

FT = Float64
atmos_boundary_space = TestHelper.create_space(FT)
land_sim = ClimaLandSimulation(
    FT;
    tspan = (FT(0), FT(99999)),
    dt = FT(1000),
    atmos_boundary_space,
    n_vertical_elements = 10,
    start_date = Dates.DateTime(2008),
    area_fraction = nothing,
);
Interfacer.step!(land_sim, FT(1));
