import Test: @test, @testset
import Statistics: cor, mean
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC
import Oceananigans as OC
import ClimaOcean as CO
import Adapt
import KernelAbstractions
import ConservativeRegridding
import ClimaSeaIce
import ClimaCoupler
import ClimaCoupler: Interfacer, Input

CMIPExt = Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt)
@assert !isnothing(CMIPExt)

@testset "land_fraction_from_ocean" begin
    FT = Float64
    context = ClimaComms.context()
    ClimaComms.init(context)
    arch = OC.CPU()

    boundary_space = CC.CommonSpaces.CubedSphereSpace(
        FT;
        radius = FT(6.371e6),
        n_quad_points = 4,
        h_elem = 4,
        context,
    )

    underlying_grid = OC.TripolarGrid(
        arch;
        size = (36, 18, 1),
        southernmost_latitude = -80,
        north_poles_latitude = 55,
        first_pole_longitude = 70,
        fold_topology = OC.RightCenterFolded,
        z = (-100.0, 0.0),
        halo = (4, 4, 4),
    )
    bottom_height = CO.regrid_bathymetry(
        underlying_grid;
        minimum_depth = 30,
        interpolation_passes = 5,
        major_basins = 1,
    )
    grid = OC.ImmersedBoundaryGrid(
        underlying_grid,
        OC.GridFittedBottom(bottom_height);
        active_cells_map = false,
    )

    remapping = CMIPExt.construct_remapper(grid, boundary_space)
    area_fraction = ones(boundary_space)
    component_area_fraction = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    surface_mask = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    CMIPExt.fill_ocean_wet_mask!(surface_mask, grid)
    ice_concentration = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    ocean = (model = (grid = grid,),)
    ocean_sim = CMIPExt.OceananigansSimulation(
        ocean,
        area_fraction,
        component_area_fraction,
        surface_mask,
        (;),
        remapping,
        ice_concentration,
        1800.0,
    )

    land_fraction = CMIPExt.land_fraction_from_ocean(
        ocean_sim,
        boundary_space;
        binary_area_fraction = true,
    )

    land_arr = CC.Fields.field2array(land_fraction)
    @test minimum(land_arr) >= 0 && maximum(land_arr) <= 1
    @test 0 < count(land_arr .> 0.5) < length(land_arr)
    @test 0 < count(land_arr .< 0.5) < length(land_arr)
    # Earth-like global land coverage (not ~100% land from a misapplied binary threshold)
    @test 0.15 < mean(land_arr) < 0.55

    etopo_land = Input.get_land_fraction(
        boundary_space,
        context;
        land_fraction_source = "etopo",
        binary_area_fraction = true,
    )
    etopo_arr = CC.Fields.field2array(etopo_land)
    mismatch = sum(abs.(land_arr .- etopo_arr))
    @test mismatch > 0
    inverted_mismatch = sum(abs.(land_arr .- (1 .- etopo_arr)))
    @test mismatch <= inverted_mismatch
    # Same sign as ETOPO (both increase over continents)
    @test cor(land_arr, etopo_arr) > 0.2

    land_fraction_binary_flag = CMIPExt.land_fraction_from_ocean(
        ocean_sim,
        boundary_space;
        binary_area_fraction = false,
    )
    @test land_fraction_binary_flag ≈ land_fraction

    ocean_sim.area_fraction .= FT(1) .- land_fraction
    CMIPExt.update_component_area_fraction!(ocean_sim)
    comp_arr = OC.interior(ocean_sim.component_area_fraction, :, :, 1)
    wet_arr = OC.interior(ocean_sim.surface_mask, :, :, 1)
    # ice concentration is zero at init, so component fraction matches the wet mask
    @test comp_arr ≈ wet_arr
    @test all(isfinite, comp_arr)
    @test extrema(comp_arr) == (0, 1)

    @test hasproperty(ocean_sim.remapping, :fv_wet_mask_1d)
    @test extrema(ocean_sim.remapping.fv_wet_mask_1d) == (0, 1)
end
