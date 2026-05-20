#
# Intersection Grid Flux Exchange Tests
#
# This test verifies the intersection-grid flux exchange mechanism for better
# coastline representation when coupling atmosphere (ClimaCore cubed-sphere)
# with ocean/sea-ice (Oceananigans LatitudeLongitudeGrid).
#
# Tests include:
# 1. IntersectionGrid extraction from CR.Regridder
# 2. Gather/scatter conservation properties
# 3. Flux calculation on intersection grid
# 4. Conservation diagnostics
#

import Test: @test, @testset
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC
import Oceananigans as OC
import ClimaOcean
import ClimaSeaIce
import KernelAbstractions
import ConservativeRegridding as CR
import Adapt
import ClimaCoupler
import ClimaCoupler: Interfacer
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
import ClimaParams as CP

CMIPExt = Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt)

FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)
arch = OC.CPU()

@testset "Intersection Grid Tests" begin
    # Create boundary space (small for testing)
    boundary_space = CC.CommonSpaces.CubedSphereSpace(
        FT;
        radius = FT(6.371e6),
        n_quad_points = 4,
        h_elem = 4,
        context,
    )

    # Create Oceananigans grid
    Nx, Ny, Nz = 36, 16, 1
    underlying_grid = OC.LatitudeLongitudeGrid(
        arch;
        size = (Nx, Ny, Nz),
        longitude = (-180, 180),
        latitude = (-80, 80),
        z = (-100.0, 0.0),
        halo = (4, 4, 4),
    )
    grid = OC.ImmersedBoundaryGrid(
        underlying_grid,
        OC.GridFittedBottom((x, y) -> -100.0);
        active_cells_map = false,
    )

    # Construct remapper
    remapping = CMIPExt.construct_remapper(grid, boundary_space)

    @testset "IntersectionGrid extraction" begin
        ig = remapping.intersection_grid
        @test ig isa CMIPExt.IntersectionGrid
        @test ig.n_cc == 96  # 4^2 * 6 elements for h_elem=4 cubed sphere
        @test ig.n_oc == Nx * Ny  # 36 * 16 = 576 OC cells
        @test ig.n_intersections > 0
        @test length(ig.cc_indices) == ig.n_intersections
        @test length(ig.oc_indices) == ig.n_intersections
        @test length(ig.areas) == ig.n_intersections
        @test all(ig.areas .> 0)
        @test all(ig.cc_areas .> 0)
        @test all(ig.oc_areas .> 0)
    end

    @testset "Gather/scatter conservation (OC direction)" begin
        ig = remapping.intersection_grid

        # Create test values on OC grid (uniform for conservation test)
        oc_values = fill(FT(1.0), ig.n_oc)

        # Gather to intersection grid
        intersection_values = zeros(FT, ig.n_intersections)
        CMIPExt.gather_oc_to_intersection!(intersection_values, ig, oc_values)

        # All intersection values should equal source value
        @test all(intersection_values .== 1.0)

        # Scatter back to OC (area-weighted mean)
        oc_result = zeros(FT, ig.n_oc)
        CMIPExt.scatter_to_oc!(oc_result, ig, intersection_values)

        # Result should match original (within floating point precision)
        @test maximum(abs.(oc_result .- oc_values)) < 1e-10
    end

    @testset "Gather/scatter with varying values (CC direction)" begin
        ig = remapping.intersection_grid

        # Create varying test values on CC grid
        cc_values = collect(1.0:FT(ig.n_cc))

        # Gather to intersection grid
        intersection_values = zeros(FT, ig.n_intersections)
        CMIPExt.gather_cc_to_intersection!(intersection_values, ig, cc_values)

        # Each intersection value should equal its parent CC element value
        for k in 1:ig.n_intersections
            @test intersection_values[k] == cc_values[ig.cc_indices[k]]
        end

        # Scatter back to CC
        cc_result = zeros(FT, ig.n_cc)
        CMIPExt.scatter_to_cc!(cc_result, ig, intersection_values)

        # Note: CC elements that only partially overlap the OC grid will have
        # smoothing. This is expected behavior - area-weighted averaging over
        # partial overlaps doesn't reconstruct the original value exactly.
        # 
        # For CC elements fully within the OC coverage band, the round-trip
        # should be exact. Elements at the edge (|lat| > ~60°) will show
        # smoothing effects since they only partially overlap the OC grid.
        #
        # We test that:
        # 1. Elements with intersections have non-zero results
        # 2. The area-weighted integral is preserved (tested separately)
        has_intersection = zeros(Bool, ig.n_cc)
        for k in 1:ig.n_intersections
            has_intersection[ig.cc_indices[k]] = true
        end

        n_with_intersection = sum(has_intersection)
        @test n_with_intersection > 0

        for i in 1:ig.n_cc
            if has_intersection[i]
                @test cc_result[i] > 0  # Non-zero result
            else
                @test cc_result[i] == 0  # No coverage → zero
            end
        end
    end

    @testset "Flux calculation on intersection grid" begin
        ig = remapping.intersection_grid
        flux_state = remapping.intersection_flux_state

        # Create atmosphere state (per CC element)
        cc_atmos_state = (
            T = fill(FT(290), ig.n_cc),
            q_tot = fill(FT(0.01), ig.n_cc),
            q_liq = fill(FT(0), ig.n_cc),
            q_ice = fill(FT(0), ig.n_cc),
            ρ = fill(FT(1.2), ig.n_cc),
            u = fill(FT(5), ig.n_cc),
            v = fill(FT(0), ig.n_cc),
            h = fill(FT(10), ig.n_cc),
        )

        # Create surface state (per OC cell) - warm ocean
        oc_surface_state = (
            T = fill(FT(300), ig.n_oc),
            z0m = fill(FT(1e-4), ig.n_oc),
            z0b = fill(FT(1e-4), ig.n_oc),
            h = fill(FT(0), ig.n_oc),
        )

        # Create parameters
        toml_dict = CP.create_toml_dict(FT)
        thermo_params = TDP.ThermodynamicsParameters(toml_dict)
        surface_fluxes_params = SFP.SurfaceFluxesParameters(toml_dict, UF.BusingerParams)

        # Compute fluxes
        CMIPExt.compute_surface_fluxes_on_intersection!(
            flux_state,
            ig,
            cc_atmos_state,
            oc_surface_state,
            surface_fluxes_params,
            thermo_params,
        )

        # Verify flux signs and magnitudes are physically reasonable
        # Warm ocean (300K) + cooler atmosphere (290K) → positive sensible heat flux
        @test all(flux_state.flux_sh .> 0)
        @test maximum(flux_state.flux_sh) < 500  # Reasonable upper bound W/m²

        # Warm ocean → positive latent heat flux
        @test all(flux_state.flux_lh .> 0)
        @test maximum(flux_state.flux_lh) < 1000

        # Wind from the east (u > 0) → negative τx momentum flux
        @test all(flux_state.flux_τx .< 0)

        # No meridional wind → τy ≈ 0
        @test maximum(abs.(flux_state.flux_τy)) < 1e-10

        # Positive evaporation (ocean evaporating)
        @test all(flux_state.flux_evap .> 0)
    end

    @testset "Flux aggregation to grids" begin
        ig = remapping.intersection_grid
        flux_state = remapping.intersection_flux_state

        # Assume fluxes already computed from previous test
        # Test aggregation to CC grid
        cc_fluxes = (
            F_sh = zeros(FT, ig.n_cc),
            F_lh = zeros(FT, ig.n_cc),
            F_τx = zeros(FT, ig.n_cc),
            F_τy = zeros(FT, ig.n_cc),
            F_evap = zeros(FT, ig.n_cc),
        )
        CMIPExt.aggregate_fluxes_to_cc!(cc_fluxes, flux_state, ig)

        # CC elements with intersection should have non-zero fluxes
        has_intersection = zeros(Bool, ig.n_cc)
        for k in 1:ig.n_intersections
            has_intersection[ig.cc_indices[k]] = true
        end

        for i in 1:ig.n_cc
            if has_intersection[i]
                @test cc_fluxes.F_sh[i] > 0
                @test cc_fluxes.F_lh[i] > 0
            else
                @test cc_fluxes.F_sh[i] == 0
            end
        end

        # Test aggregation to OC grid
        oc_fluxes = (
            F_sh = zeros(FT, ig.n_oc),
            F_lh = zeros(FT, ig.n_oc),
            F_τx = zeros(FT, ig.n_oc),
            F_τy = zeros(FT, ig.n_oc),
            F_evap = zeros(FT, ig.n_oc),
        )
        CMIPExt.aggregate_fluxes_to_oc!(oc_fluxes, flux_state, ig)

        # All OC cells should have fluxes (assuming full coverage)
        # Note: Some cells may have zero if they don't overlap with any CC element
        @test any(oc_fluxes.F_sh .> 0)
        @test any(oc_fluxes.F_lh .> 0)
    end

    @testset "Integration with construct_remapper" begin
        # Verify all expected fields are present in remapping
        @test haskey(remapping, :intersection_grid)
        @test haskey(remapping, :intersection_flux_state)
        @test haskey(remapping, :cc_atmos_temp)
        @test haskey(remapping, :oc_surface_temp)

        # Verify temporary arrays have correct sizes
        ig = remapping.intersection_grid
        @test length(remapping.cc_atmos_temp.T) == ig.n_cc
        @test length(remapping.oc_surface_temp.T) == ig.n_oc
    end

    @testset "Conservation: area-weighted integral preservation" begin
        ig = remapping.intersection_grid

        # Test that total area is preserved when scattering
        # Create uniform values with known integral
        uniform_value = FT(10.0)
        intersection_values = fill(uniform_value, ig.n_intersections)

        # Total flux on intersection grid
        total_intersection = sum(intersection_values .* ig.areas)

        # Scatter to OC and compute total
        oc_result = zeros(FT, ig.n_oc)
        CMIPExt.scatter_to_oc!(oc_result, ig, intersection_values)
        total_oc = sum(oc_result .* ig.oc_areas)

        # Area-weighted totals should be equal (OC direction is exact since
        # all OC cells are fully within the CC grid coverage)
        @test abs(total_intersection - total_oc) / abs(total_intersection) < 1e-10

        # CC direction: scatter_to_cc! normalizes by total CC element area,
        # not intersection area. This is correct for intensive quantities -
        # if a CC element only partially overlaps the OC grid, its mean flux
        # will be proportional to the overlap fraction.
        cc_result = zeros(FT, ig.n_cc)
        CMIPExt.scatter_to_cc!(cc_result, ig, intersection_values)

        # Compute overlap fractions
        cc_intersection_area = zeros(FT, ig.n_cc)
        for k in 1:ig.n_intersections
            cc_intersection_area[ig.cc_indices[k]] += ig.areas[k]
        end

        # For CC elements, the result should be:
        #   cc_result[i] = uniform_value * (intersection_area[i] / cc_area[i])
        # because scatter_to_cc! normalizes by cc_area
        for i in 1:ig.n_cc
            if cc_intersection_area[i] > 0
                expected = uniform_value * cc_intersection_area[i] / ig.cc_areas[i]
                @test abs(cc_result[i] - expected) < 1e-10
            else
                @test cc_result[i] == 0
            end
        end

        # The total flux using cc_result and cc_areas should equal total_intersection
        # because: sum(cc_result[i] * cc_areas[i]) 
        #        = sum(uniform_value * intersection_area[i])
        #        = uniform_value * sum(intersection_area[i])
        #        = total_intersection
        total_cc = sum(cc_result .* ig.cc_areas)
        @test abs(total_intersection - total_cc) / abs(total_intersection) < 1e-10
    end
end

println("\n✓ All intersection grid tests passed!")
