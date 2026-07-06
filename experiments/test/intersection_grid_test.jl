#
# Intersection Grid Flux Exchange Tests
#
# This test verifies the intersection-grid flux exchange mechanism for better
# coastline representation when coupling atmosphere (ClimaCore cubed-sphere)
# with ocean/sea-ice (Oceananigans `TripolarGrid`).
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
import Statistics: std, mean

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

    # Create Oceananigans grid — small `TripolarGrid`. The 
    # displaced north poles sit at `north_poles_latitude = 55°N` (over
    # Greenland / Eurasia in the production run); the southern boundary is
    # a regular latitude circle at `southernmost_latitude = -80°`. 
    Nx, Ny, Nz = 36, 18, 1
    underlying_grid = OC.TripolarGrid(
        arch;
        size = (Nx, Ny, Nz),
        southernmost_latitude = -80,
        north_poles_latitude = 55,
        first_pole_longitude = 70,
        fold_topology = OC.RightCenterFolded,
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
        @test ig.n_oc == Nx * Ny  # 36 * 18 = 648 tripolar OC cells
        @test ig.n_intersections > 0
        @test length(ig.cc_indices) == ig.n_intersections
        @test length(ig.oc_indices) == ig.n_intersections
        @test length(ig.areas) == ig.n_intersections
        @test all(ig.areas .> 0)
        @test all(ig.cc_areas .> 0)
        # CR's `PaddedTreeWrapper` (used by `TripolarGrid` under the
        # `RightCenterFolded` topology) returns the folded half of the top
        # row as zero-area "shadow" cells: their indices still exist in
        # the flat `Nx*Ny` layout (because the live OC tracer arrays use
        # that layout) but they contribute no intersection geometry. We
        # require strict positivity on every other cell.
        n_shadow = count(==(0), ig.oc_areas)
        @test n_shadow == Nx ÷ 2  # 18 = top-row fold shadows
        @test all(ig.oc_areas[ig.oc_areas .> 0] .> 0)
        @test ig.n_nodes == 16 * ig.n_cc  # Nq² = 16 for n_quad_points = 4
        @test !isempty(ig.node_gather_polygon)
        @test length(ig.node_gather_polygon) == length(ig.node_gather_node)
        @test length(ig.node_gather_polygon) == length(ig.node_gather_weight)
    end

    @testset "Nodal polygon gather (sin-cos SEM field)" begin
        ig = remapping.intersection_grid
        CRExt = CMIPExt.get_ConservativeRegriddingCCExt()
        Nq = 4

        coords = CC.Fields.coordinate_field(boundary_space)
        test_field = CC.Fields.zeros(boundary_space)
        λ = deg2rad.(coords.long)
        φ = deg2rad.(coords.lat)
        test_field .= sin.(2.0 .* λ) .* cos.(φ)

        nodal = CRExt.se_field_to_vec(test_field)

        gathered = zeros(FT, ig.n_intersections)
        CMIPExt.gather_cc_nodal_to_intersection!(gathered, ig, nodal)

        # Non-trivial spatial structure (not a constant field).
        @test std(gathered) > FT(0.01)
        @test extrema(gathered) != extrema(nodal)

        # Independent COO matvec matches the gather kernel.
        reference = zeros(FT, ig.n_intersections)
        @inbounds for idx in eachindex(ig.node_gather_polygon)
            k = ig.node_gather_polygon[idx]
            n = ig.node_gather_node[idx]
            reference[k] += ig.node_gather_weight[idx] * nodal[n]
        end
        @test all(isapprox.(gathered, reference; rtol = 0, atol = 1e-12))

        # Analytic bounds of sin(2λ) cos(φ) on the sphere (polygon average is
        # not necessarily within the parent element's nodal min/max).
        @test all(-1.0 .- 1e-6 .<= gathered .<= 1.0 .+ 1e-6)

        # Differs from broadcasting the per-element nodal mean (old gather path).
        element_means = [
            mean(nodal[((e - 1) * Nq^2 + 1):(e * Nq^2)]) for e in 1:ig.n_cc
        ]
        element_gathered = zeros(FT, ig.n_intersections)
        CMIPExt.gather_cc_to_intersection!(element_gathered, ig, element_means)
        @test maximum(abs.(gathered .- element_gathered)) > FT(1e-6)
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

        # Round-trip identity holds exactly on every OC cell that
        # carries non-zero intersection area. Fold-shadow cells receive
        # no intersection contribution by construction, so we mask them
        # out here (they will read their live model value at runtime via
        # `extract_oc_surface_state!`, not the scatter result).
        live = ig.oc_areas .> 0
        @test maximum(abs.(oc_result[live] .- oc_values[live])) < 1e-10
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

        @test all(cc_result .> 0)

        @test all(minimum(cc_values) .<= cc_result .<= maximum(cc_values))
    end

    @testset "Flux calculation on intersection grid" begin
        ig = remapping.intersection_grid
        flux_state = remapping.intersection_flux_state

        # Nodal atmosphere with a smooth sin-cos pattern on the cubed sphere
        coords = CC.Fields.coordinate_field(boundary_space)
        CRExt = CMIPExt.get_ConservativeRegriddingCCExt()
        λ = deg2rad.(coords.long)
        φ = deg2rad.(coords.lat)
        trig = sin.(2.0 .* λ) .* cos.(φ)

        T_field = CC.Fields.zeros(boundary_space)
        T_field .= FT(280) .+ FT(5) .* trig
        u_field = CC.Fields.zeros(boundary_space)
        u_field .= FT(5) .+ FT(2) .* cos.(φ)

        n_nodes = ig.n_nodes
        cc_atmos_state = (
            T = CRExt.se_field_to_vec(T_field),
            q_tot = fill(FT(0.01), n_nodes),
            q_liq = fill(FT(0), n_nodes),
            q_ice = fill(FT(0), n_nodes),
            ρ = fill(FT(1.2), n_nodes),
            u = CRExt.se_field_to_vec(u_field),
            v = fill(FT(0), n_nodes),
            h = fill(FT(10), n_nodes),
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

        # Verify flux signs and magnitudes are physically reasonable.
        # T(λ, φ) = 280 + 5 sin(2λ) cos(φ) stays below the 300 K ocean.
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

        @testset "SEM nodal remapping via ConservativeRegridding" begin
            # Constant polygon flux → constant at every GLL node.
            flux_state_const = CMIPExt.IntersectionFluxState(FT, ig.n_intersections)
            const_flux = FT(42.0)
            fill!(flux_state_const.flux_sh, const_flux)
            fluxes_const = CMIPExt.intersection_fluxes_to_boundary_fields(
                boundary_space,
                remapping,
                flux_state_const,
            )
            nodal_const = CRExt.se_field_to_vec(fluxes_const.F_sh)
            @test all(isapprox.(nodal_const, const_flux; rtol = 0, atol = 1e-8))

            # Spatially varying fluxes → nodal structure within SE elements.
            fluxes = CMIPExt.intersection_fluxes_to_boundary_fields(
                boundary_space,
                remapping,
                flux_state,
            )
            p = parent(CC.Fields.field_values(fluxes.F_sh))
            intra_element_std = [std(vec(p[:, :, 1, e])) for e in 1:size(p, 4)]
            @test any(intra_element_std .> FT(0.01))

            # Differs from the legacy per-element constant broadcast.
            cc_F_sh = zeros(FT, ig.n_cc)
            CMIPExt.aggregate_fluxes_to_cc!(
                (; F_sh = cc_F_sh, F_lh = zeros(FT, ig.n_cc), F_τx = zeros(FT, ig.n_cc),
                   F_τy = zeros(FT, ig.n_cc), F_evap = zeros(FT, ig.n_cc)),
                flux_state,
                ig,
            )
            F_sh_legacy = CC.Fields.zeros(boundary_space)
            CMIPExt._element_values_to_se_field!(F_sh_legacy, cc_F_sh, boundary_space)
            @test maximum(abs.(fluxes.F_sh .- F_sh_legacy)) > FT(1e-6)
        end
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

        # With a `TripolarGrid` every CC element is tiled by intersection
        # polygons, so the warm-ocean / cool-atmosphere setup produces
        # strictly positive sensible and latent heat fluxes on every CC
        # element.
        @test all(cc_fluxes.F_sh .> 0)
        @test all(cc_fluxes.F_lh .> 0)

        # Test aggregation to OC grid
        oc_fluxes = (
            F_sh = zeros(FT, ig.n_oc),
            F_lh = zeros(FT, ig.n_oc),
            F_τx = zeros(FT, ig.n_oc),
            F_τy = zeros(FT, ig.n_oc),
            F_evap = zeros(FT, ig.n_oc),
        )
        CMIPExt.aggregate_fluxes_to_oc!(oc_fluxes, flux_state, ig)

        # On a tripolar grid every *live* underlying OC cell (i.e. every
        # cell with non-zero intersection area; fold-shadow cells on the
        # folded top row are explicitly excluded) is tiled by at least
        # one intersection polygon, so the warm-ocean / cool-atmosphere
        # setup produces strictly positive turbulent fluxes everywhere
        # the FV grid actually carries physics.
        live = ig.oc_areas .> 0
        @test all(oc_fluxes.F_sh[live] .> 0)
        @test all(oc_fluxes.F_lh[live] .> 0)
    end

    @testset "Integration with construct_remapper" begin
        # Verify all expected fields are present in remapping
        @test haskey(remapping, :intersection_grid)
        @test haskey(remapping, :intersection_flux_state)
        @test haskey(remapping, :cc_atmos_temp)
        @test haskey(remapping, :oc_surface_temp)

        # Verify temporary arrays have correct sizes
        ig = remapping.intersection_grid
        @test length(remapping.cc_atmos_temp.T) == ig.n_nodes
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

        # CC direction: `scatter_to_cc!` normalizes by the total CC element
        # area. With a `TripolarGrid` the OC mesh fully tiles every CC
        # element, so `cc_intersection_area[i] == ig.cc_areas[i]` for all
        # `i` and the area-weighted mean of a uniform field reproduces the
        # uniform value exactly. (For a lat-long grid truncated at ±80° we
        # would instead have `cc_result[i] = uniform_value * intersection_area[i] / cc_areas[i]`
        # in the polar band — that branch is no longer needed.)
        cc_result = zeros(FT, ig.n_cc)
        CMIPExt.scatter_to_cc!(cc_result, ig, intersection_values)

        cc_intersection_area = zeros(FT, ig.n_cc)
        for k in 1:ig.n_intersections
            cc_intersection_area[ig.cc_indices[k]] += ig.areas[k]
        end
        @test maximum(abs.(cc_intersection_area .- ig.cc_areas) ./ ig.cc_areas) < 1e-10
        @test maximum(abs.(cc_result .- uniform_value)) < 1e-10

        # The total flux on the CC grid must equal the total on the
        # intersection grid: `sum(cc_result * cc_areas) = uniform_value * 4πR²`.
        total_cc = sum(cc_result .* ig.cc_areas)
        @test abs(total_intersection - total_cc) / abs(total_intersection) < 1e-10
    end
end

println("\n✓ All intersection grid tests passed!")
