#=
# Exchange (intersection) grid tests

Unit tests for the exchange grid used to couple the ClimaCore cubed-sphere
spectral-element boundary space with an Oceananigans `TripolarGrid`:
geometry/weight construction (`build_exchange_grid`), gather/scatter
conservation, area fractions, and per-polygon fluxes.

Run with the CMIP project environment, which provides Oceananigans and
ConservativeRegridding:

    julia --project=experiments/CMIP experiments/test/intersection_grid_test.jl
=#

import Test: @test, @testset
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC
import Oceananigans as OC
import ClimaOcean
import ClimaSeaIce
import ClimaAtmos
import KernelAbstractions
import ConservativeRegridding as CR
import Adapt
import ClimaCoupler
import Statistics: mean

CMIPExt = Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt)
@assert !isnothing(CMIPExt)

const FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)

radius = FT(6.371e6)
boundary_space = CC.CommonSpaces.CubedSphereSpace(
    FT;
    radius,
    n_quad_points = 4,
    h_elem = 4,
    context,
)
Nq = 4
n_elem = 6 * 4^2

# Small `TripolarGrid` matching the production layout: displaced north poles
# at 55°N, regular southern boundary at 80°S, folded top row.
Nx, Ny, Nz = 36, 18, 1
arch = OC.CPU()
underlying_tripolar = OC.TripolarGrid(
    arch;
    size = (Nx, Ny, Nz),
    southernmost_latitude = -80,
    north_poles_latitude = 55,
    first_pole_longitude = 70,
    fold_topology = OC.RightCenterFolded,
    z = (-100.0, 0.0),
    halo = (4, 4, 4),
)

# Fully wet grid: bottom below the deepest z everywhere.
wet_grid = OC.ImmersedBoundaryGrid(
    underlying_tripolar,
    OC.GridFittedBottom((x, y) -> -200.0);
    active_cells_map = false,
)

# Half-dry grid: an idealized land hemisphere. TripolarGrid longitudes are
# not restricted to [-180, 180), so wrap before comparing.
coastal_grid = OC.ImmersedBoundaryGrid(
    underlying_tripolar,
    OC.GridFittedBottom((x, y) -> mod(x, 360) < 180 ? 100.0 : -200.0);
    active_cells_map = false,
)

# Ground-truth wet area of an immersed grid: the sum of the cell areas of
# non-immersed cells, excluding the duplicated (shadow) half of the tripolar
# fold row, which Oceananigans stores but ConservativeRegridding treats as
# degenerate.
function expected_wet_area(grid, shadow_cells)
    grid_cpu = OC.on_architecture(OC.CPU(), grid)
    total = 0.0
    for j in 1:Ny, i in 1:Nx
        c = (j - 1) * Nx + i
        c in shadow_cells && continue
        OC.ImmersedBoundaries.immersed_cell(i, j, Nz, grid_cpu) && continue
        total += OC.Operators.Azᶜᶜᶜ(i, j, Nz, grid_cpu.underlying_grid)
    end
    return total
end

@testset "Exchange grid" begin
    eg = CMIPExt.build_exchange_grid(boundary_space, wet_grid)

    @testset "geometry (fully wet)" begin
        @test eg isa CMIPExt.ExchangeGrid{FT}
        @test eg.n_elem == n_elem
        @test eg.n_oc == Nx * Ny
        @test eg.n_nodes == Nq^2 * n_elem
        @test eg.n_poly > 0
        @test length(eg.area) == eg.n_poly
        @test length(eg.elem_of_poly) == eg.n_poly
        @test length(eg.oc_of_poly) == eg.n_poly
        @test all(eg.area .>= 0)

        # The polygons tile the ocean grid: the grid ends at 80°S, so the
        # total is slightly less than the sphere area but matches the sum of
        # the wet, non-shadow FV cell areas from the Oceananigans metrics.
        n_shadow = count(==(0), eg.oc_wet_area)
        @test n_shadow == Nx ÷ 2
        shadow_cells = Set(findall(==(0), eg.oc_wet_area))
        @test all(c -> (c - 1) ÷ Nx + 1 == Ny, shadow_cells)
        @test sum(eg.area) < 4π * radius^2
        @test sum(eg.area) ≈ expected_wet_area(wet_grid, shadow_cells) rtol = 1e-6

        # Every wet, non-shadow FV cell is exactly tiled by its polygons.
        grid_cpu = OC.on_architecture(OC.CPU(), wet_grid)
        for c in 1:(eg.n_oc)
            c in shadow_cells && continue
            i, j = mod1(c, Nx), (c - 1) ÷ Nx + 1
            Az = OC.Operators.Azᶜᶜᶜ(i, j, Nz, grid_cpu.underlying_grid)
            @test eg.oc_wet_area[c] ≈ Az rtol = 1e-6
        end

        # Per-element polygon areas match the element areas (the integral of
        # the nodal quadrature weights Jw over each element) wherever the
        # element lies fully inside the ocean grid (north of 80°S).
        CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
        Jw = CRExt.se_node_weights(boundary_space)
        coords = CC.Fields.coordinate_field(boundary_space)
        node_lat = CRExt.flat_nodal_data(CC.Fields.field_values(coords.lat))
        elem_area_from_Jw = zeros(FT, n_elem)
        elem_min_lat = fill(FT(90), n_elem)
        for n in 1:(eg.n_nodes)
            e = (n - 1) ÷ Nq^2 + 1
            elem_area_from_Jw[e] += Jw[n]
            elem_min_lat[e] = min(elem_min_lat[e], node_lat[n])
        end
        elem_area_from_polys = zeros(FT, n_elem)
        for k in 1:(eg.n_poly)
            elem_area_from_polys[eg.elem_of_poly[k]] += eg.area[k]
        end
        for e in 1:n_elem
            if elem_min_lat[e] > -75
                @test elem_area_from_polys[e] ≈ elem_area_from_Jw[e] rtol = 1e-4
            end
            # No element is over-covered by polygons.
            @test elem_area_from_polys[e] <= elem_area_from_Jw[e] * (1 + 1e-4)
        end

        # Per-cell wet areas are consistent with the polygon areas.
        @test sum(eg.oc_wet_area) ≈ sum(eg.area) rtol = 1e-14

        # CSR structure: every polygon couples to at most Nq² nodes of
        # exactly one element; the OC-major CSR lists each polygon once.
        for k in 1:(eg.n_poly)
            rng = eg.gpoly_ptr[k]:(eg.gpoly_ptr[k + 1] - 1)
            @test 0 < length(rng) <= Nq^2
            elems = unique((eg.gnode[rng] .- 1) .÷ Nq^2 .+ 1)
            @test elems == [eg.elem_of_poly[k]]
        end
        @test sort(eg.soc_poly) == 1:(eg.n_poly)
        for c in 1:(eg.n_oc)
            rng = eg.soc_ptr[c]:(eg.soc_ptr[c + 1] - 1)
            @test all(k -> eg.oc_of_poly[k] == c, eg.soc_poly[rng])
        end
    end

    @testset "partition of unity (fully wet)" begin
        # Gather rows sum to 1 exactly (they are normalized by Σₙ B_kn).
        row_sums = zeros(FT, eg.n_poly)
        for k in 1:(eg.n_poly)
            for p in eg.gpoly_ptr[k]:(eg.gpoly_ptr[k + 1] - 1)
                row_sums[k] += eg.gweight[p]
            end
        end
        @test all(s -> isapprox(s, 1, atol = 1e-12), row_sums)

        # On a fully wet grid every polygon is retained, so the wet and
        # geometric coverages coincide identically and their ratio — the wet
        # fraction — is exactly 1 wherever the ocean grid reaches. The raw
        # coverages themselves carry quadrature and quad-approximation error
        # (up to ~15% at fold/pole nodes on this deliberately coarse
        # fixture), which is exactly why the ratio is the meaningful
        # quantity.
        @test eg.node_cov == eg.node_cov_total
        CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
        coords = CC.Fields.coordinate_field(boundary_space)
        node_lat = CRExt.flat_nodal_data(CC.Fields.field_values(coords.lat))
        for n in 1:(eg.n_nodes)
            if node_lat[n] > -75
                # Interior nodes: geometric coverage is ≈ 1 to within the
                # coarse-fixture quad-approximation error.
                @test eg.node_cov_total[n] ≈ 1 rtol = 0.6
            end
        end
        # Nodes south of the ocean grid boundary are only partially covered.
        @test any(n -> node_lat[n] < -80 && eg.node_cov_total[n] < 0.9, 1:(eg.n_nodes))

        # node_cov is consistent with the scatter weights.
        cov_from_sweight = zeros(FT, eg.n_nodes)
        for n in 1:(eg.n_nodes)
            for p in eg.snode_ptr[n]:(eg.snode_ptr[n + 1] - 1)
                cov_from_sweight[n] += eg.sweight[p]
            end
        end
        @test cov_from_sweight ≈ eg.node_cov rtol = 1e-14
    end

    @testset "immersed mask (idealized land hemisphere)" begin
        eg_c = CMIPExt.build_exchange_grid(boundary_space, coastal_grid)

        # The retained polygon area matches the wet cell area exactly, and is
        # roughly half the sphere.
        shadow_cells = Set(findall(==(0), eg.oc_wet_area))
        @test sum(eg_c.area) ≈ expected_wet_area(coastal_grid, shadow_cells) rtol = 1e-6
        @test 0.4 * 4π * radius^2 < sum(eg_c.area) < 0.6 * 4π * radius^2

        # Dry cells have no polygons and no wet area.
        grid_cpu = OC.on_architecture(OC.CPU(), coastal_grid)
        for c in 1:(eg_c.n_oc)
            i, j = mod1(c, Nx), (c - 1) ÷ Nx + 1
            if OC.ImmersedBoundaries.immersed_cell(i, j, Nz, grid_cpu)
                @test eg_c.oc_wet_area[c] == 0
            end
        end

        # The wet fraction is the ratio of wet to geometric coverage:
        # exactly 1 at nodes all of whose polygons are wet, exactly 0 at
        # interior dry nodes, intermediate (with Gibbs overshoot, clamped by
        # consumers) at the coastline.
        wet_frac = [
            t > 0 ? w / t : zero(w) for
            (w, t) in zip(eg_c.node_cov, eg_c.node_cov_total)
        ]
        @test count(==(1), wet_frac) > 0.2 * length(wet_frac)
        @test count(==(0), wet_frac) > 0.2 * length(wet_frac)
        @test any(f -> 0.05 < f < 0.95, wet_frac)
        @test all(f -> -1 <= f <= 2, wet_frac)
        clamped = clamp.(wet_frac, 0, 1)
        @test any(f -> 0.05 < f < 0.95, clamped)
    end
end

import Random

@testset "Gather/scatter operations" begin
    eg = CMIPExt.build_exchange_grid(boundary_space, coastal_grid)
    CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
    Jw = CRExt.se_node_weights(boundary_space)
    rng = Random.Xoshiro(42)

    poly_scratch = zeros(FT, eg.n_poly)
    nodal_scratch = zeros(FT, eg.n_nodes)
    cell_scratch = zeros(FT, eg.n_oc)

    @testset "constant preservation" begin
        # Gather: rows sum to 1, so a constant nodal field is exact.
        CMIPExt.gather_nodes_to_polys!(poly_scratch, eg, ones(FT, eg.n_nodes))
        @test all(v -> isapprox(v, 1, atol = 1e-12), poly_scratch)

        # Raw scatter of a constant reproduces the nodal wet coverage.
        CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, ones(FT, eg.n_poly))
        @test nodal_scratch ≈ eg.node_cov rtol = 1e-12

        # Normalized scatter of a constant is exactly the constant on
        # covered nodes and 0 elsewhere.
        CMIPExt.scatter_polys_to_nodes_normalized!(
            nodal_scratch,
            eg,
            ones(FT, eg.n_poly),
            FT(1e-3),
        )
        for n in 1:(eg.n_nodes)
            if eg.node_cov[n] > 1e-3
                @test nodal_scratch[n] ≈ 1 rtol = 1e-12
            else
                @test nodal_scratch[n] == 0
            end
        end

        # Cell scatter of a constant is exact on wet cells, 0 on dry/shadow
        # cells; the fold mirror then fills shadow partners.
        CMIPExt.scatter_polys_to_cells!(cell_scratch, eg, ones(FT, eg.n_poly))
        for c in 1:(eg.n_oc)
            if eg.oc_wet_area[c] > 0
                @test cell_scratch[c] ≈ 1 rtol = 1e-12
            else
                @test cell_scratch[c] == 0
            end
        end
        shadow_cells = [c for c in ((Ny - 1) * Nx + 1):(Ny * Nx) if eg.oc_wet_area[c] == 0]
        CMIPExt.mirror_fold_partners!(cell_scratch, coastal_grid)
        for c in shadow_cells
            j = Ny
            i = c - (Ny - 1) * Nx
            partner = (Ny - 1) * Nx + (Nx + 1 - i)
            @test cell_scratch[c] == cell_scratch[partner]
        end
    end

    @testset "cell gather is direct indexing" begin
        cell_values = rand(rng, FT, eg.n_oc)
        CMIPExt.gather_cells_to_polys!(poly_scratch, eg, cell_values)
        for k in 1:(eg.n_poly)
            @test poly_scratch[k] == cell_values[eg.oc_of_poly[k]]
        end
    end

    @testset "conservation" begin
        F = rand(rng, FT, eg.n_poly)

        # FV side: the scatter preserves the area integral exactly.
        CMIPExt.scatter_polys_to_cells!(cell_scratch, eg, F)
        @test sum(eg.area .* F) ≈ sum(eg.oc_wet_area .* cell_scratch) rtol = 1e-12

        # SE side: the L2 scatter preserves the quadrature-area integral
        # exactly (before DSS).
        CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, F)
        @test sum(Jw .* nodal_scratch) ≈ sum(eg.b_area .* F) rtol = 1e-12

        # And the quadrature areas agree with the geometric areas.
        @test sum(eg.b_area) ≈ sum(eg.area) rtol = 1e-3
    end

    @testset "allocation-free application (CPU)" begin
        nodal = rand(rng, FT, eg.n_nodes)
        polyv = rand(rng, FT, eg.n_poly)
        CMIPExt.gather_nodes_to_polys!(poly_scratch, eg, nodal)
        CMIPExt.gather_cells_to_polys!(poly_scratch, eg, cell_scratch)
        CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, polyv)
        CMIPExt.scatter_polys_to_nodes_normalized!(nodal_scratch, eg, polyv, FT(1e-3))
        CMIPExt.scatter_polys_to_cells!(cell_scratch, eg, polyv)
        CMIPExt.mirror_fold_partners!(cell_scratch, coastal_grid)
        @test @allocated(CMIPExt.gather_nodes_to_polys!(poly_scratch, eg, nodal)) == 0
        @test @allocated(CMIPExt.gather_cells_to_polys!(poly_scratch, eg, cell_scratch)) ==
              0
        @test @allocated(CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, polyv)) == 0
        @test @allocated(
            CMIPExt.scatter_polys_to_nodes_normalized!(
                nodal_scratch,
                eg,
                polyv,
                FT(1e-3),
            )
        ) == 0
        @test @allocated(CMIPExt.scatter_polys_to_cells!(cell_scratch, eg, polyv)) == 0
        @test @allocated(CMIPExt.mirror_fold_partners!(cell_scratch, coastal_grid)) == 0
    end

    if ClimaComms.device() isa ClimaComms.CUDADevice
        @testset "CPU == GPU" begin
            arch_gpu = OC.GPU()
            eg_d = CMIPExt.on_device(arch_gpu, eg)
            to_d(x) = CMIPExt.on_device(arch_gpu, x)

            nodal = rand(rng, FT, eg.n_nodes)
            polyv = rand(rng, FT, eg.n_poly)
            cellv = rand(rng, FT, eg.n_oc)
            rt = sqrt(eps(FT))

            CMIPExt.gather_nodes_to_polys!(poly_scratch, eg, nodal)
            @test Array(CMIPExt.gather_nodes_to_polys!(to_d(poly_scratch), eg_d, to_d(nodal))) ≈
                  poly_scratch rtol = rt

            CMIPExt.gather_cells_to_polys!(poly_scratch, eg, cellv)
            @test Array(CMIPExt.gather_cells_to_polys!(to_d(poly_scratch), eg_d, to_d(cellv))) ==
                  poly_scratch

            CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, polyv)
            @test Array(
                CMIPExt.scatter_polys_to_nodes!(to_d(nodal_scratch), eg_d, to_d(polyv)),
            ) ≈ nodal_scratch rtol = rt

            CMIPExt.scatter_polys_to_nodes_normalized!(
                nodal_scratch,
                eg,
                polyv,
                FT(1e-3),
            )
            @test Array(
                CMIPExt.scatter_polys_to_nodes_normalized!(
                    to_d(nodal_scratch),
                    eg_d,
                    to_d(polyv),
                    FT(1e-3),
                ),
            ) ≈ nodal_scratch rtol = rt

            CMIPExt.scatter_polys_to_cells!(cell_scratch, eg, polyv)
            cell_d = CMIPExt.scatter_polys_to_cells!(to_d(cell_scratch), eg_d, to_d(polyv))
            @test Array(cell_d) ≈ cell_scratch rtol = rt

            CMIPExt.mirror_fold_partners!(cell_scratch, coastal_grid)
            gpu_coastal = OC.on_architecture(arch_gpu, coastal_grid)
            @test Array(CMIPExt.mirror_fold_partners!(cell_d, gpu_coastal)) ≈ cell_scratch rtol =
                rt
        end
    end
end

@testset "Wet-ocean fraction field" begin
    eg_wet = CMIPExt.build_exchange_grid(boundary_space, wet_grid)
    eg_c = CMIPExt.build_exchange_grid(boundary_space, coastal_grid)

    CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
    coords = CC.Fields.coordinate_field(boundary_space)
    node_lat = CRExt.flat_nodal_data(CC.Fields.field_values(coords.lat))
    node_lon = CRExt.flat_nodal_data(CC.Fields.field_values(coords.long))

    # `topography_damping_factor = 1` means zero smoothing iterations
    # (maxiter = round(log(1)/0.05) = 0): the raw, DSS'd, clamped fraction.
    f_raw = CMIPExt.wet_ocean_fraction_field(
        boundary_space,
        eg_c;
        topography_damping_factor = 1,
    )
    f_smooth = CMIPExt.wet_ocean_fraction_field(
        boundary_space,
        eg_c;
        topography_damping_factor = 5,
    )
    raw = CRExt.se_field_to_vec(f_raw)
    smooth = CRExt.se_field_to_vec(f_smooth)

    # Bounds and non-triviality.
    @test all(v -> 0 <= v <= 1, raw)
    @test all(v -> 0 <= v <= 1, smooth)
    @test maximum(abs.(smooth .- raw)) > 0

    # The diffusion is (weak-form) conservative; with the clamp the global
    # integral moves only slightly. `sum` of a ClimaCore field is the
    # area-weighted integral.
    @test sum(f_smooth) ≈ sum(f_raw) rtol = 2e-2

    # Interior values: fully ocean far from the idealized coastlines at
    # wrapped longitudes 0 and 180, fully land in the middle of the land
    # hemisphere. (Land: mod(lon, 360) ∈ [0, 180).)
    for n in 1:(eg_c.n_nodes)
        lonw = mod(node_lon[n], 360)
        abs(node_lat[n]) > 45 && continue
        if 60 < lonw < 120 # deep inside land
            @test smooth[n] < 0.1
        elseif 240 < lonw < 300 # deep inside ocean
            @test smooth[n] > 0.9
        end
    end

    # Fully wet grid: fraction ≈ 1 everywhere the ocean grid reaches, and
    # the smoothing must not degrade it.
    f_wet = CMIPExt.wet_ocean_fraction_field(boundary_space, eg_wet)
    wet_vals = CRExt.se_field_to_vec(f_wet)
    for n in 1:(eg_wet.n_nodes)
        if node_lat[n] > -60
            @test wet_vals[n] > 0.95
        end
    end
end

import SurfaceFluxes as SF
import Thermodynamics.Parameters as TDP
import ClimaParams as CP
import ClimaCoupler: FluxCalculator, Utilities

@testset "Per-polygon ocean fluxes" begin
    eg = CMIPExt.build_exchange_grid(boundary_space, coastal_grid)
    CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)

    thermo_params = TDP.ThermodynamicsParameters(FT)
    surface_fluxes_params = SF.Parameters.SurfaceFluxesParameters(FT, SF.UniversalFunctions.BusingerParams)
    config = SF.SurfaceFluxConfig(
        SF.COARE3RoughnessParams{FT}(),
        SF.ConstantGustinessSpec(FT(1)),
    )

    # Uniform atmospheric and surface state: cool, moist air over a warmer
    # ocean with a westerly wind.
    state = (;
        T_atmos = FT(285),
        q_tot = FT(0.008),
        q_liq = FT(0),
        q_ice = FT(0),
        ρ_atmos = FT(1.2),
        u_atmos = FT(10),
        v_atmos = FT(0),
        height_int = FT(30),
        height_sfc = FT(0),
        T_sfc = FT(290),
    )

    @testset "UV <-> CT12 round trip" begin
        import Random
        rng = Random.Xoshiro(7)
        ct1 = CC.Fields.zeros(boundary_space)
        ct2 = CC.Fields.zeros(boundary_space)
        parent(ct1) .= randn(rng, FT, size(parent(ct1)))
        parent(ct2) .= randn(rng, FT, size(parent(ct2)))
        ct1_orig = copy(parent(ct1))
        ct2_orig = copy(parent(ct2))
        uv = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)
        CMIPExt.contravariant_to_cartesian!(uv, ct1, ct2)
        CMIPExt.cartesian_to_contravariant!(ct1, ct2, uv)
        @test parent(ct1) ≈ ct1_orig rtol = 1e-10
        @test parent(ct2) ≈ ct2_orig rtol = 1e-10
    end

    @testset "polygon kernel matches scalar SurfaceFluxes" begin
        fs = CMIPExt.ExchangeFluxState{FT}(OC.CPU(), eg.n_poly)
        fs.T_atmos .= state.T_atmos
        fs.q_tot .= state.q_tot
        fs.q_liq .= state.q_liq
        fs.q_ice .= state.q_ice
        fs.ρ_atmos .= state.ρ_atmos
        fs.u_atmos .= state.u_atmos
        fs.v_atmos .= state.v_atmos
        fs.height_int .= state.height_int
        fs.height_sfc .= state.height_sfc
        fs.T_sfc .= state.T_sfc
        fs.sic .= 0

        CMIPExt.compute_ocean_polygon_fluxes!(fs, surface_fluxes_params, thermo_params, config)
        fs.n_acc[] += 1

        reference = CMIPExt._polygon_surface_fluxes(
            surface_fluxes_params,
            thermo_params,
            config,
            state...,
        )
        @test all(fs.F_sh .== reference.F_sh)
        @test all(fs.F_lh .== reference.F_lh)
        @test all(fs.F_moisture .== reference.F_turb_moisture)
        @test all(fs.F_τu .== reference.F_turb_ρτxz)
        @test all(fs.F_τv .== reference.F_turb_ρτyz)

        # Physical sanity: warm ocean under cool air gives upward heat and
        # moisture fluxes; westerly wind gives negative (downward eastward
        # momentum) zonal stress and negligible meridional stress.
        @test reference.F_sh > 0
        @test reference.F_lh > 0
        @test reference.F_turb_moisture > 0
        @test reference.F_turb_ρτxz < 0
        @test abs(reference.F_turb_ρτyz) < abs(reference.F_turb_ρτxz) * 1e-6

        # Accumulators hold one contribution.
        @test fs.n_acc[] == 1
        @test all(fs.acc_F_sh .== fs.F_sh)

        @testset "scatter to boundary space" begin
            flux_scratch = (;
                F_turb_ρτxz = CC.Fields.zeros(boundary_space),
                F_turb_ρτyz = CC.Fields.zeros(boundary_space),
                F_sh = CC.Fields.zeros(boundary_space),
                F_lh = CC.Fields.zeros(boundary_space),
                F_turb_moisture = CC.Fields.zeros(boundary_space),
            )
            flux_dss_buffer = Utilities.init_dss_buffer(flux_scratch.F_sh)
            node_cov_dss = CC.Fields.zeros(boundary_space)
            CRExt.vec_to_se_field!(node_cov_dss, eg.node_cov)
            Utilities.apply_dss!(node_cov_dss, flux_dss_buffer)
            temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)
            weight_cov_scratch = CC.Fields.zeros(boundary_space)
            remapping = (;
                flux_scratch,
                flux_dss_buffer,
                node_cov_dss,
                weight_cov_scratch,
                temp_uv_vec,
                exchange_grid = eg,
            )

            # Open-water weighting (SIC = 0 here, so weight ≡ 1).
            fs.scratch2 .= 1 .- fs.sic
            CMIPExt.scatter_poly_fluxes_to_boundary!(remapping, eg, fs, fs.scratch2)

            # Scalar fluxes: uniform per-polygon values must be reproduced on
            # well-covered nodes and zeroed on uncovered ones.
            sh = CRExt.se_field_to_vec(flux_scratch.F_sh)
            covd = CRExt.se_field_to_vec(node_cov_dss)
            for n in 1:(eg.n_nodes)
                if covd[n] > 0.5
                    @test sh[n] ≈ reference.F_sh rtol = 1e-10
                elseif covd[n] <= 1e-3
                    @test sh[n] == 0
                end
            end

            # Momentum: converting the CT12 result back to UV must recover
            # the uniform UV stress on well-covered nodes.
            CMIPExt.contravariant_to_cartesian!(
                temp_uv_vec,
                flux_scratch.F_turb_ρτxz,
                flux_scratch.F_turb_ρτyz,
            )
            τu = CRExt.se_field_to_vec(temp_uv_vec.components.data.:1)
            τv = CRExt.se_field_to_vec(temp_uv_vec.components.data.:2)
            for n in 1:(eg.n_nodes)
                if covd[n] > 0.5
                    @test τu[n] ≈ reference.F_turb_ρτxz rtol = 1e-8
                    @test abs(τv[n]) < abs(reference.F_turb_ρτxz) * 1e-6
                end
            end
        end
    end
end

@testset "Per-polygon ice fluxes" begin
    eg = CMIPExt.build_exchange_grid(boundary_space, coastal_grid)
    thermo_params = TDP.ThermodynamicsParameters(FT)
    surface_fluxes_params = SF.Parameters.SurfaceFluxesParameters(FT, SF.UniversalFunctions.BusingerParams)
    config = SF.SurfaceFluxConfig(
        SF.ConstantRoughnessParams(FT(5.8e-5), FT(5.8e-5)),
        SF.ConstantGustinessSpec(FT(1)),
    )
    σ = FT(5.67e-8)
    ϵ = FT(1)
    α_albedo = FT(0.7)
    T_melt = FT(273.15)

    is = CMIPExt.IceExchangeState{FT}(OC.CPU(), eg.n_poly)
    fs = is.fluxes
    # Cold air over ice near the melting point; weak wind; half the polygons
    # ice covered.
    fs.T_atmos .= FT(260)
    fs.q_tot .= FT(0.002)
    fs.q_liq .= 0
    fs.q_ice .= 0
    fs.ρ_atmos .= FT(1.3)
    fs.u_atmos .= FT(5)
    fs.v_atmos .= FT(0)
    fs.height_int .= FT(30)
    fs.height_sfc .= FT(0)
    fs.T_sfc .= FT(268)
    fs.sic .= 0
    fs.sic[1:2:end] .= FT(0.8)
    is.R .= FT(1.0) / FT(2.0) # 1 m of ice at conductivity 2 W/m/K
    is.T_i .= FT(271.2)
    is.SW_d .= FT(50)
    is.LW_d .= FT(200)
    is.T_sfc_new .= 0

    CMIPExt.compute_ice_polygon_fluxes!(
        is,
        surface_fluxes_params,
        thermo_params,
        config,
        σ,
        ϵ,
        α_albedo,
        T_melt,
    )
    fs.n_acc[] += 1

    # Ice-free polygons short-circuit to zero flux and keep the T_sfc guess.
    for k in 1:(eg.n_poly)
        if fs.sic[k] == 0
            @test fs.F_sh[k] == 0
            @test fs.F_lh[k] == 0
            @test fs.F_moisture[k] == 0
            @test fs.F_τu[k] == 0
            @test is.T_sfc_new[k] == fs.T_sfc[k]
        else
            # The diagnosed surface temperature is bounded by melt and stays
            # physical; the diagnosed fluxes are finite.
            @test is.T_sfc_new[k] <= T_melt + sqrt(eps(FT))
            @test FT(200) < is.T_sfc_new[k] < FT(280)
            @test isfinite(fs.F_sh[k]) && isfinite(fs.F_lh[k])
            @test fs.F_τu[k] < 0 # westerly wind drag
        end
    end
    # Uniform inputs: all icy polygons agree exactly.
    icy = findall(>(0), fs.sic)
    @test length(unique(fs.F_sh[icy])) == 1
    @test length(unique(is.T_sfc_new[icy])) == 1
    @test all(fs.acc_F_sh .== fs.F_sh)

    @testset "ice-weighted scatter" begin
        flux_scratch = (;
            F_turb_ρτxz = CC.Fields.zeros(boundary_space),
            F_turb_ρτyz = CC.Fields.zeros(boundary_space),
            F_sh = CC.Fields.zeros(boundary_space),
            F_lh = CC.Fields.zeros(boundary_space),
            F_turb_moisture = CC.Fields.zeros(boundary_space),
        )
        flux_dss_buffer = Utilities.init_dss_buffer(flux_scratch.F_sh)
        weight_cov_scratch = CC.Fields.zeros(boundary_space)
        temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)
        remapping = (;
            flux_scratch,
            flux_dss_buffer,
            weight_cov_scratch,
            temp_uv_vec,
            exchange_grid = eg,
        )
        CMIPExt.scatter_poly_fluxes_to_boundary!(remapping, eg, fs, fs.sic)

        # The sic-weighted average of a field that is uniform on icy polygons
        # is exactly that value wherever any ice coverage exists.
        CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
        sh = CRExt.se_field_to_vec(flux_scratch.F_sh)
        cov = CRExt.se_field_to_vec(weight_cov_scratch)
        F_ice = fs.F_sh[icy[1]]
        for n in 1:(eg.n_nodes)
            if cov[n] > 0.1
                @test sh[n] ≈ F_ice rtol = 1e-8
            elseif cov[n] <= 1e-3
                @test sh[n] == 0
            end
        end
    end
end

@testset "Full ocean flux driver allocations (CPU)" begin
    eg = CMIPExt.build_exchange_grid(boundary_space, coastal_grid)
    CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
    thermo_params = TDP.ThermodynamicsParameters(FT)
    surface_fluxes_params = SF.Parameters.SurfaceFluxesParameters(FT, SF.UniversalFunctions.BusingerParams)
    config = SF.SurfaceFluxConfig(
        SF.COARE3RoughnessParams{FT}(),
        SF.ConstantGustinessSpec(FT(1)),
    )

    # Fake coupler fields and remapping bundle with everything the driver
    # touches.
    make_field(v) = (f = CC.Fields.zeros(boundary_space); f .= FT(v); f)
    csf = (;
        u_int = make_field(8),
        v_int = make_field(2),
        T_atmos = make_field(285),
        q_tot_atmos = make_field(0.008),
        q_liq_atmos = make_field(0),
        q_ice_atmos = make_field(0),
        ρ_atmos = make_field(1.2),
        height_int = make_field(30),
        height_sfc = make_field(0),
    )
    flux_scratch = (;
        F_turb_ρτxz = CC.Fields.zeros(boundary_space),
        F_turb_ρτyz = CC.Fields.zeros(boundary_space),
        F_sh = CC.Fields.zeros(boundary_space),
        F_lh = CC.Fields.zeros(boundary_space),
        F_turb_moisture = CC.Fields.zeros(boundary_space),
    )
    remapping = (;
        flux_scratch,
        flux_dss_buffer = Utilities.init_dss_buffer(flux_scratch.F_sh),
        weight_cov_scratch = CC.Fields.zeros(boundary_space),
        temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space),
        exchange_grid = eg,
    )
    fs = CMIPExt.ExchangeFluxState{FT}(OC.CPU(), eg.n_poly)
    fs.T_sfc .= FT(290)
    fs.sic .= 0

    function ocean_driver!()
        CMIPExt.gather_atmos_state_to_polys!(fs, eg, csf, remapping.temp_uv_vec)
        CMIPExt.compute_ocean_polygon_fluxes!(
            fs,
            surface_fluxes_params,
            thermo_params,
            config,
        )
        fs.scratch2 .= 1 .- fs.sic
        CMIPExt.scatter_poly_fluxes_to_boundary!(remapping, eg, fs, fs.scratch2)
        return nothing
    end

    ocean_driver!()
    allocated = @allocated ocean_driver!()
    # The steady-state driver may only allocate small wrappers (array views,
    # broadcast objects) — no O(n_poly) or O(n_nodes) buffers.
    @test allocated < 50_000
    @test all(isfinite, parent(flux_scratch.F_sh))
end
