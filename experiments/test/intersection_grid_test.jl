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
to_cpu(x::Array) = x
to_cpu(x) = Array(x)
exchange_arch() = ClimaComms.device() isa ClimaComms.CUDADevice ? OC.GPU() : OC.CPU()
on_exchange_device(x) = CMIPExt.on_device(exchange_arch(), x)

radius = FT(6.371e6)
boundary_space =
    CC.CommonSpaces.CubedSphereSpace(FT; radius, n_quad_points = 4, h_elem = 4, context)
Nq = 4
n_elem = 6 * 4^2

# Small `TripolarGrid` matching the production layout: displaced north poles
# at 55°N, regular southern boundary at 80°S, folded top row.
Nx, Ny, Nz = 36, 18, 1
arch = OC.CPU()
tripolar_grid_on(arch) = OC.TripolarGrid(
    arch;
    size = (Nx, Ny, Nz),
    southernmost_latitude = -80,
    north_poles_latitude = 55,
    first_pole_longitude = 70,
    fold_topology = OC.RightCenterFolded,
    z = (-100.0, 0.0),
    halo = (4, 4, 4),
)
underlying_tripolar = tripolar_grid_on(arch)

# Fully wet grid: bottom below the deepest z everywhere.
wet_grid = OC.ImmersedBoundaryGrid(
    underlying_tripolar,
    OC.GridFittedBottom((x, y) -> -200.0);
    active_cells_map = false,
)

# Half-dry grid: an idealized land hemisphere. TripolarGrid longitudes are
# not restricted to [-180, 180), so wrap before comparing. Built per
# architecture so a test that mixes ocean fields with coupler fields can put
# both in the same memory space.
coastal_bottom(x, y) = mod(x, 360) < 180 ? 100.0 : -200.0
coastal_grid_on(arch) = OC.ImmersedBoundaryGrid(
    tripolar_grid_on(arch),
    OC.GridFittedBottom(coastal_bottom);
    active_cells_map = false,
)
coastal_grid = coastal_grid_on(arch)
uv_basis = CMIPExt.uv_basis_coefficients(boundary_space)

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
        Jw = to_cpu(CRExt.se_node_weights(boundary_space))
        coords = CC.Fields.coordinate_field(boundary_space)
        node_lat = to_cpu(CRExt.flat_nodal_data(CC.Fields.field_values(coords.lat)))
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

        # The scatter weights are the non-negative control-volume weights.
        @test all(>=(0), eg.sweight)

        # On a fully wet grid every polygon is retained. The control volumes are
        # sized so the weights sum to exactly 1 over the polygons tiling a node's
        # control volume, so the coverage is 1 at every node of an element the
        # ocean grid fully covers (north of its 80°S boundary).
        CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
        coords = CC.Fields.coordinate_field(boundary_space)
        node_lat = to_cpu(CRExt.flat_nodal_data(CC.Fields.field_values(coords.lat)))
        elem_min_lat =
            [minimum(node_lat[((e - 1) * Nq ^ 2 + 1):(e * Nq ^ 2)]) for e in 1:n_elem]
        for n in 1:(eg.n_nodes)
            if elem_min_lat[(n - 1) ÷ Nq ^ 2 + 1] > -75
                @test eg.node_cov[n] ≈ 1 atol = 1e-8
            end
            @test -1e-12 <= eg.node_cov[n] <= 1 + 1e-8
        end
        # Nodes south of the ocean grid boundary are only partially covered.
        @test any(n -> node_lat[n] < -80 && eg.node_cov[n] < 0.9, 1:(eg.n_nodes))

        # Scatter column sums (the control-volume polygon areas) match the
        # geometric areas closely.
        @test eg.b_area ≈ eg.area rtol = 0.05

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

        # The wet fraction is the wet control-volume coverage: 1 at nodes whose
        # control volume lies entirely over wet cells, exactly 0 where it lies
        # entirely over dry cells, intermediate at the coastline. The weights
        # are non-negative, so it never leaves [0, 1].
        wet_frac = eg_c.node_cov
        @test count(f -> isapprox(f, 1; atol = 1e-8), wet_frac) > 0.2 * length(wet_frac)
        @test count(==(0), wet_frac) > 0.2 * length(wet_frac)
        @test any(f -> 0.05 < f < 0.95, wet_frac)
        @test all(f -> 0 <= f <= 1 + 1e-8, wet_frac)
    end
end

import Random

@testset "Gather/scatter operations" begin
    eg = CMIPExt.build_exchange_grid(boundary_space, coastal_grid)
    CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
    Jw = to_cpu(CRExt.se_node_weights(boundary_space))
    rng = Random.Xoshiro(42)

    poly_scratch = zeros(FT, eg.n_poly)
    nodal_scratch = zeros(FT, eg.n_nodes)
    cell_scratch = zeros(FT, eg.n_oc)

    @testset "constant preservation" begin
        # Gather: rows sum to 1, so a constant nodal field is exact.
        CMIPExt.gather_nodes_to_polys!(poly_scratch, eg, ones(FT, eg.n_nodes))
        @test all(v -> isapprox(v, 1, atol = 1e-12), poly_scratch)
        CMIPExt.gather_nodes_to_polys!(poly_scratch, eg, ones(FT, eg.n_nodes), eg.tweight)
        @test all(v -> isapprox(v, 1, atol = 1e-12), poly_scratch)

        # Raw scatter of a constant reproduces the nodal wet coverage.
        CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, ones(FT, eg.n_poly))
        @test nodal_scratch ≈ eg.node_cov rtol = 1e-12

        # Normalized scatter of a constant is exactly the constant on
        # covered nodes and 0 elsewhere.
        CMIPExt.scatter_polys_to_nodes_normalized!(nodal_scratch, eg, ones(FT, eg.n_poly))
        for n in 1:(eg.n_nodes)
            if eg.node_cov[n] > 0
                @test nodal_scratch[n] ≈ 1 rtol = 1e-12
            else
                @test nodal_scratch[n] == 0
            end
        end

        # Area-weighted cell mean of a constant is exact on wet cells, 0 on
        # dry/shadow cells; the fold mirror then fills shadow partners.
        CMIPExt.scatter_polys_to_cells!(
            cell_scratch,
            eg,
            ones(FT, eg.n_poly),
            eg.area,
            eg.oc_wet_area,
        )
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

        # FV side: integrated with the FV grid's own cell areas, the scatter
        # gives the SE-side integral exactly.
        CMIPExt.scatter_polys_to_cells!(cell_scratch, eg, F)
        grid_cpu = OC.on_architecture(OC.CPU(), coastal_grid.underlying_grid)
        Az = [
            OC.Operators.Azᶜᶜᶜ(mod1(c, Nx), (c - 1) ÷ Nx + 1, Nz, grid_cpu) for
            c in 1:(eg.n_oc)
        ]
        @test sum(Az .* cell_scratch) ≈ sum(eg.b_area .* F) rtol = 1e-12

        # SE side: the scatter preserves the integral over the control-volume
        # areas exactly (before DSS).
        CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, F)
        @test sum(Jw .* nodal_scratch) ≈ sum(eg.b_area .* F) rtol = 1e-12

        # And the control-volume areas agree with the geometric areas.
        @test sum(eg.b_area) ≈ sum(eg.area) rtol = 1e-3

        # Transpose gather: weighted by any w, the polygons receive the
        # atmosphere's share Σ_n Jw_n α_n g_n, with α the scatter of w.
        g, w = rand(rng, FT, eg.n_nodes), rand(rng, FT, eg.n_poly)
        CMIPExt.gather_nodes_to_polys!(poly_scratch, eg, g, eg.tweight)
        CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, w)
        atmosphere_share = sum(Jw .* nodal_scratch .* g)
        @test sum(eg.b_area .* w .* poly_scratch) ≈ atmosphere_share rtol = 1e-12
    end

    @testset "bounded normalized scatter" begin
        # A masked, normalized transfer must stay within the range of its
        # inputs; a constant-field test cannot detect a violation. Values in
        # [-1, 1], unweighted and with a random (sea-ice-like) weight.
        F = 2 .* rand(rng, FT, eg.n_poly) .- 1
        CMIPExt.scatter_polys_to_nodes_normalized!(nodal_scratch, eg, F)
        @test all(v -> -1 - 1e-12 <= v <= 1 + 1e-12, nodal_scratch)
        w = rand(rng, FT, eg.n_poly)
        num, den = zeros(FT, eg.n_nodes), zeros(FT, eg.n_nodes)
        CMIPExt.scatter_polys_to_nodes!(num, eg, w .* F)
        CMIPExt.scatter_polys_to_nodes!(den, eg, w)
        @test all(
            n -> den[n] == 0 || -1 - 1e-12 <= num[n] / den[n] <= 1 + 1e-12,
            1:(eg.n_nodes),
        )

        # The transpose gather is a convex combination of nodal values.
        g = 2 .* rand(rng, FT, eg.n_nodes) .- 1
        CMIPExt.gather_nodes_to_polys!(poly_scratch, eg, g, eg.tweight)
        @test all(v -> -1 - 1e-12 <= v <= 1 + 1e-12, poly_scratch)
    end

    @testset "no O(n) allocations (CPU)" begin
        # KernelAbstractions CPU launches allocate small constant-size
        # wrappers; the invariant is that no O(n_poly)/O(n_nodes) buffers
        # appear (those would be ≥ 8 * n ≳ 100 kB here).
        alloc_budget = 10_000
        nodal = rand(rng, FT, eg.n_nodes)
        polyv = rand(rng, FT, eg.n_poly)
        CMIPExt.gather_nodes_to_polys!(poly_scratch, eg, nodal)
        CMIPExt.gather_cells_to_polys!(poly_scratch, eg, cell_scratch)
        CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, polyv)
        CMIPExt.scatter_polys_to_nodes_normalized!(nodal_scratch, eg, polyv)
        CMIPExt.scatter_polys_to_cells!(cell_scratch, eg, polyv)
        CMIPExt.mirror_fold_partners!(cell_scratch, coastal_grid)
        @test @allocated(CMIPExt.gather_nodes_to_polys!(poly_scratch, eg, nodal)) <
              alloc_budget
        @test @allocated(CMIPExt.gather_cells_to_polys!(poly_scratch, eg, cell_scratch)) <
              alloc_budget
        @test @allocated(CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, polyv)) <
              alloc_budget
        @test @allocated(
            CMIPExt.scatter_polys_to_nodes_normalized!(nodal_scratch, eg, polyv)
        ) < alloc_budget
        @test @allocated(CMIPExt.scatter_polys_to_cells!(cell_scratch, eg, polyv)) <
              alloc_budget
        @test @allocated(CMIPExt.mirror_fold_partners!(cell_scratch, coastal_grid)) <
              alloc_budget
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
            @test Array(
                CMIPExt.gather_nodes_to_polys!(to_d(poly_scratch), eg_d, to_d(nodal)),
            ) ≈ poly_scratch rtol = rt

            CMIPExt.gather_cells_to_polys!(poly_scratch, eg, cellv)
            @test Array(
                CMIPExt.gather_cells_to_polys!(to_d(poly_scratch), eg_d, to_d(cellv)),
            ) == poly_scratch

            CMIPExt.scatter_polys_to_nodes!(nodal_scratch, eg, polyv)
            @test Array(
                CMIPExt.scatter_polys_to_nodes!(to_d(nodal_scratch), eg_d, to_d(polyv)),
            ) ≈ nodal_scratch rtol = rt

            CMIPExt.scatter_polys_to_nodes_normalized!(nodal_scratch, eg, polyv)
            @test Array(
                CMIPExt.scatter_polys_to_nodes_normalized!(
                    to_d(nodal_scratch),
                    eg_d,
                    to_d(polyv),
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
    node_lat = to_cpu(CRExt.flat_nodal_data(CC.Fields.field_values(coords.lat)))
    node_lon = to_cpu(CRExt.flat_nodal_data(CC.Fields.field_values(coords.long)))

    f = CMIPExt.wet_ocean_fraction_field(boundary_space, eg_c)
    vals = to_cpu(CRExt.se_field_to_vec(f))

    # Bounds and non-triviality.
    @test all(v -> 0 <= v <= 1, vals)
    @test any(!=(0), vals) && any(!=(1), vals)

    # Interior values: fully ocean far from the idealized coastlines at
    # wrapped longitudes 0 and 180, fully land in the middle of the land
    # hemisphere. (Land: mod(lon, 360) ∈ [0, 180).)
    for n in 1:(eg_c.n_nodes)
        lonw = mod(node_lon[n], 360)
        abs(node_lat[n]) > 45 && continue
        if 60 < lonw < 120 # deep inside land
            @test vals[n] < 0.1
        elseif 240 < lonw < 300 # deep inside ocean
            @test vals[n] > 0.9
        end
    end

    # Fully wet grid: fraction ≈ 1 everywhere the ocean grid reaches.
    f_wet = CMIPExt.wet_ocean_fraction_field(boundary_space, eg_wet)
    wet_vals = to_cpu(CRExt.se_field_to_vec(f_wet))
    for n in 1:(eg_wet.n_nodes)
        if node_lat[n] > -60
            @test wet_vals[n] > 0.95
        end
    end
end

import SurfaceFluxes as SF
import Thermodynamics.Parameters as TDP
import ClimaParams as CP
import ClimaCoupler: FieldExchanger, FluxCalculator, Interfacer, Utilities

@testset "Per-polygon ocean fluxes" begin
    eg = CMIPExt.build_exchange_grid(boundary_space, coastal_grid)
    CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)

    thermo_params = TDP.ThermodynamicsParameters(FT)
    surface_fluxes_params =
        SF.Parameters.SurfaceFluxesParameters(FT, SF.UniversalFunctions.BusingerParams)
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
        copyto!(parent(ct1), randn(rng, FT, size(parent(ct1))))
        copyto!(parent(ct2), randn(rng, FT, size(parent(ct2))))
        ct1_orig = to_cpu(parent(ct1))
        ct2_orig = to_cpu(parent(ct2))
        uv = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)
        CMIPExt.contravariant_to_cartesian!(uv, ct1, ct2)
        CMIPExt.cartesian_to_contravariant!(ct1, ct2, uv)
        @test to_cpu(parent(ct1)) ≈ ct1_orig rtol = 1e-10
        @test to_cpu(parent(ct2)) ≈ ct2_orig rtol = 1e-10
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

        CMIPExt.compute_ocean_polygon_fluxes!(
            fs,
            surface_fluxes_params,
            thermo_params,
            config,
        )

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
            # Fields live on `boundary_space`'s device; gather/scatter kernels
            # take the destination's backend, so the exchange grid and flux
            # state must match (GPU CI uses CUDA Fields).
            eg_d = on_exchange_device(eg)
            fs_d = on_exchange_device(fs)
            flux_scratch = (;
                F_turb_ρτxz = CC.Fields.zeros(boundary_space),
                F_turb_ρτyz = CC.Fields.zeros(boundary_space),
                F_sh = CC.Fields.zeros(boundary_space),
                F_lh = CC.Fields.zeros(boundary_space),
                F_turb_moisture = CC.Fields.zeros(boundary_space),
            )
            flux_dss_buffer = Utilities.init_dss_buffer(flux_scratch.F_sh)
            node_cov_dss = CC.Fields.zeros(boundary_space)
            CRExt.vec_to_se_field!(node_cov_dss, on_exchange_device(eg.node_cov))
            Utilities.apply_dss!(node_cov_dss, flux_dss_buffer)
            temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)
            weight_cov_scratch = CC.Fields.zeros(boundary_space)
            remapping = (;
                flux_scratch,
                flux_dss_buffer,
                node_cov_dss,
                weight_cov_scratch,
                temp_uv_vec,
                uv_basis,
                exchange_grid = eg_d,
            )

            # Open-water weighting (SIC = 0 here, so weight ≡ 1).
            fs_d.scratch2 .= 1 .- fs_d.sic
            CMIPExt.scatter_poly_fluxes_to_boundary!(remapping, eg_d, fs_d, fs_d.scratch2)

            # The flux-weight coverage is the DSS'd wet coverage (weight ≡ 1).
            covd = to_cpu(CRExt.se_field_to_vec(node_cov_dss))
            @test to_cpu(CRExt.se_field_to_vec(weight_cov_scratch)) ≈ covd rtol = 1e-12

            # Scalar fluxes: the normalized scatter is a convex combination,
            # so uniform per-polygon values are reproduced at every node with
            # any coverage and zeroed at uncovered ones.
            sh = to_cpu(CRExt.se_field_to_vec(flux_scratch.F_sh))
            for n in 1:(eg.n_nodes)
                if covd[n] > 0
                    @test sh[n] ≈ reference.F_sh rtol = 1e-10
                else
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
            τu = to_cpu(CRExt.se_field_to_vec(temp_uv_vec.components.data.:1))
            τv = to_cpu(CRExt.se_field_to_vec(temp_uv_vec.components.data.:2))
            for n in 1:(eg.n_nodes)
                if covd[n] > 0
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
    surface_fluxes_params =
        SF.Parameters.SurfaceFluxesParameters(FT, SF.UniversalFunctions.BusingerParams)
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
    fs.SW_d .= FT(50)
    fs.LW_d .= FT(200)
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
        eg_d = on_exchange_device(eg)
        fs_d = on_exchange_device(fs)
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
            uv_basis,
            exchange_grid = eg_d,
        )
        CMIPExt.scatter_poly_fluxes_to_boundary!(remapping, eg_d, fs_d, fs_d.sic)

        # The sic-weighted average of a field that is uniform on icy polygons
        # is exactly that value wherever any ice coverage exists.
        CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
        sh = to_cpu(CRExt.se_field_to_vec(flux_scratch.F_sh))
        cov = to_cpu(CRExt.se_field_to_vec(weight_cov_scratch))
        F_ice = fs.F_sh[icy[1]]
        for n in 1:(eg.n_nodes)
            if cov[n] > 0
                @test sh[n] ≈ F_ice rtol = 1e-8
            else
                @test sh[n] == 0
            end
        end
    end
end

@testset "Area fraction × flux conserves the polygon integral" begin
    # The coupler hands the atmosphere `area_fraction × flux`. With the area
    # fraction equal to the coverage of the flux weight, that product integrates
    # exactly to the polygon flux integral over the control-volume areas, and
    # every nodal flux stays within the range of the polygon fluxes.
    eg = CMIPExt.build_exchange_grid(boundary_space, coastal_grid)
    eg_d = on_exchange_device(eg)
    CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
    Jw = to_cpu(CRExt.se_node_weights(boundary_space))
    rng = Random.Xoshiro(11)

    fs = CMIPExt.ExchangeFluxState{FT}(OC.CPU(), eg.n_poly)
    for f in (fs.F_sh, fs.F_lh, fs.F_moisture, fs.F_τu, fs.F_τv)
        f .= 100 .+ 50 .* (2 .* rand(rng, FT, eg.n_poly) .- 1)
    end
    # Ice-free, partly and fully ice-covered polygons.
    fs.sic .= clamp.(2 .* rand(rng, FT, eg.n_poly) .- 0.5, 0, 1)
    fs_d = on_exchange_device(fs)

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
        downward_flux_scratch = CC.Fields.zeros(boundary_space),
        temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space),
        uv_basis,
        exchange_grid = eg_d,
    )
    coverage = CC.Fields.zeros(boundary_space)

    for weight in (1 .- fs.sic, copy(fs.sic)) # open water, sea ice
        weight_d = on_exchange_device(weight)
        CMIPExt.scatter_poly_fluxes_to_boundary!(remapping, eg_d, fs_d, weight_d)
        α = to_cpu(CRExt.se_field_to_vec(remapping.weight_cov_scratch))
        @test all(a -> 0 <= a <= 1 + 1e-8, α)

        # A downward flux whose element copies of a shared node disagree, as the
        # atmosphere's radiation can, still reaches the polygons as exactly the
        # atmosphere's share Σ_n Jw_n α_n g_n.
        g = rand(rng, FT, eg.n_nodes)
        g_field = CC.Fields.zeros(boundary_space)
        CRExt.vec_to_se_field!(g_field, on_exchange_device(g))
        poly_d = on_exchange_device(zeros(FT, eg.n_poly))
        CMIPExt.gather_downward_flux!(poly_d, remapping, g_field)
        @test sum(eg.b_area .* weight .* to_cpu(poly_d)) ≈ sum(Jw .* α .* g) rtol = 1e-12

        for (name, F_poly) in
            ((:F_sh, fs.F_sh), (:F_lh, fs.F_lh), (:F_turb_moisture, fs.F_moisture))
            F = to_cpu(CRExt.se_field_to_vec(getproperty(flux_scratch, name)))
            @test sum(Jw .* α .* F) ≈ sum(eg.b_area .* weight .* F_poly) rtol = 1e-10
            @test sum(Jw .* α .* F) ≈ sum(eg.area .* weight .* F_poly) rtol = 1e-3
            @test all(v -> v == 0 || 50 - 1e-8 <= v <= 150 + 1e-8, F)
            @test all(n -> α[n] == 0 || F[n] != 0, eachindex(F))
        end

        # The emission temperature: area fraction × T⁴ integrates to the
        # polygon integral of T⁴, and T stays within the polygon range.
        T⁴ = (250 .+ 50 .* rand(rng, FT, eg.n_poly)) .^ 4
        emitted = sum(eg.b_area .* weight .* T⁴)
        T_field = CC.Fields.zeros(boundary_space)
        CMIPExt.emission_temperature!(T_field, remapping, on_exchange_device(T⁴), weight_d)
        T = to_cpu(CRExt.se_field_to_vec(T_field))
        @test sum(Jw .* α .* T .^ 4) ≈ emitted rtol = 1e-10
        @test all(n -> α[n] == 0 || 250 - 1e-8 <= T[n] <= 300 + 1e-8, eachindex(T))
    end

    # The open-ocean fraction the coupler builds, wet − ice, is the coverage of
    # `1 - sic`: the denominator of the open-ocean flux scatter.
    wet =
        to_cpu(CRExt.se_field_to_vec(CMIPExt.wet_ocean_fraction_field(boundary_space, eg)))
    buffer = remapping.flux_dss_buffer
    function coverage_of(weight)
        CMIPExt.flux_weight_coverage!(coverage, eg_d, on_exchange_device(weight), buffer)
        return copy(to_cpu(CRExt.se_field_to_vec(coverage)))
    end
    @test wet .- coverage_of(copy(fs.sic)) ≈ coverage_of(1 .- fs.sic) atol = 1e-12
end

@testset "Full ocean flux driver allocations (CPU)" begin
    eg_cpu = CMIPExt.build_exchange_grid(boundary_space, coastal_grid)
    # Gather/scatter kernels follow the destination backend; Fields are on
    # `boundary_space`'s device, so the exchange grid and flux state match it.
    eg = on_exchange_device(eg_cpu)
    thermo_params = TDP.ThermodynamicsParameters(FT)
    surface_fluxes_params =
        SF.Parameters.SurfaceFluxesParameters(FT, SF.UniversalFunctions.BusingerParams)
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
        uv_basis,
        exchange_grid = eg,
    )
    fs = CMIPExt.ExchangeFluxState{FT}(exchange_arch(), eg.n_poly)
    fs.T_sfc .= FT(290)
    fs.sic .= 0

    function ocean_driver!()
        CMIPExt.gather_atmos_state_to_polys!(
            fs,
            eg,
            csf,
            remapping.temp_uv_vec,
            remapping.uv_basis,
        )
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
    if exchange_arch() isa OC.CPU
        allocated = @allocated ocean_driver!()
        @test allocated < 100_000
    end
    @test all(isfinite, to_cpu(parent(flux_scratch.F_sh)))
end

# Stand-in for an Oceananigans model in the `update_sim!` radiation path, which
# reads only the grid, the tracers (and their top flux BCs) and the float type.
struct FakeOceanModel{G, T}
    grid::G
    tracers::T
end
Base.eltype(::FakeOceanModel) = FT

# Σ_{k ∈ c} b_area_k / Az_c with the grid's own cell area: the factor the cell
# scatter applies to a uniform flux.
function cell_area_ratio(eg, grid)
    grid_cpu = OC.on_architecture(OC.CPU(), grid.underlying_grid)
    ratio = zeros(FT, eg.n_oc)
    for k in 1:(eg.n_poly)
        ratio[eg.oc_of_poly[k]] += eg.b_area[k]
    end
    for c in 1:(eg.n_oc)
        ratio[c] /= OC.Operators.Azᶜᶜᶜ(mod1(c, Nx), (c - 1) ÷ Nx + 1, Nz, grid_cpu)
    end
    return ratio
end

@testset "Exchange-grid radiation and precipitation" begin
    # Make sure both the boundary space and the ocean grid are on the same device
    ocean_arch = exchange_arch()
    grid = coastal_grid_on(ocean_arch)
    eg_cpu = CMIPExt.build_exchange_grid(boundary_space, grid)
    eg = on_exchange_device(eg_cpu)

    # Only the tracers, their top flux BCs, the grid and `eltype` are read, so
    # a full `HydrostaticFreeSurfaceModel` is unnecessary.
    tracer(flux) = OC.CenterField(
        grid,
        boundary_conditions = OC.FieldBoundaryConditions(
            grid,
            (OC.Center(), OC.Center(), OC.Center());
            top = OC.FluxBoundaryCondition(flux),
        ),
    )
    surface_field() = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    T_flux, S_flux = surface_field(), surface_field()
    model = FakeOceanModel(grid, (; T = tracer(T_flux), S = tracer(S_flux)))
    T_sfc, S_sfc = FT(12), FT(35) # [°C], [psu]
    OC.set!(model.tracers.T, T_sfc)
    OC.set!(model.tracers.S, S_sfc)

    ocean_properties = (;
        reference_density = FT(1020),
        heat_capacity = FT(3995),
        σ = FT(5.67e-8),
        C_to_K = FT(273.15),
    )
    remapping = (;
        use_exchange_grid = true,
        exchange_grid = eg,
        ocean_flux_state = CMIPExt.ExchangeFluxState{FT}(ocean_arch, eg_cpu.n_poly),
        scratch_field_oc3 = surface_field(),
        weight_cov_scratch = CC.Fields.zeros(boundary_space),
        downward_flux_scratch = CC.Fields.zeros(boundary_space),
        flux_dss_buffer = Utilities.init_dss_buffer(CC.Fields.zeros(boundary_space)),
        exchange_fields = (;
            SW_d = CC.Fields.zeros(boundary_space),
            LW_d = CC.Fields.zeros(boundary_space),
            P_liq = CC.Fields.zeros(boundary_space),
            P_snow = CC.Fields.zeros(boundary_space),
            scalar_temp = CC.Fields.zeros(boundary_space),
        ),
    )
    sim = CMIPExt.OceananigansSimulation(
        (; model),
        nothing,
        ocean_properties,
        remapping,
        surface_field(), # ice concentration: ice-free
        nothing,
        nothing,
    )

    # Uniform forcing: the node gather preserves constants exactly, so each wet
    # cell receives the single-column answer times its `cell_area_ratio`.
    uniform(v) = (f = CC.Fields.zeros(boundary_space); f .= FT(v); f)
    SW_d, LW_d, P_liq, P_snow = FT(300), FT(350), FT(-1e-5), FT(-2e-6)
    csf = (;
        SW_d = uniform(SW_d),
        LW_d = uniform(LW_d),
        P_liq = uniform(P_liq),
        P_snow = uniform(P_snow),
    )
    FieldExchanger.update_sim!(sim, csf)

    (; reference_density, heat_capacity, σ, C_to_K) = ocean_properties
    α = FT(0.011) # `Val(:surface_direct_albedo)`
    ϵ = FT(0.97)  # `Val(:emissivity)`
    expected_T =
        (-(1 - α) * SW_d - ϵ * (LW_d - σ * (T_sfc + C_to_K)^4)) /
        (reference_density * heat_capacity)
    expected_S = -S_sfc * (P_liq + P_snow) / reference_density

    cell_values(f) = vec(OC.interior(OC.on_architecture(OC.CPU(), f), :, :, 1))
    flux_T = cell_values(CMIPExt.surface_flux(model.tracers.T))
    flux_S = cell_values(CMIPExt.surface_flux(model.tracers.S))
    r = cell_area_ratio(eg_cpu, grid)
    n_wet = 0
    # Skip the fold row, whose shadow cells are filled by mirroring rather than
    # by the scatter.
    for c in 1:((Ny - 1) * Nx)
        if eg_cpu.oc_wet_area[c] > 0
            n_wet += 1
            @test flux_T[c] ≈ expected_T * r[c] rtol = 1e-12
            @test flux_S[c] ≈ expected_S * r[c] rtol = 1e-12
        else
            # Dry cells receive nothing: the property this path exists for.
            @test flux_T[c] == 0
            @test flux_S[c] == 0
        end
    end
    @test n_wet > 0

    # The emission temperature: the open-ocean area fraction times T⁴ integrates to
    # the polygon integral of (1 − SIC) T⁴ of the top-layer temperature in kelvin.
    OC.set!(model.tracers.T, (λ, φ, z) -> 10 + 10 * cosd(φ))
    OC.set!(sim.ice_concentration, (λ, φ) -> φ < 0 ? FT(0.3) : FT(0))
    T_field = CC.Fields.zeros(boundary_space)
    Interfacer.get_field!(T_field, sim, Val(:surface_temperature))
    CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
    Jw = to_cpu(CRExt.se_node_weights(boundary_space))
    coverage = to_cpu(CRExt.se_field_to_vec(remapping.weight_cov_scratch))
    T = to_cpu(CRExt.se_field_to_vec(T_field))
    T_cells = vec(OC.interior(OC.on_architecture(OC.CPU(), model.tracers.T), :, :, Nz))
    k_cell = eg_cpu.oc_of_poly
    sic = cell_values(sim.ice_concentration)[k_cell]
    emitted = sum(eg_cpu.b_area .* (1 .- sic) .* (T_cells[k_cell] .+ C_to_K) .^ 4)
    @test sum(Jw .* coverage .* T .^ 4) ≈ emitted rtol = 1e-10
end

# Stand-in for a ClimaSeaIce model in the `update_sim!` exchange-grid path
struct FakeIceModel{G, C, H, S}
    grid::G
    ice_concentration::C
    external_heat_fluxes::H
    snowfall::S
end
Base.eltype(::FakeIceModel) = FT

@testset "Exchange-grid radiation and snowfall on sea ice" begin
    ice_arch = exchange_arch()
    grid = coastal_grid_on(ice_arch)
    eg_cpu = CMIPExt.build_exchange_grid(boundary_space, grid)
    eg = on_exchange_device(eg_cpu)

    surface_field() = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    ice_concentration, top, snowfall = surface_field(), surface_field(), surface_field()
    OC.set!(ice_concentration, (λ, φ) -> φ < 0 ? FT(0.8) : FT(0))
    model = FakeIceModel(grid, ice_concentration, (; top), snowfall)
    remapping = (;
        use_exchange_grid = true,
        exchange_grid = eg,
        ice_flux_state = CMIPExt.IceExchangeState{FT}(ice_arch, eg_cpu.n_poly),
        scratch_field_oc1 = surface_field(),
        weight_cov_scratch = CC.Fields.zeros(boundary_space),
        downward_flux_scratch = CC.Fields.zeros(boundary_space),
        flux_dss_buffer = Utilities.init_dss_buffer(CC.Fields.zeros(boundary_space)),
        exchange_fields = (;
            SW_d = CC.Fields.zeros(boundary_space),
            LW_d = CC.Fields.zeros(boundary_space),
            P_snow = CC.Fields.zeros(boundary_space),
            scalar_temp = CC.Fields.zeros(boundary_space),
        ),
    )
    sim = CMIPExt.ClimaSeaIceSimulation(
        (; model),
        nothing,
        remapping,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )

    uniform(v) = (f = CC.Fields.zeros(boundary_space); f .= FT(v); f)
    SW_d, LW_d, P_snow = FT(300), FT(250), FT(-2e-6) # P_snow is downward
    csf = (; SW_d = uniform(SW_d), LW_d = uniform(LW_d), P_snow = uniform(P_snow))
    FieldExchanger.update_sim!(sim, csf)

    # Wet cells with ice receive the absorbed radiation as a cell mean and the
    # downward snow flux as positive snowfall, each times the cell's
    # `cell_area_ratio`; all other cells receive neither.
    α = FT(0.7) # `Val(:surface_direct_albedo)`
    ϵ = FT(1)   # `Val(:emissivity)`
    cell_values(f) = vec(OC.interior(OC.on_architecture(OC.CPU(), f), :, :, 1))
    ℵ, heat, snow = cell_values(ice_concentration), cell_values(top), cell_values(snowfall)
    r = cell_area_ratio(eg_cpu, grid)
    n_icy = 0
    for c in 1:((Ny - 1) * Nx) # Skip the fold row
        if eg_cpu.oc_wet_area[c] > 0 && ℵ[c] > 0
            n_icy += 1
            @test heat[c] ≈ ℵ[c] * (-(1 - α) * SW_d - ϵ * LW_d) * r[c] rtol = 1e-12
            @test snow[c] ≈ -P_snow * r[c] rtol = 1e-12
        else
            @test heat[c] == 0
            @test snow[c] == 0
        end
    end
    @test n_icy > 0

    # The emission temperature: the ice area fraction times T⁴ integrates to the
    # polygon integral of SIC T⁴ of the polygons' last diagnosed skin temperatures.
    T_sfc_new = 240 .+ 30 .* rand(Random.Xoshiro(3), FT, eg_cpu.n_poly)
    copyto!(remapping.ice_flux_state.T_sfc_new, T_sfc_new)
    T_field = CC.Fields.zeros(boundary_space)
    Interfacer.get_field!(T_field, sim, Val(:surface_temperature))
    CRExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
    Jw = to_cpu(CRExt.se_node_weights(boundary_space))
    coverage = to_cpu(CRExt.se_field_to_vec(remapping.weight_cov_scratch))
    T = to_cpu(CRExt.se_field_to_vec(T_field))
    emitted = sum(eg_cpu.b_area .* ℵ[eg_cpu.oc_of_poly] .* T_sfc_new .^ 4)
    @test sum(Jw .* coverage .* T .^ 4) ≈ emitted rtol = 1e-10
end
