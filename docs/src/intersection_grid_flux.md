# Intersection Grid Flux Exchange

This page demonstrates how flux exchange works on the intersection grid between
the ClimaCore cubed-sphere (CC) atmosphere grid and the Oceananigans
LatitudeLongitudeGrid (OC) ocean/sea-ice grid.

## Motivation: Why Intersection-Grid Fluxes?

Traditional flux exchange computes fluxes on one grid (typically the atmosphere's
boundary space) and then remaps the results to other grids. This causes:

1. **Gradient smoothing at coastlines**: Sharp temperature gradients between land
   and ocean get averaged when surface properties are remapped before flux computation.

2. **Area fraction approximation**: Partial cells use area fractions that approximate
   the true intersection geometry.

The intersection-grid approach computes fluxes directly on the polygons formed by
overlapping CC elements and OC cells, preserving sharp gradients and using exact areas.

```@example intersection_flux
import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore as CC
import Oceananigans as OC
import ClimaOcean, ClimaSeaIce, KernelAbstractions, ConservativeRegridding, Adapt
import ClimaCoupler
import ClimaCoupler: Interfacer

CMIPExt = Base.get_extension(ClimaCoupler, :ClimaCouplerCMIPExt)
@assert !isnothing(CMIPExt)

FT = Float64
context = ClimaComms.context()
ClimaComms.init(context)
arch = OC.CPU()

# Create a small cubed-sphere boundary space (atmosphere grid)
boundary_space = CC.CommonSpaces.CubedSphereSpace(
    FT;
    radius = FT(6.371e6),
    n_quad_points = 4,
    h_elem = 16,  # 6 elements per cube face edge
    context,
)

# Create Oceananigans lat-lon grid (ocean grid)
Nx, Ny, Nz = 360, 180, 1
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
    OC.GridFittedBottom((x, y) -> -50.0);
    active_cells_map = false,
)

# Construct remapper and extract intersection grid
remapping = CMIPExt.construct_remapper(grid, boundary_space)
ig = remapping.intersection_grid

println("Grid Statistics:")
println("  CC elements: ", ig.n_cc)
println("  OC cells: ", ig.n_oc)
println("  Intersection polygons: ", ig.n_intersections)
println("  Polygons per CC element (avg): ", round(ig.n_intersections / ig.n_cc, digits=1))
println("  Polygons per OC cell (avg): ", round(ig.n_intersections / ig.n_oc, digits=1))
```

## Visualizing the Grids

```@example intersection_flux
using CairoMakie

# Get CC element coordinates (centroids)
coords = CC.Fields.coordinate_field(boundary_space)
lats_cc = parent(coords.lat)
lons_cc = parent(coords.long)

# Get OC cell center coordinates
lons_oc = collect(OC.λnodes(grid, OC.Center(), OC.Center(), OC.Center()))
lats_oc = collect(OC.φnodes(grid, OC.Center(), OC.Center(), OC.Center()))

fig = Figure(size = (1000, 500))

ax1 = Axis(fig[1, 1];
    title = "CC Elements (Cubed-Sphere)",
    xlabel = "Longitude (°)",
    ylabel = "Latitude (°)",
)
scatter!(ax1, vec(lons_cc), vec(lats_cc); 
    markersize = 3, 
    color = :blue, 
    alpha = 0.6,
    label = "CC element centers"
)
xlims!(ax1, -180, 180)
ylims!(ax1, -90, 90)

ax2 = Axis(fig[1, 2];
    title = "OC Cells (LatLon Grid)",
    xlabel = "Longitude (°)",
    ylabel = "Latitude (°)",
)
# Create grid of OC cell centers
oc_lon_grid = repeat(lons_oc, 1, length(lats_oc))
oc_lat_grid = repeat(lats_oc', length(lons_oc), 1)
scatter!(ax2, vec(oc_lon_grid), vec(oc_lat_grid);
    markersize = 4,
    color = :red,
    alpha = 0.6,
    label = "OC cell centers"
)
xlims!(ax2, -180, 180)
ylims!(ax2, -90, 90)

fig
```

## Visualizing the Intersection Polygons

So far we've only looked at where CC elements and OC cells *live*. We can also
draw the intersection polygons themselves — one for each row of `(cc_indices,
oc_indices, areas)` in `remapping.intersection_grid` — by re-using the
[`ConservativeRegridding`](https://github.com/CliMA/ConservativeRegridding.jl)
trees that the regridder builds internally, then clipping each pair with
`GeometryOps`.

The helpers below reconstruct the polygons in `(lon°, lat°)` coordinates and
do the two pieces of map-plotting work that flat (lon, lat) rendering needs:

1. **Pole-vertex resolution.** Cube-face elements that touch a pole have one
   vertex with `lat = ±90°`, where longitude is degenerate.  `atand` returns
   an arbitrary value (e.g. `-180°`), which makes the edge from a non-pole
   neighbour to the pole vertex appear to span up to 360° of longitude.  We
   replace each pole vertex with up to two vertices at `lat = ±90°` using the
   longitudes of its non-pole neighbours, so the polar cap is drawn as a
   proper meridional wedge.
2. **Antimeridian cutting.** Polygons that genuinely cross `lon = ±180°` are
   first *unwrapped* into a single continuous (but possibly out-of-range) ring
   by shifting successive vertices by `±360°`, then cut at each meridian
   `±180°, ±540°, …` with `GeometryOps.cut`.  Each resulting piece is shifted
   back into `[-180°, 180°]`.

```@example intersection_flux
import GeometryOps as GO
import GeoInterface as GI
import GeoMakie

# --- low-level conversions -------------------------------------------------
function _xyz_of(pt)
    v = hasproperty(pt, :v) ? pt.v : pt
    return (v[1], v[2], v[3])
end
function _usp_to_lonlat(pt)
    x, y, z = _xyz_of(pt)
    return (atand(y, x), asind(clamp(z, -1.0, 1.0)))
end
_point_to_lonlat(pt::Tuple) = (Float64(pt[1]), Float64(pt[2]))
_point_to_lonlat(pt)        = _usp_to_lonlat(pt)
function _polygon_lonlat_ring(poly)
    pts = collect(GI.getpoint(GI.getexterior(poly)))
    return [_point_to_lonlat(p) for p in pts]
end

const POLE_EPS = 0.05  # treat |lat| > 90 - POLE_EPS as "at the pole"

# Drop the closing duplicate from a ring if present.
function _open_ring(ring)
    return (length(ring) >= 2 && ring[1] == ring[end]) ? ring[1:end-1] : collect(ring)
end

# Replace each pole vertex with up-to-two vertices at lat = ±90° whose
# longitudes match the nearest non-pole neighbours.
function _resolve_pole_vertices(pts)
    n = length(pts)
    n == 0 && return pts
    is_pole = i -> abs(pts[i][2]) > 90 - POLE_EPS
    any_non_pole = findfirst(i -> !is_pole(i), 1:n)
    out = Tuple{Float64,Float64}[]
    for i in 1:n
        cur = pts[i]
        if !is_pole(i)
            push!(out, cur); continue
        end
        prev_lon = nothing; j = i
        for _ in 1:n
            j = mod1(j - 1, n)
            !is_pole(j) && (prev_lon = pts[j][1]; break)
        end
        next_lon = nothing; j = i
        for _ in 1:n
            j = mod1(j + 1, n)
            !is_pole(j) && (next_lon = pts[j][1]; break)
        end
        pole_lat = sign(cur[2]) * 90.0
        if isnothing(prev_lon) || isnothing(next_lon)
            push!(out, ((any_non_pole === nothing ? cur[1] : pts[any_non_pole][1]), pole_lat))
        elseif prev_lon == next_lon
            push!(out, (prev_lon, pole_lat))
        else
            push!(out, (prev_lon, pole_lat))
            push!(out, (next_lon, pole_lat))
        end
    end
    return out
end

# Make consecutive lons differ by no more than 180° by shifting by ±360°.
function _unwrap_lons(pts)
    n = length(pts)
    n == 0 && return pts
    out = [(Float64(pts[1][1]), Float64(pts[1][2]))]
    prev = out[1][1]
    for i in 2:n
        l, la = Float64(pts[i][1]), Float64(pts[i][2])
        while l - prev > 180; l -= 360; end
        while l - prev < -180; l += 360; end
        push!(out, (l, la))
        prev = l
    end
    return out
end

# Cut the unwrapped polygon at meridians ±180°, ±540°, … and shift each
# resulting piece back into [-180°, 180°].
function _cut_at_antimeridians(pts)
    pts_closed = vcat(pts, [pts[1]])
    poly = GI.Polygon([GI.LinearRing(pts_closed)])
    lons = [p[1] for p in pts]
    lo, hi = extrema(lons)

    cut_lons = Float64[]
    for m in (-540.0, -180.0, 180.0, 540.0, 900.0)
        lo < m < hi && push!(cut_lons, m)
    end
    sort!(cut_lons)

    pieces = Any[poly]
    for m in cut_lons
        new_pieces = Any[]
        line = GI.Line([(m, -90.0), (m, 90.0)])
        for p in pieces
            res = try
                GO.cut(p, line)
            catch
                Any[p]
            end
            append!(new_pieces, res)
        end
        pieces = new_pieces
    end

    rings = Vector{Vector{Tuple{Float64,Float64}}}()
    for p in pieces
        ring = [(Float64(GI.x(pt)), Float64(GI.y(pt))) for pt in GI.getpoint(GI.getexterior(p))]
        ring_open = (ring[1] == ring[end]) ? ring[1:end-1] : ring
        ls = [r[1] for r in ring_open]
        centre = (minimum(ls) + maximum(ls)) / 2
        shift = if centre < -180
            360.0 * cld(-180 - centre, 360)
        elseif centre > 180
            -360.0 * cld(centre - 180, 360)
        else
            0.0
        end
        push!(rings, [(r[1] + shift, r[2]) for r in ring_open])
    end
    return rings
end

# End-to-end: spherical polygon -> list of map-safe (lon, lat) rings.
function _process_ring(ring)
    pts = _open_ring(ring)
    length(pts) < 3 && return Vector{Vector{Tuple{Float64,Float64}}}()
    pts = _resolve_pole_vertices(pts)
    pts = _unwrap_lons(pts)
    return _cut_at_antimeridians(pts)
end

function intersection_polygons_lonlat(remapping, boundary_space, grid)
    ig = remapping.intersection_grid

    boundary_space_cpu = CC.Adapt.adapt(Array, boundary_space)
    grid_underlying_cpu = OC.on_architecture(OC.CPU(), grid.underlying_grid)

    FT_cc = CC.Spaces.undertype(boundary_space_cpu)
    R = CC.Spaces.topology(boundary_space_cpu).mesh.domain.radius
    manifold = GO.Spherical(; radius = FT_cc(R))

    dst_tree = ConservativeRegridding.Trees.treeify(manifold, boundary_space_cpu)
    src_tree = ConservativeRegridding.Trees.treeify(manifold, grid_underlying_cpu)

    rings = Vector{Vector{Tuple{Float64,Float64}}}()
    areas = Float64[]
    for k in 1:ig.n_intersections
        cc_poly = ConservativeRegridding.Trees.getcell(dst_tree, ig.cc_indices[k])
        oc_poly = ConservativeRegridding.Trees.getcell(src_tree, ig.oc_indices[k])
        ipolys_raw = GO.intersection(
            GO.ConvexConvexSutherlandHodgman(manifold),
            cc_poly, oc_poly;
            target = GI.PolygonTrait(),
        )
        # `GO.intersection` returns either a single Polygon or a Vector of them.
        ipolys = ipolys_raw isa AbstractVector ? ipolys_raw : (ipolys_raw,)
        for poly in ipolys
            for branch in _process_ring(_polygon_lonlat_ring(poly))
                length(branch) >= 3 || continue
                push!(rings, branch)
                push!(areas, ig.areas[k])
            end
        end
    end
    return rings, areas, dst_tree
end

rings, poly_areas, dst_tree = intersection_polygons_lonlat(remapping, boundary_space, grid)
println("Reconstructed $(length(rings)) map-safe intersection polygons.")
```

We can now overlay the intersection polygons (filled, colored by area) on a
GeoMakie projection, with the parent CC elements drawn as black outlines so
you can see how each cubed-sphere element is carved into pieces by the
LatitudeLongitudeGrid:

```@example intersection_flux
# CC element outlines, processed the same way as the intersection polygons.
cc_outline_rings = Vector{Vector{Tuple{Float64,Float64}}}()
for idx in 1:ig.n_cc
    poly = ConservativeRegridding.Trees.getcell(dst_tree, idx)
    for branch in _process_ring(_polygon_lonlat_ring(poly))
        length(branch) >= 3 && push!(cc_outline_rings, branch)
    end
end

fig_poly = Figure(size = (1300, 700))
ax_poly = GeoMakie.GeoAxis(fig_poly[1, 1];
    title  = "Intersection polygons (n = $(length(rings))), colored by area [m²]",
    dest   = "+proj=eqearth",
)
GeoMakie.lines!(ax_poly, GeoMakie.coastlines();
    color = :black, linewidth = 0.5,
)
pl = poly!(ax_poly,
    [Point2f.(r) for r in rings];
    color       = poly_areas,
    colormap    = :viridis,
    strokewidth = 0.2,
    strokecolor = (:black, 0.5),
)
poly!(ax_poly,
    [Point2f.(r) for r in cc_outline_rings];
    color       = (:white, 0.0),
    strokewidth = 0.8,
    strokecolor = (:black, 0.8),
)
Colorbar(fig_poly[1, 2], pl; label = "intersection area [m²]")
fig_poly
```

Each colored polygon is one row of `remapping.intersection_grid`; together they
tile the globe.  Polygons near the poles look larger in the equal-earth
projection because the spherical-cap geometry — but the colorbar shows that
their actual area is similar to (or smaller than) the tropical polygons,
which is exactly what the area-weighted gather/scatter operations rely on.

A standalone script with the same plotting helper is available at
`experiments/test/plot_intersection_polygons.jl`. From a REPL that has an
`OceananigansSimulation` in scope you can simply do

```julia
include("experiments/test/plot_intersection_polygons.jl")
fig = plot_intersection_polygons(ocean_sim;
    output_file = "intersection_polygons.png",
    color_by    = :area,        # or :cc_index / :oc_index / :random
    draw_cc     = true,
)
```

## Continents on the Intersection Plot

In a CMIP coupled simulation, the Oceananigans grid is an
`ImmersedBoundaryGrid` whose bathymetry mask defines where the ocean model
actually solves equations.  Continents show up as "dry" (immersed) ocean cells
that get skipped during the time step.  The intersection grid still contains
polygons over those dry cells and represent CC-OC geometric overlaps — so
we colour each intersection polygon by whether its parent OC cell is "wet"
or "dry", and immediately see where the ocean model is and isn't active.

To make this self-contained inside the doc we rebuild the OC grid with
[`ClimaOcean.regrid_bathymetry`](https://clima.github.io/ClimaOcean.jl/dev/) —
the same call `OceananigansSimulation` uses internally to interpolate ETOPO
onto its lat-lon grid — and then read the wet/dry status of each surface cell
directly from the resulting `ImmersedBoundaryGrid`** via
`Oceananigans.ImmersedBoundaries.immersed_cell`.

```@example intersection_flux
bottom_height = ClimaOcean.regrid_bathymetry(
    underlying_grid;
    minimum_depth        = 30,
    interpolation_passes = 5,
    major_basins         = 1,
)
continent_grid = OC.ImmersedBoundaryGrid(
    underlying_grid,
    OC.GridFittedBottom(bottom_height);
    active_cells_map = false,
)
continent_remapping = CMIPExt.construct_remapper(continent_grid, boundary_space)
continent_ig = continent_remapping.intersection_grid

# Per-OC-cell wet/dry status, read straight from the ocean model's grid.
# `immersed_cell(i, j, k, grid)` returns `true` for cells the model skips.
# Flat index matches the regridder's column-major order, `(j - 1) * Nx + i`.
function ocean_grid_dry_mask(grid, Nx, Ny)
    is_dry = falses(Nx * Ny)
    for j in 1:Ny, i in 1:Nx
        is_dry[(j - 1) * Nx + i] =
            OC.ImmersedBoundaries.immersed_cell(i, j, 1, grid)
    end
    return is_dry
end
is_oc_cell_dry = ocean_grid_dry_mask(continent_grid, Nx, Ny)
n_land_cells  = count(is_oc_cell_dry)
n_ocean_cells = Nx * Ny - n_land_cells
println("OC grid: $n_ocean_cells wet cells, $n_land_cells dry (land) cells.")
```

Now reconstruct the intersection polygons against the new immersed grid, and
tag each polygon's parent OC cell as wet or dry:

```@example intersection_flux
src_tree_cont = ConservativeRegridding.Trees.treeify(
    GO.Spherical(; radius = CC.Spaces.undertype(CC.Adapt.adapt(Array, boundary_space))(
        CC.Spaces.topology(boundary_space).mesh.domain.radius,
    )),
    OC.on_architecture(OC.CPU(), continent_grid.underlying_grid),
)

land_rings  = Vector{Vector{Tuple{Float64,Float64}}}()
ocean_rings = Vector{Vector{Tuple{Float64,Float64}}}()
ocean_areas = Float64[]
for k in 1:continent_ig.n_intersections
    cc_idx = continent_ig.cc_indices[k]
    oc_idx = continent_ig.oc_indices[k]
    cc_poly = ConservativeRegridding.Trees.getcell(dst_tree, cc_idx)
    oc_poly = ConservativeRegridding.Trees.getcell(src_tree_cont, oc_idx)
    ipolys_raw = GO.intersection(
        GO.ConvexConvexSutherlandHodgman(
            GO.Spherical(; radius = CC.Spaces.topology(boundary_space).mesh.domain.radius),
        ),
        cc_poly, oc_poly; target = GI.PolygonTrait(),
    )
    ipolys = ipolys_raw isa AbstractVector ? ipolys_raw : (ipolys_raw,)
    for poly in ipolys
        for branch in _process_ring(_polygon_lonlat_ring(poly))
            length(branch) >= 3 || continue
            if is_oc_cell_dry[oc_idx]
                push!(land_rings, branch)
            else
                push!(ocean_rings, branch)
                push!(ocean_areas, continent_ig.areas[k])
            end
        end
    end
end
println("Intersection polygons: $(length(ocean_rings)) over ocean, ",
        "$(length(land_rings)) over land (no ocean computation).")
```

Plot them together — land polygons in tan (no ocean computation), ocean
polygons coloured by intersection area, plus the CC element outlines for
context:

```@example intersection_flux
fig_land = Figure(size = (1300, 700))
ax_land = GeoMakie.GeoAxis(fig_land[1, 1];
    title = "Intersection polygons over ImmersedBoundaryGrid: " *
            "$(length(ocean_rings)) ocean, $(length(land_rings)) land",
    dest = "+proj=eqearth",
)
# Land first, so ocean polygons stroke over the boundary.
poly!(ax_land,
    [Point2f.(r) for r in land_rings];
    color       = (:tan, 0.85),
    strokewidth = 0.15,
    strokecolor = (:saddlebrown, 0.6),
    label       = "land (no ocean equations solved)",
)
pl_ocean = poly!(ax_land,
    [Point2f.(r) for r in ocean_rings];
    color       = ocean_areas,
    colormap    = :deep,
    strokewidth = 0.15,
    strokecolor = (:black, 0.4),
    label       = "ocean polygon (area-coloured)",
)
GeoMakie.lines!(ax_land, GeoMakie.coastlines();
    color = (:black, 0.45), linewidth = 0.5,
)
poly!(ax_land,
    [Point2f.(r) for r in cc_outline_rings];
    color       = (:white, 0.0),
    strokewidth = 0.8,
    strokecolor = (:black, 0.7),
)
Colorbar(fig_land[1, 2], pl_ocean; label = "ocean intersection area [m²]")
fig_land
```

A few things worth noting from this plot:

* The continents are taken from the ETOPO2022 bathymetry that
  `ClimaOcean.regrid_bathymetry` interpolates onto the OC grid, then encoded
  as the immersed boundary of an `ImmersedBoundaryGrid`.  Coastlines follow
  the OC cell edges at the resolution of the underlying lat-lon grid — so at
  `Nx × Ny = 36 × 20`, large continents like Eurasia and South America are
  well resolved while small islands and narrow channels are merged into the
  surrounding land or ocean. Increasing `Nx`, `Ny` would give finer
  coastlines.
* The wet/dry classification is read **straight from the ocean model's
  grid** via `OC.ImmersedBoundaries.immersed_cell(i, j, k, grid)`, so what
  you see is exactly the set of cells Oceananigans skips on each time step
  — no separate analytical landmask is computed.
* A single cubed-sphere CC element can contain both ocean and land
  intersection polygons.  Those *coastal CC elements* are exactly the ones
  that benefit most from intersection-grid flux exchange: instead of running
  a single SurfaceFluxes calculation with an averaged surface state, each
  intersection polygon uses the surface state from its actual parent OC
  cell — land flux states from the land CC neighbour, ocean state from the
  underlying wet OC cell.
* Land polygons have area in the IntersectionGrid but contribute zero ocean
  fluxes; in a real coupled run their contribution is replaced by the land
  model's surface flux through the area-weighted scatter step.

## The Flux Computation Flow

### Step 1: Gather State to Intersection Grid

For each intersection polygon, we look up:
- **Atmosphere state** from its parent CC element (T, q, ρ, u, v, height)
- **Surface state** from its parent OC cell (T_sfc, roughness)

```@example intersection_flux
# Demonstrate gather operation with a test field

# Create a CC field with latitude-dependent values (simulating temperature)
cc_test = CC.Fields.zeros(boundary_space)
cc_test .= 280.0 .+ 20.0 .* cos.(deg2rad.(coords.lat))

# Extract per-element values
CRExt = CMIPExt.get_ConservativeRegriddingCCExt()
cc_per_element = zeros(FT, ig.n_cc)
field_ones = CC.Fields.ones(boundary_space)
CRExt.get_value_per_element!(cc_per_element, cc_test, field_ones)

# Gather to intersection grid
intersection_values = zeros(FT, ig.n_intersections)
CMIPExt.gather_cc_to_intersection!(intersection_values, ig, cc_per_element)

println("CC field statistics:")
println("  Range: ", extrema(cc_per_element))
println("Intersection values (gathered from CC):")
println("  Range: ", extrema(intersection_values))
println("  Same range confirms correct gather!")
```

### Step 2: Compute Fluxes on Intersection Polygons

Each intersection polygon gets its own flux calculation using the atmosphere
state from its CC element and surface state from its OC cell:

```@example intersection_flux
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import SurfaceFluxes.UniversalFunctions as UF
import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
import ClimaParams as CP

# Create parameters
toml_dict = CP.create_toml_dict(FT)
thermo_params = TDP.ThermodynamicsParameters(toml_dict)
surface_fluxes_params = SFP.SurfaceFluxesParameters(toml_dict, UF.BusingerParams)

# Create atmosphere state (per CC element) - cooler atmosphere
cc_atmos_state = (
    T = fill(FT(285), ig.n_cc),      # 285 K atmosphere
    q_tot = fill(FT(0.008), ig.n_cc), # 8 g/kg humidity
    q_liq = fill(FT(0), ig.n_cc),
    q_ice = fill(FT(0), ig.n_cc),
    ρ = fill(FT(1.2), ig.n_cc),
    u = fill(FT(5), ig.n_cc),         # 5 m/s wind
    v = fill(FT(2), ig.n_cc),
    h = fill(FT(10), ig.n_cc),        # 10 m reference height
)

# Create surface state (per OC cell) - latitude-varying SST
# Warm tropics, cold high latitudes
sst_by_cell = zeros(FT, ig.n_oc)
for j in 1:Ny
    for i in 1:Nx
        cell_idx = (j-1) * Nx + i
        lat = lats_oc[j]
        # SST pattern: warm tropics (303K), cold poles (273K)
        sst_by_cell[cell_idx] = 288.0 + 15.0 * cos(deg2rad(lat))
    end
end

oc_surface_state = (
    T = sst_by_cell,
    z0m = fill(FT(1e-4), ig.n_oc),
    z0b = fill(FT(1e-4), ig.n_oc),
    h = fill(FT(0), ig.n_oc),
)

# Allocate flux state
flux_state = CMIPExt.IntersectionFluxState(FT, ig.n_intersections)

# Compute fluxes on intersection grid
CMIPExt.compute_surface_fluxes_on_intersection!(
    flux_state,
    ig,
    cc_atmos_state,
    oc_surface_state,
    surface_fluxes_params,
    thermo_params,
)

println("\nFlux Statistics on Intersection Grid:")
println("  Sensible heat: min=$(round(minimum(flux_state.flux_sh), digits=1)), max=$(round(maximum(flux_state.flux_sh), digits=1)) W/m²")
println("  Latent heat: min=$(round(minimum(flux_state.flux_lh), digits=1)), max=$(round(maximum(flux_state.flux_lh), digits=1)) W/m²")
```

### Step 3: Scatter Fluxes Back to Component Grids

After computing fluxes on the intersection grid, we aggregate them back:
- To **CC elements**: area-weighted mean for the atmosphere
- To **OC cells**: area-weighted mean for the ocean

```@example intersection_flux
# Scatter sensible heat flux to both grids
cc_flux_sh = zeros(FT, ig.n_cc)
oc_flux_sh = zeros(FT, ig.n_oc)

CMIPExt.scatter_to_cc!(cc_flux_sh, ig, flux_state.flux_sh)
CMIPExt.scatter_to_oc!(oc_flux_sh, ig, flux_state.flux_sh)

println("\nScattered Flux Statistics:")
println("  CC grid SH flux: min=$(round(minimum(filter(x->x!=0, cc_flux_sh)), digits=1)), max=$(round(maximum(cc_flux_sh), digits=1)) W/m²")
println("  OC grid SH flux: min=$(round(minimum(oc_flux_sh), digits=1)), max=$(round(maximum(oc_flux_sh), digits=1)) W/m²")
```

## Visualizing the Flux Pattern

```@example intersection_flux

sst_matrix     = reshape(sst_by_cell, Nx, Ny)
oc_flux_matrix = reshape(oc_flux_sh,  Nx, Ny)

fig2 = Figure(size = (1200, 500))

ax1 = Axis(fig2[1, 1];
    title = "Sea Surface Temperature (Input)",
    xlabel = "Longitude (°)",
    ylabel = "Latitude (°)",
)
hm1 = heatmap!(ax1, lons_oc, lats_oc, sst_matrix;
    colormap = :thermal,
    colorrange = (273, 303),
)
Colorbar(fig2[1, 2], hm1; label = "SST (K)")

ax2 = Axis(fig2[1, 3];
    title = "Sensible Heat Flux (Output)",
    xlabel = "Longitude (°)",
    ylabel = "Latitude (°)",
)
hm2 = heatmap!(ax2, lons_oc, lats_oc, oc_flux_matrix;
    colormap = :RdBu,
    colorrange = (-100, 100),
)
Colorbar(fig2[1, 4], hm2; label = "SH Flux (W/m²)")

fig2
```

## Conservation Properties

The intersection-grid approach preserves area-weighted integrals:

```@example intersection_flux
# Total flux on intersection grid
total_intersection = sum(flux_state.flux_sh .* ig.areas)

# Total flux on OC grid (using OC cell areas)
total_oc = sum(oc_flux_sh .* ig.oc_areas)

# Total flux on CC grid (using intersection areas per CC element)
cc_intersection_area = zeros(FT, ig.n_cc)
for k in 1:ig.n_intersections
    cc_intersection_area[ig.cc_indices[k]] += ig.areas[k]
end
total_cc = sum(cc_flux_sh .* cc_intersection_area)

println("\nConservation Check (total sensible heat flux):")
println("  Intersection grid: $(round(total_intersection / 1e12, digits=3)) TW")
println("  OC grid:           $(round(total_oc / 1e12, digits=3)) TW")
println("  CC grid:           $(round(total_cc / 1e12, digits=3)) TW")
println("  OC relative error: $(abs(total_oc - total_intersection) / abs(total_intersection))")
```

## Comparison: Traditional vs Intersection-Grid Approach

### Traditional Approach
1. Remap surface T from OC → CC (smoothing occurs)
2. Compute fluxes on CC grid
3. Remap fluxes from CC → OC (more smoothing)

### Intersection-Grid Approach  
1. Gather CC atmos state → intersection (no smoothing, just lookup)
2. Gather OC surface state → intersection (no smoothing, just lookup)
3. Compute fluxes on intersection grid (preserves gradients)
4. Scatter fluxes → CC and OC (area-weighted, conservative)

The key difference is that **no field is remapped before flux computation**.
Each intersection polygon uses the exact values from its parent cells, preserving
sharp gradients at coastlines and grid boundaries.

## Performance Considerations

The intersection grid typically has more elements than either source grid:

```@example intersection_flux
println("Element counts:")
println("  CC elements: ", ig.n_cc)
println("  OC cells: ", ig.n_oc)
println("  Intersection polygons: ", ig.n_intersections)
println("  Ratio (intersections / max(CC, OC)): ", round(ig.n_intersections / max(ig.n_cc, ig.n_oc), digits=2))
```

The overhead of computing fluxes on more points is offset by:
1. Eliminating two remapping operations (surface state to CC, fluxes to OC)
2. Better physical fidelity at coastlines
3. Exact area weighting instead of approximate area fractions

For GPU execution, the gather operations are embarrassingly parallel. The scatter
operations require atomic additions or segmented reductions (see implementation notes).
