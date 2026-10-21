# Exchange (intersection) grid

When the coupled simulation uses an Oceananigans ocean (the CMIP
configuration), ClimaCoupler can compute atmosphere-ocean and atmosphere-ice
turbulent fluxes on the *exchange grid*: the set of polygons formed where the
elements of the ClimaCore cubed-sphere spectral-element (SE) boundary space
intersect the finite-volume (FV) cells of the ocean's `TripolarGrid`. This is
enabled by default and controlled by the `use_intersection_grid` configuration
option.

## Why

Computing fluxes on the boundary (atmosphere) space and remapping them to the
ocean loses coastline detail: the ocean's land-sea mask (its immersed
bathymetry) does not coincide with the coupler's land fraction, so flux
weights and the region where the ocean actually has wet cells disagree
(issue [#1838](https://github.com/CliMA/ClimaCoupler.jl/issues/1838)). The
exchange grid solves both problems at once, because the *same* intersection
areas define

1. the surface area fractions, and
2. the weights with which per-polygon fluxes are aggregated to both sides.

## Geometry and weights

For a more in-depth derivation of the regridding operations, see [the 
ConservativeRegridding.jl documentation](https://juliageo.org/ConservativeRegridding.jl/stable/extensions/climacore/).
`build_exchange_grid` constructs the polygons once on the CPU (in `Float64`)
via ConservativeRegridding.jl's intersection-operator API, with polygons over
dry (immersed) cells dropped. Each polygon ``\\Omega_k`` in SE element ``e``
carries two sets of weights.

- **Gather** (SE nodes → FV polygons) uses the SEM basis integrals
  ``B_{kn} = \\int_{\\Omega_k} \\phi_n \\, dA``:
  ``\\bar f_k = \\sum_n (B_{kn}/\\sum_n B_{kn}) f_n``. Rows sum to 1, so
  constants are preserved exactly, and the atmospheric state is interpolated
  to high order.
- **Scatter** (FV polygons → SE nodes) uses non-negative *control-volume*
  weights ``\\tilde B_{kn}``: the area of ``\\Omega_k`` inside the control
  volume of node ``n``. The control volumes are boxes in reference
  coordinates that tile the element, close to the GLL subcells (whose 1D
  widths are the GLL weights). Their edges are placed element by element so
  that each box's area under the interpolated Jacobian is exactly
  ``(Jw)_n``. Then
  ``F_n = \\sum_k (\\tilde B_{kn}/(Jw)_n) F_k``, and the weights at a node
  sum to exactly 1 over the polygons tiling its control volume. The scatter
  is followed by `weighted_dss!`.
- **Downward fluxes** (radiation and precipitation, SE nodes → FV polygons)
  use the transpose of the scatter:
  ``\\bar g_k = \\sum_n (\\tilde B_{kn}/\\sum_n \\tilde B_{kn}) g_n``. Each
  polygon then receives a convex combination of the nodal values whose
  control volumes it overlaps, so these fluxes are resolved at the scale of
  the control volumes. The polygons together receive exactly the
  atmosphere's share ``\\sum_n (Jw)_n \\alpha_n g_n`` for an area fraction
  ``\\alpha_n`` built from the scatter weights and DSS'd, as the surface
  fractions are. The atmosphere's radiation and precipitation can
  differ between the element copies of a shared node, so they are DSS'd
  (weighted, which keeps ``\\sum_n (Jw)_n g_n``) before the gather.
- **Nodal coverage**: `node_cov` (sum over wet polygons) is the fraction of
  each node's area over wet ocean, in ``[0, 1]`` by construction. The
  wet-ocean fraction is the DSS'd `node_cov`.

The scatter cannot use ``B_{kn}``, because the GLL Lagrange basis is
sign-indefinite for ``N_q \\ge 3``. Summed over all polygons, ``B_{kn}`` is
close to the nodal weight. Summed over a subset (wet polygons, or
polygons weighted by sea-ice concentration) it can cancel. Then a transfer
normalized by that sum is unbounded, and the sum is not an area fraction. The
control-volume weights have the same nodal masses as ClimaCore's GLL
quadrature, and each normalized scatter is a convex combination of the
polygon values.

On the FV side, each polygon belongs to exactly one cell, and gathers are
direct reads. A flux is scattered as
``F_c = \\sum_{k \\in c} (\\sum_n \\tilde B_{kn}) F_k / A_c``, with ``A_c`` the
FV grid's own cell area, so the FV integral ``\\sum_c A_c F_c`` equals the
SE-side integral exactly even where the two grids disagree on a cell's area.
A state, such as the diagnosed sea-ice surface temperature, is scattered as
the mean weighted by the geometric polygon areas.

## Surface fractions


Each coupling step, `FieldExchanger.align_surface_fractions!` composes

```
ice   = clamp(ice_coverage, 0, wet)
ocean = wet - ice
land  = 1 - wet
```

which sum to 1 identically (asserted at runtime). `ice_coverage` is the DSS'd
control-volume scatter of the per-polygon sea-ice concentration. That is exactly the
denominator of the ice flux scatter; `wet - ice` is the denominator of the
open-ocean one. So the fractions and the fluxes share one set of
weights, and `area_fraction × flux` integrates exactly to
``\\sum_k b_k w_k F_k``, the polygon flux integral over the control-volume
areas ``b_k = \\sum_n \\tilde B_{kn}``.

## Per-polygon fluxes

For the ocean and for sea ice, `FluxCalculator.compute_surface_fluxes!` runs
one `SurfaceFluxes.surface_fluxes` evaluation per polygon (the ice version
includes the skin-temperature diagnosis; ice-free polygons short-circuit).
Atmospheric inputs are gathered from the SE nodes with velocities converted
to the global UV (east/north) basis; surface inputs (SST, sea-ice
concentration, ice thickness, interface temperature) are read directly from
the owning FV cell.

Aggregation:

- **To the ocean/ice**: the sea-ice-concentration weighting is applied *per
  polygon* — `(1 - sic)` for ocean heat/salinity/momentum, per-ice-area
  fluxes for ClimaSeaIce — before the conservative scatter to the cells, so
  the ice/ocean partition is resolved at polygon level rather than after
  remapping.
- **To the atmosphere**: the raw (coverage-weighted) control-volume scatter is DSS'd
  in the UV basis and divided by the DSS'd weighted coverage. Ocean fluxes
  are `(1 - sic)`-weighted averages and ice fluxes `sic`-weighted ones.
  Momentum is converted back to the local CT1/CT2 components only at the end
  (the UV basis is global, so its components are DSS-safe scalars).
- **Surface temperature for radiation**: the ocean and sea-ice temperatures
  from which the atmosphere computes the surface emission are sent up the
  same way using the weighted average of ``T^4``. At each exchange, the emitted 
  longwave the atmosphere computes then integrates to what the ocean and sea ice emit.

Slow surfaces time-average on per-polygon accumulators
(`FluxCalculator.push_and_reset!` overrides), so the averaged flux pushed to
the ocean/ice retains full polygon resolution.

## Fallbacks

The exchange grid requires a process-local `SpectralElementSpace2D` boundary
space; the column (`PointSpace`) setup and `use_intersection_grid: false`
fall back to the regridder-based boundary-space flux path and the legacy
fraction update. Distributed runs take the same fallback, which fails under
MPI; see [Multi-GPU coupled runs](@ref).

## API

```@docs
    FieldExchanger.align_surface_fractions!
```
