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

`build_exchange_grid` constructs the polygons once on the CPU (in `Float64`)
via ConservativeRegridding.jl's intersection-operator API, with polygons over
dry (immersed) cells dropped. For each polygon ``\\Omega_k`` in SE element
``e``, the SEM basis integrals ``B_{kn} = \\int_{\\Omega_k} \\phi_n \\, dA``
provide:

- a **gather** (SE nodes → polygons): ``\\bar f_k = \\sum_n (B_{kn}/\\sum_n
  B_{kn}) f_n`` — rows sum to 1, constants preserved exactly;
- a **scatter** (polygons → SE nodes): the per-element L2 projection
  ``F_n = \\sum_k (B_{kn}/(Jw)_n) F_k``, followed by `weighted_dss!`;
- the **nodal wet coverage** `node_cov` (the L2 projection of the ocean wet
  mask) and its unmasked counterpart `node_cov_total`; their ratio cancels the
  systematic error of approximating strongly curved FV cells (e.g. along the
  tripolar fold) by four-corner spherical quads.

On the FV side each polygon belongs to exactly one cell; gathers are direct
reads and scatters are area-weighted averages, followed by a mirror of the
tripolar fold's duplicated (shadow) cells. All couplings are stored in CSR
form, so every per-step operation is a race-free segmented reduction with no
atomics — a serial loop on CPU and a KernelAbstractions kernel on GPU, with
zero (CPU) or bounded (GPU) per-step allocations.

## Surface fractions

The wet-ocean fraction on the boundary space is the clamped ratio
`node_cov / node_cov_total`, made C0-continuous with `weighted_dss!` and then
**smoothed with the same diffusion recipe ClimaAtmos applies to its Earth
orography** (`ClimaCore.Hypsography.diffuse_surface_elevation!` with
``\\kappa = 0.05\\,\\Delta h^2``, ``\\mathrm{maxiter} =
\\mathrm{round}(\\log(d)/0.05)`` where ``d`` is `topography_damping_factor`).
The coupler therefore never sees land-sea contrasts sharper than the
atmosphere's own smoothed topography. The `topography_damping_factor`
configuration option must match the ClimaAtmos option of the same name
(default 5).

Each coupling step, `FieldExchanger.align_surface_fractions!` composes

```
ice   = clamp(ice_concentration, 0, wet)
ocean = wet - ice
land  = 1 - wet
```

which sum to 1 identically and supersede the legacy ETOPO-based update.

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
  fluxes for ClimaSeaIce — before the conservative area-weighted scatter, so
  the ice/ocean partition is resolved at polygon level rather than after
  remapping.
- **To the atmosphere**: the raw (coverage-weighted) L2 scatter is DSS'd in
  the UV basis and divided by the DSS'd weighted coverage — ocean fluxes are
  `(1 - sic)`-weighted and ice fluxes `sic`-weighted averages, consistent
  with the open-ocean and ice area fractions `update_flux_fields!` multiplies
  them by. Momentum is converted back to the local CT1/CT2 components only at
  the end (the UV basis is global, so its components are DSS-safe scalars).

Slow surfaces time-average on per-polygon accumulators
(`FluxCalculator.push_and_reset!` overrides), so the averaged flux pushed to
the ocean/ice retains full polygon resolution.

## Fallbacks

The exchange grid requires a process-local `SpectralElementSpace2D` boundary
space; column (`PointSpace`) and distributed setups automatically fall back
to the regridder-based boundary-space flux path and the legacy fraction
update, as does `use_intersection_grid: false`.

## API

```@docs
    FieldExchanger.align_surface_fractions!
```
