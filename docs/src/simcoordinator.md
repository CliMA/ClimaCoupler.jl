# SimCoordinator

This module contains functions for coordinating the execution of coupled simulations,
including stepping through time and running full simulations. The `SimCoordinator` module
provides the main control flow for advancing coupled simulations forward in time.

## Overview

The `SimCoordinator` module provides two key functions for running coupled simulations:

- **`step!`**: Advances the simulation by one coupling timestep, coordinating all component
  models and field exchanges.
- **`run!`**: Executes the full simulation from start to finish, handling precompilation,
  timing, and cleanup.

These functions orchestrate the interaction between component models (atmosphere, land, ocean,
sea ice) through the coupler, ensuring proper field exchanges, flux calculations, and
conservation checks at each timestep.

## Concurrent component stepping

By default `step!` advances the component models one after another. Two
configuration options change that; both are off by default.

`step_concurrently` steps sea ice then ocean in a task of their own while the
calling task steps land then atmosphere. Land and atmosphere are implicitly
coupled, and the sea ice model holds views into the ocean's surface velocity and
salinity, so the two cannot step in parallel.
A component joins the ocean/sea ice group through `Interfacer.is_overlapped`.
This only pays on a GPU, where each component submits to its own CUDA stream; on
a CPU device every component fans out over all threads through
KernelAbstractions, and the coupler warns.

`overlap_slow_surfaces`, which requires `step_concurrently`, runs each ocean/sea
ice step as a task across the `k = dt_ocean / dt_cpl` coupling steps it spans,
instead of within one. The step is launched at the end of the coupling step that
assembles its forcing, and joined after the components step in the coupling step
that needs its result. While it is in flight the coupler holds the ocean and sea
ice state it last read: `update_surface_fractions!`, `update_sim!`, the
accumulator push and `ocean_seaice_fluxes!` stand down for them, their share of
the blended surface fields and turbulent fluxes is reused, and their progress
reports wait for the join. Their turbulent fluxes are computed at the join and
held for the window, while radiative and precipitation fluxes keep accumulating,
so the ocean and sea ice receive the window mean of those.

The cadence for `k = 3`:

```
k = dt_ocean / dt_cpl = 3.  A = atmosphere + land,  Oₙ = the nth ocean/sea ice
step.  Each column is one coupling step; ──▶ marks a step spanning several.

                        │  1  │  2  │  3  │  4  │  5  │  6  │  per window
                        ├─────┼─────┼─────┼─────┼─────┼─────┤
  sequential            │  A  │  A  │ A O₁│  A  │  A  │ A O₂│  3A + O
                        ├─────┼─────┼─────┼─────┼─────┼─────┤
  step_concurrently     │  A  │  A  │  A  │  A  │  A  │  A  │  2A + max(A,O)
                        │     │     │  O₁ │     │     │  O₂ │
                        ├─────┼─────┼─────┼─────┼─────┼─────┤
  overlap_slow_surfaces │  A  │  A  │  A  │  A  │  A  │  A  │  max(3A, O)
                        │     │     │ O₁──────────▶   │ O₂──▶
                        │ atmos uses O₀               │ O₁  │
                        └─────┴─────┴─────┴─────┴─────┴─────┘
```

The first two give the same answers. Under the third the atmosphere sees each
ocean step `k − 1` coupling steps later, or one step later when `k = 1`.

## Choosing between the concurrency options

Write `A` for the cost of one coupling step of atmosphere, land and coupler work,
and `O` for the cost of one ocean and sea ice step. Over one window the
sequential path costs `k·A + O`. `step_concurrently` hides the ocean behind the
coupling step it shares, saving `min(A, O)`; `overlap_slow_surfaces` hides it
behind the whole window, saving `min(k·A, O)`. Overlapping therefore only
improves on `step_concurrently` when `O > A`, and its saving peaks when
`O ≈ k·A`. For `k = 1` it saves no more than `step_concurrently` and adds a
coupling step of lag.

## SimCoordinator API

```@docs
SimCoordinator.run!
SimCoordinator.step!
```
