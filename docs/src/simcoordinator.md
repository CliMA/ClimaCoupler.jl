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
configuration options change that.

`step_concurrently` splits them into two groups that advance at the same time:
land then atmosphere in one task, sea ice then ocean in another. The groupings
are not arbitrary. Land and atmosphere are implicitly coupled, so they have to
run in sequence. Sea ice and ocean have to run in sequence for a different
reason: the sea ice model is constructed holding views into the ocean's surface
velocity and salinity fields, so stepping them in parallel would let the ice
read state the ocean is writing. Group membership is decided by type, through
`Interfacer.is_overlapped`.

This is only worth enabling on a GPU, where each component occupies a single
Julia thread and submits to its own CUDA stream. On a CPU device every component
fans out over all available threads through KernelAbstractions, so the groups
oversubscribe each other and nothing is gained; the coupler warns if it is asked
to do this.

`overlap_slow_surfaces` goes further. Because the ocean and sea ice typically
take a timestep several coupling steps long, `step_concurrently` can only hide
one coupling step of atmosphere and land work behind them. This option instead
launches their step as a task at the coupling step where their forcing is fully
assembled and joins it several steps later, so the whole window overlaps. While
that step is in flight the coupler must not read or write ice or ocean state, so
`step_model_sims!`, `update_surface_fractions!`, `exchange!`,
`update_model_sims!`, the turbulent flux accumulator push and
`ocean_seaice_fluxes!` all stand down for the slow components. The atmosphere
still needs a blended surface temperature and albedo every coupling step, so
`combine_surfaces!` sums the fast and slow surfaces separately and reuses the
slow sum for the rest of the window — land keeps contributing live values.

The cost is that the atmosphere sees an ocean one slow step older than it
otherwise would. `prime_slow_surfaces` removes that by advancing the slow group
one step during initialization, so each overlapped step integrates the window
that is about to happen rather than the one that just did and lands in time to
be used on schedule. The lag does not vanish; it moves to the ocean's forcing,
which is integrated from the previous window's accumulated fluxes. That is the
side better able to absorb it, since a component whose timestep spans several
coupling steps already integrates under forcing held constant across a whole
window.

Under priming the slow components' own diagnostics are written on their own
clocks, so those files carry times that lead coupler time by up to one slow
step. Coupler diagnostics remain on coupler time and hold the surface state the
atmosphere actually saw.

## SimCoordinator API

```@docs
SimCoordinator.run!
SimCoordinator.step!
```
