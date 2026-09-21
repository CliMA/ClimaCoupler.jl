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
are not arbitrary. Land and atmosphere are implicitly coupled, and ocean and
sea ice must step together to properly pass the frazil heat flux.
Group membership is decided by type, through
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

`prime_slow_surfaces` shifts that schedule one window earlier, by taking a slow
step during initialization. The cadence of all four cases, for `k = 3`:

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
                        │ O₁──────────▶   │ O₂──────────▶   │
                        │ atmos uses O₀   │ atmos uses O₁   │  one window old
                        ├─────┼─────┼─────┼─────┼─────┼─────┤
  + prime_slow_surfaces │  A  │  A  │  A  │  A  │  A  │  A  │  max(3A, O)
               O₁──────▶│ O₂──────────▶   │ O₃──────────▶   │
                        │ atmos uses O₁   │ atmos uses O₂   │  current
                        └─────┴─────┴─────┴─────┴─────┴─────┘
               ↑ extra step taken during initialization
```

The first two leave the same answers: the join happens inside the coupling step,
so the atmosphere sees ocean state of the same age either way. The third does
not — the atmosphere sees an ocean one window older. The fourth restores the
timing at the cost of integrating the ocean under the previous window's forcing.

How much this buys depends on how many coupling steps one slow step spans. Where
the slow timestep equals `dt_cpl`, so that one slow step spans a single coupling
step, there is no idle atmosphere and land work to hide the ocean behind and the
compute saving is the same as `step_concurrently` alone. The option is still not
idle in that case: `step_concurrently` joins both groups before the coupler's
exchange and flux work, leaving that work — regridding between component grids,
and on several nodes the communication that goes with it — running with nothing
behind it, whereas an overlapped slow step is still in flight across it. So a
single-step window trades a coupling step of lag for overlapped communication
rather than for overlapped computation. Configurations that give the ocean and
sea ice a timestep several coupling steps long get both.

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

## Choosing between the concurrency options

Both options are off by default, and which one is worth enabling depends on the
relative cost of the two groups rather than on preference. Write `A` for the
cost of one coupling step of atmosphere, land and coupler work, `O` for the cost
of one ocean and sea ice step, and `k` for the number of coupling steps one slow
step spans (`dt_ocean / dt_cpl`). Over one window the sequential path costs
`k·A + O`. `step_concurrently` can hide the ocean behind the single coupling
step it shares, so it saves `min(A, O)`. `overlap_slow_surfaces` spreads that
step across the whole window, so it saves `min(k·A, O)`.

Two consequences follow. `overlap_slow_surfaces` only improves on
`step_concurrently` when `O > A` — when the ocean step does not already fit
inside one coupling step of atmosphere work. And the saving peaks when
`O ≈ k·A`, falling away in both directions: a cheap ocean leaves little to hide,
and an ocean costing far more than a window cannot be hidden by it. Ocean
resources are therefore best sized to approach that balance rather than
minimised.

Measured for the standard CMIP configuration (prognostic EDMF, `h_elem` 16,
`Float64`, ORCA ocean) on an A100 with `dt_ocean` at 1800 s, so `k = 5`:
`A = 0.47 s`, `O = 0.27 s`, and hourly radiation adding 3.06 s every ten
coupling steps. That puts `O/(k·A)` at about 0.11 — far below the balance point,
and with `O < A` — so at this resolution `step_concurrently` alone already
captures the whole available saving of roughly 6%, and `overlap_slow_surfaces`
adds only its lag. Radiation, at over a third of runtime, bounds what any
concurrency scheme can achieve here.

That changes as the ocean is refined. Holding `dt_ocean` fixed, `O` grows with
the ocean cell count while `A` barely moves, which carries the ratio into the
range where overlapping is the larger win:

| ocean resolution | `O` | `O/(k·A)` | `step_concurrently` | `overlap_slow_surfaces` |
|:--|--:|--:|--:|--:|
| 1° (measured) | 0.27 s | 0.11 | 6% | 6% |
| 1/4° (measured) | 2.03 s | 0.79 | 21% | 33% |
| 1/12° (extrapolated) | ~18 s | ~7 | 9% | 18% |

The middle row is what the option exists for: overlapping is worth roughly 1.6
times plain concurrency there. The last row is extrapolated linearly from the
1/4° measurement and shows the far side of the balance point, where the run has
become ocean-bound and the remedy is more ocean resource rather than more
overlapping.

Two things are worth noting about how the ratio moves. Refining the ocean
sixteen-fold raised `O` by only about eight, because at 1° the ocean is too
small to keep a modern GPU busy and refinement buys back some of that
efficiency; do not assume `O` scales with cell count. And the coupler's exchange
cost is not resolution-independent — regridding between the atmosphere and a
sixteen-times finer ocean cost close to four times as much here, which is part
of why `A` rises slightly between the rows.

This all assumes the ocean is refined harder than the atmosphere, which is the
usual case since the oceanic Rossby radius is much smaller. Under *proportional*
refinement the atmosphere grows faster, because shortening its timestep for
stability costs an extra factor of resolution, and the ratio moves the other
way.

Neither option is answer-neutral in the strict sense, since the model is not
bitwise reproducible on a GPU, so any comparison has to be read against a noise
floor measured from repeated identical runs — and from *several*, since a single
pair gives one sample of a distribution rather than a threshold. The figures
below come from five identical baselines over three simulated days, giving ten
within-noise pairs to place each option against.

`step_concurrently` is indistinguishable from sequential stepping, as expected:
it changes when the groups run, not what they read.

`overlap_slow_surfaces` separates cleanly from the noise on the ocean and sea
ice surface fractions, at about four times the largest difference the baselines
produce among themselves, with no overlap between the two sets at all. That much
is close to definitional, since `update_surface_fractions!` stands down for the
slow components while a slow step is in flight, so the fractions are held across
the window by construction. It also carries a smaller but systematic signal into
surface temperature and the fluxes tied to it — upward longwave and sensible
heat — at between 1.05 and 1.4 times the noise. The turbulent energy flux stayed
within the noise range.

`prime_slow_surfaces` did not reduce either effect, and was marginally larger on
both. That is consistent with the mechanism rather than surprising: priming
changes *when* the slow group steps, not whether the surface fractions are
refreshed inside a window, so it cannot address the part of the lag that comes
from holding them.

None of this is a stability result. Three days is simply the longest window over
which these differences can be separated from other known behaviour; it says
nothing about longer integrations.

## SimCoordinator API

```@docs
SimCoordinator.run!
SimCoordinator.step!
```
