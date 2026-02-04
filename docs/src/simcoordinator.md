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

## SimCoordinator API

```@docs
SimCoordinator.run!
SimCoordinator.step!
```
