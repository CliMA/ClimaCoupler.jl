# Checkpointer

This module contains general functions for logging the model states and restarting simulations. The `Checkpointer` uses `ClimaCore.InputOutput` infrastructure, which allows it to handle arbitrarily distributed logging and restart setups.

## Checkpointer API

```@docs
    ClimaCoupler.Checkpointer.get_model_prog_state
    ClimaCoupler.Checkpointer.restart_model_state!
    ClimaCoupler.Checkpointer.checkpoint_model_state
```
