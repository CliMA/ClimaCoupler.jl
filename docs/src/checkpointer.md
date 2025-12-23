# Checkpointer

## How to save and restart from checkpoints

`ClimaCoupler` supports saving and reading simulation checkpoints. This is
useful to split a long simulation into smaller, more manageable chunks.

Checkpoints are a mix of HDF5 and JLD2 files and are typically saved in a
`checkpoints` folder in the simulation output. See
[`Utilities.setup_output_dirs`](@ref) for more information.

!!! warning "Known limitations"

    - The number of MPI processes has to remain the same across checkpoints
    - Restart files are generally not portable across machines, julia versions,
      and package versions
    - Adding/changing new component models will probably require adding/changing code

### Saving checkpoints

If you are running a model (such as AMIP), chances are that you can enable
checkpointing just by setting a command-line argument; The `checkpoint_dt`
option controls how frequently a checkpoint should be produced.

If your model does not come with this option already, you can checkpoint the
simulation by adding a callback that calls the
[`Checkpointer.checkpoint_sims`](@ref) function.

For example, to add a callback to checkpoint every hour of simulated time,
assuming you have a `start_date`
```julia
import Dates

import ClimaCoupler: Checkpointer, TimeManager
import ClimaDiagnostics.Schedules: EveryCalendarDtSchedule

schedule = EveryCalendarDtSchedule(Dates.Hour(1); start_date)
checkpoint_callback = TimeManager.Callback(schedule_checkpoint, Checkpointer.checkpoint_sims)

# In the coupling loop:
TimeManager.maybe_trigger_callback(checkpoint_callback, coupled_simulation, time)
```

### Reading checkpoints

There are two ways to restart a simulation from checkpoints. By default,
`ClimaCoupler` tries finding suitable checkpoints and automatically use them.
Alternatively, you can specify a directory `restart_dir` and a simulation time
`restart_t` and restart from files saved in the given directory at the given
time. If the model you are running supports writing checkpoints via command-line
argument, it will probably also support reading them. In this case, the
arguments `restart_dir` and `restart_t` identify the path of the top level
directory containing all the checkpoint files and the simulated times in second.

If the model does not support directly reading a checkpoint, the `Checkpointer`
module provides a straightforward way to add this feature.
[`Checkpointer.restart!`](@ref) takes a coupled simulation, a `restart_dir`, and
a `restart_t` and overwrites the content of the coupled simulation with what is
in the checkpoint.

## Developer notes

In theory, the state of the component models should fully determine the state of
the coupled simulation and one should be able to restart a coupled simulation
just by using the states of the component models. Unfortunately, this is
currently not the case in `ClimaCoupler`. The main reason for this is the
complex interdependencies between component models and within `ClimaAtmos` which
make the initialization step inconsistent. For example, in a coupled simulation,
the surface albedo should be determined by the surface models and used by the
atmospheric model for radiation transfer, but `ClimaAtmos` also tries to set the
surface albedo (since it has to do so when run in standalone mode). In addition
to this, `ClimaAtmos` has a large cache that has internal interdependencies that
are hard to disentangle, and changing a field might require changing some other
field in a different part of the cache. As a result, it is not easy for
`ClimaCoupler` to consistently do initialization from a cold state. To conclude,
restarting a simulation exclusively using the states of the component models is
currently impossible.

Given that restarting a simulation from the state is impossible, `ClimaCoupler`
needs to save the states and the caches. Let us review how we use
`ClimaCore.InputOutput` and `JLD2` package to accomplish this.

`ClimaCore.InputOutput` provides a loss-less way to save the content of certain
`ClimaCore` objects to HDF5 files. Objects saved in this way are not tied to a
particular computing device or configuration. When running with MPI,
`ClimaCore.InputOutput` are also efficiently written in parallel.

Unfortunately, `ClimaCore.InputOutput` only supports certain objects, such as
`Field`s and `Space`s, but the cache in component models is more complex than
this and contains complex objects with highly stateful quantities (e.g., C
pointers). Because of this, model states are saved to HDF5 but caches must be
saved to JLD2 files.

`JLD2` allows us to save more complex objects without writing specific
serialization methods for every struct. `JLD2` allows us to take a big step
forward, but there are still several challenges that need to be solved:
1. `JLD2` does not support CUDA natively. To go around this, we have to move
   everything onto the CPU first. Then, when the data is read back, we have to
   move it back to the GPU.
2. `JLD2` does not support MPI natively. To go around this, each process writes
   its `jld2` checkpoint and reads it back. This introduces the constraint that
   the number of MPI processes cannot change across restarts.
3. Some quantities are best not saved and read (for example, anything with
   pointers). For this, we write a recursive function that traverses the cache
   and only restores quantities of a certain type (typically, `ClimaCore`
   objects)

Point 3. adds significant amount of code and requires component models to
specify how their cache has to be restored.

If you are adding a component model, you have to extend the methods.
```
Checkpointer.get_model_prog_state
Checkpointer.get_model_cache
Checkpointer.restore_cache!
```

`ClimaCoupler` moves objects to the CPU with `Adapt(Array, x)`. `Adapt`
traverses the object recursively, and proper `Adapt` methods have to be defined
for every object involved in the chain. The easiest way to do this is using the
`Adapt.@adapt_structure` macro, which defines a recursive Adapt for the given
object.

Types to watch for:
- `MPI` related objects (e.g., `MPICommsContext`)
- `TimeVaryingInputs` (because they contain `NCDatasets`, which contain pointers
  to files)

!!! warning "Adapt and references"
    For objects that contain multiple fields referencing the same object, using
    the `Adapt.@adapt_structure` macro leads to unnecessary copies of the same
    object. This happens because `Adapt.@adapt_structure` defines a recursive
    `Adapt` that does not account for the possibility that multiple fields could
    be referencing the same object. As a result, this means that the same object
    is recreated over and over again when calling `Adapt` on the cache. This can
    easily make the file size of the saved cache much bigger than it needs to
    be. Because of this, we've implemented a `CacheIterator` object - please see
    the section below for details.

### `CacheIterator`

Instead of defining a proper `Adapt` method for the cache, an alternative
approach is to recursively iterate over the cache fields and selectively save
only the parts that need to be saved. This recursive iteration is performed by
the `CacheIterator`. To initialize a `CacheIterator` for a component model, you
must implement `get_cache_ignore`.

Using the `CacheIterator` allows `adapt` to be called on each individual field
instead of on the entire cache. Furthermore, the file size can be reduced by
avoiding duplicate saves of fields that reference the same memory. This is
accomplished by tracking object IDs and storing references to objects instead of
creating copies when the same object is encountered multiple times.

This approach allows for a signficant reducation in the file size of the cache.

## Checkpointer API

```@docs
    Checkpointer.get_model_prog_state
    Checkpointer.get_model_cache
    Checkpointer.get_model_cache_to_checkpoint
    Checkpointer.restart!
    Checkpointer.checkpoint_sims
    Checkpointer.t_start_from_checkpoint
    Checkpointer.restore!
```
