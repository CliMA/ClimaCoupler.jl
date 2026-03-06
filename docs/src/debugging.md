# Debugging Tips

Coupled simulations can fail or produce incorrect results for several broad reasons:

- **Package incompatibilities or versioning errors**: A dependency conflict prevents the
  simulation from starting at all or causes unexpected behavior.
- **Numerical instabilities**: The simulation starts but crashes during the run
  due to unphysical values propagating through the component models.
- **Software errors**: A bug in the coupling code causes incorrect dispatch, wrong field
  updates, or silent failures.

Here we suggest some approaches for debugging each of these cases.

## Package incompatibilities

Julia resolves package versions at environment instantiation time. Because ClimaCoupler.jl
depends on a large number of packages from the CliMA ecosystem and beyond, it can be
difficult to keep all compat entries simultaneously satisfied.

**Tips for maintaining compatibility:**

- Keep `[compat]` entries in `Project.toml` files as broad as possible in ClimaCoupler.jl and
  upstream packages. Only restrict a version range when a specific breaking change, bug, or
  API incompatibility requires it. Overly narrow compat entries are a common source of
  resolution failures across the ecosystem.

**Debugging a broken environment:**

If Julia cannot resolve or load your environment, try the following steps in order:

```julia
# 1. Make sure the environment is fully instantiated
import Pkg
Pkg.instantiate()

# 2. Try resolving to enforce version constraints from scratch
Pkg.resolve()

# 3. Update all packages to the latest compatible versions
Pkg.update()
```

If the above steps don't help, the issue is often caused by one specific package being
pinned to an incompatible version. A somewhat-reliable workaround is to remove the
offending package and re-add it so the resolver can pick a fresh compatible version:

```julia
Pkg.rm("OffendingPackage")
Pkg.add("OffendingPackage")
```

Sometimes two packages constrain each other in conflicting ways; in that case, remove
both at once before re-adding them. You can also run `Pkg.status()` to inspect version
pins and `Pkg.resolve()` to see the full error message from the resolver.

## Numerical instabilities

A common symptom is a simulation that starts successfully but crashes at a particular
timestep with a NaN or domain error. A common debugging workflow for this is:

1. Run the simulation once to identify the failing timestep.
2. Start a new simulation interactively and advance it via `step!` until just before the
   crash, then inspect coupler fields and component model states.

```julia
# Run from the top-level ClimaCoupler.jl/ directory
include("experiments/ClimaEarth/code_loading.jl")

cs = CoupledSimulation()

# Advance to one step before the known failure
for i in 1:failing_step - 1
    step!(cs)
end

# Now inspect coupler fields for unphysical values
@info extrema(cs.fields.T_sfc)
@info extrema(cs.fields.F_turb_moisture)
@info extrema(cs.fields.F_lh)
```

Look for: negative temperatures, extremely large or small fluxes, NaN
values in any field, or values that are zero where they shouldn't be.

### Plotting ClimaCore fields

The built-in `Plotting.debug` function produces a grid of heatmaps of all coupler and
component model fields, and is the easiest first step:

```julia
# Plot all coupler fields and component model fields to a directory
Plotting.debug(cs, "debug_output")
# Produces: debug_output/debug_coupler.png, debug_output/debug_<ModelName>.png, etc.
```

You can also plot individual ClimaCore fields with ClimaCoreMakie directly:

```julia
using ClimaCoreMakie, CairoMakie, Makie

field = cs.fields.T_sfc  # any ClimaCore surface field

fig = Figure()
ax = Axis(fig[1, 1], title = "T_sfc $(extrema(field))")
hm = fieldheatmap!(ax, field)
Colorbar(fig[:, end+1], hm)
Makie.save("T_sfc.png", fig)
```

For 3D fields (e.g. from the atmosphere), extract a level first:

```julia
import ClimaCore as CC

field_3d = cs.model_sims.atmos.integrator.u.c.ρe_tot  # example 3D field
field_sfc = CC.Fields.level(field_3d, 1)               # bottom level

fig = Figure()
ax = Axis(fig[1, 1], title = "ρe_tot (level 1)")
hm = fieldheatmap!(ax, field_sfc)
Colorbar(fig[:, end+1], hm)
save("rhoe_sfc.png", fig)
```

### Plotting Oceananigans fields

When running a CMIP simulation with `OceananigansSimulation`, use `CairoMakie.heatmap!`
on a horizontal slice of the Oceananigans field. Note that this applies
to `ClimaSeaIceSimulation` fields too, since it uses Oceananigans.jl under the hood.

For a 2D field we just access the underlying array:
```julia
using CairoMakie, Makie
import Oceananigans as OC

ice = cs.model_sims.ice_sim.ice
ice_T = ice.model.ice_thermodynamics.top_surface_temperature # Oceananigans field
ice_T_array = OC.interior(ice_T, :, :, 1) # 2D underlying array

fig = Figure()
ax = Axis(fig[1, 1], title = "Ice surface T")
hm = heatmap!(ax, ice_T_array)
Colorbar(fig[:, end+1], hm)
Makie.save("ice_T_sfc.png", fig)
```

For a 3D field we need to index into the top level of the grid:
```julia
using CairoMakie, Makie
import Oceananigans as OC

ocean = cs.model_sims.ocean_sim.ocean
ocean_T_3d = ocean.model.tracers.T
# Index into the surface (top) level of the 3D field
ocean_T_surface = OC.interior(ocean_T_3d, :, :, ocean.model.grid.Nz)

fig = Figure()
ax = Axis(fig[1, 1], title = "Ocean surface T")
hm = heatmap!(ax, ocean_T_surface)
Colorbar(fig[:, end+1], hm)
Makie.save("ocean_T_sfc.png", fig)
```

For `AbstractOperation` fields (derived quantities), compute them first
then index as appropriate for the field's domain:

```julia
field_3d = OC.Field(some_operation)
OC.compute!(field_3d)
CairoMakie.heatmap!(ax, view(field_3d, :, :, derived.grid.Nz))
```

The `Plotting.debug` function handles both field types automatically when the
`ClimaCouplerCMIPMakieExt` extension is loaded.

## Software errors

### Identifying which method is being dispatched

A common source of subtle bugs is a function being called with an unexpected type,
causing Julia to dispatch to the wrong method. If you suspect this is the case, use
`@which` to confirm which method will be called:

```julia
# See exactly which method Julia will call for these argument types
@which FieldExchanger.update_sim!(sim, coupler_fields, area_fraction)

# List all defined methods of a function
methods(FieldExchanger.update_sim!)

# List all methods applicable to a specific type
import InteractiveUtils: methodswith
methodswith(typeof(sim), FieldExchanger.update_sim!)
```

This is particularly useful when adding a new component model: if a required interface
function (e.g. `Interfacer.get_field`, `Interfacer.update_field!`) is not extended for
your new type, Julia will silently dispatch to a fallback method rather than erroring.
`@which` will reveal this immediately.

### Checking return values and field updates

If a field is unexpectedly zero or unchanged after a coupling step, trace back through
the call chain by inserting `@info` statements or using the Julia debugger:

```julia
import Debugger
Debugger.@enter FieldExchanger.exchange!(cs)
```

You can also use `@show` to quickly print intermediate values without a full debugger
session:

```julia
@show extrema(cs.fields.T_sfc)
@show nameof(typeof(cs.model_sims.ocean_sim))
```

### Interactive debugging with Infiltrator.jl

[Infiltrator.jl](https://github.com/JuliaDebug/Infiltrator.jl) lets you drop into an
interactive REPL at any point in the code, even entering functions from external packages.
This is especially useful for coupled simulations where the call stack spans
multiple packages.

!!! note
    Infiltrator.jl must be installed in your default julia environment, not in the
    experiment project. Julia loads it from the default environment automatically.

#### Pausing inside the coupling loop

Place `@infiltrate` at the point you want to inspect. For example, to pause after every
coupling step:

```julia
import Infiltrator

# After loading and constructing the CoupledSimulation `cs`:
for i in 1:n_steps
    step!(cs)
    Infiltrator.@infiltrate  # drops you into a REPL here
end
```

Inside the Infiltrator REPL you have access to all local variables (`cs`, `i`, etc.)
and can run any Julia code:

```
infil> extrema(cs.fields.T_sfc)
(215.3, 318.7)
infil> cs.t[]
1200
```

Type `@continue` to resume execution, or `@exit` to abort.

#### Conditional infiltration

If the crash happens at a specific timestep, use a condition to avoid stopping at every
iteration:

```julia
for i in 1:n_steps
    step!(cs)
    Infiltrator.@infiltrate i == 42  # only pause at step 42
end
```

You can also trigger on a data condition, such as the first appearance of NaN:

```julia
for i in 1:n_steps
    step!(cs)
    Infiltrator.@infiltrate any(isnan, parent(cs.fields.T_sfc))
end
```

#### Infiltrating inside library code

If the problem is inside a function you don't own (e.g. in `FieldExchanger` or a
component model extension), you can temporarily redefine that function with an
`@infiltrate` call. For example, to pause inside `update_sim!` for a specific model type:

```julia
import Infiltrator
import ClimaCoupler: FieldExchanger, Interfacer

function FieldExchanger.update_sim!(
    sim::Interfacer.AbstractSurfaceSimulation,
    fields,
    area_fraction,
)
    Infiltrator.@infiltrate
    # ... original body (copy-paste from source) ...
end
```

From the Infiltrator REPL, use `@trace` to see the full call stack, which helps identify
which caller led to unexpected state.

#### Exfiltrating variables to the main REPL

Sometimes it's easier to pull variables out of a deep scope for post-hoc analysis rather
than debugging inline. Use `@exfiltrate` to copy local variables into
`Infiltrator.safehouse`, which persists after you `@continue`:

```julia
import Infiltrator

for i in 1:n_steps
    step!(cs)
    if any(isnan, parent(cs.fields.T_sfc))
        Infiltrator.@exfiltrate
        break
    end
end

# Back in the main REPL:
(; cs, i) = Infiltrator.safehouse
extrema(cs.fields.F_turb_moisture)
```

For more details about using Infiltrator with ClimaCore objects,
including using with ClimaCore's `DebugOnly` callback system for catching
NaNs inside low-level field operations, see the
[ClimaCore debugging guide](https://clima.github.io/ClimaCore.jl/stable/debugging/).

For more general information, see the [Infiltrator.jl docs](https://juliadebug.github.io/Infiltrator.jl/stable/).
