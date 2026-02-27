# The CoupledSimulation Object

`CoupledSimulation` is the central object in ClimaCoupler.jl. It holds everything
needed to run a coupled simulation: the component models, exchange fields,
time-stepping information, output directories, diagnostics, and callbacks. It is
created by the `CoupledSimulation` constructor and passed to `run!`, `step!`,
or `postprocess`.

Here we will describe the structure and internals of `CoupledSimulation`.

## `CoupledSimulation` Fields

| Field                 | Type | Description |
|-----------------------|------|-------------|
| `model_sims`          | `NamedTuple` | All component model simulation objects, keyed by name (see below) |
| `fields`              | `ClimaCore.Fields.Field` | Shared exchange-grid fields passed between component models (see below) |
| `t`                   | `Ref` of `ITime` or `Float64` | Current simulation time (Float64 seconds or `ITime`) |
| `tspan`               | `Tuple` of `Float64` | Start and end times of the simulation |
| `Δt_cpl`              | `Int` or `Float64` | Coupling timestep in seconds |
| `start_date`          | `Ref{Dates.DateTime}` | Calendar date corresponding to `t = t_start` |
| `dir_paths`           | `NamedTuple` of `String` | Output, checkpoint, artifacts, and per-model output directories |
| `callbacks`           | `NamedTuple` of `TimeManager.Callback` | Scheduled callbacks triggered each coupling step (e.g. checkpointing) |
| `conservation_checks` | `NamedTuple{(:energy,:water)}` or `nothing` | Energy and water conservation accumulators; `nothing` when disabled |
| `diags_handler`       | `ClimaDiagnostics.DiagnosticsHandler` or `nothing` | ClimaDiagnostics.jl handler that schedules and writes diagnostics |
| `thermo_params`       | `Thermodynamics.Parameters.ThermodynamicsParameters` or `nothing` | Thermodynamic parameters shared across component models |
| `save_cache`          | `Bool` | Whether model caches are included when writing checkpoint files |

## Component model simulations (`cs.model_sims`)

`cs.model_sims` is a `NamedTuple` containing one entry per active component model. Only
models that are present in the simulation appear in the tuple — for example, an aquaplanet
simulation has no `land_sim` or `ice_sim`. The standard keys are:

| Key | Model |
|---|---|
| `atmos_sim` | Atmosphere model |
| `land_sim` | Land model (if present) |
| `ocean_sim` | Ocean model (if present) |
| `ice_sim` | Sea ice model (if present) |

Each value is a concrete subtype of `AbstractComponentSimulation`. You can iterate over
all component models with:

```julia
for sim in cs.model_sims
    @info nameof(sim), extrema(Interfacer.get_field(sim, Val(:surface_temperature)))
end
```

See [Available simulation types](@ref) for details on which models are present in each
simulation type.

## Coupler exchange fields (`cs.fields`)

`cs.fields` is a `ClimaCore.Fields.Field` of `NamedTuple`s defined on the shared
boundary (exchange) grid. It holds the intermediate values passed between the atmosphere
and surface models at each coupling step. Individual fields are accessed as properties:

```julia
cs.fields.T_sfc        # combined surface temperature on the exchange grid
cs.fields.F_lh         # latent heat flux
extrema(cs.fields.SW_d)  # check shortwave radiation range
```

The full list of default fields and their descriptions is documented in
[`Interfacer`](@ref).
Component models can register additional fields beyond the defaults by extending
`Interfacer.add_coupler_fields!`. The full list of exchange fields of a simulation can be
inspected with `propertynames(cs.fields)`.

## Output directories (`cs.dir_paths`)

`cs.dir_paths` is a `NamedTuple` of directory paths set up at initialization,
with the following structure and keys in the `NamedTuple`:

```
output_dir_root/           → `output_dir_root`
├── clima_coupler/         → `coupler_output_dir`
├── clima_atmos/           → `atmos_output_dir`
├── clima_land/            → `land_output_dir`
├── clima_ocean/           → `ocean_output_dir`
├── clima_seaice/          → `ice_output_dir`
├── artifacts/             → `artifacts_dir`
├── checkpoints/           → `checkpoints_dir`
└── regrid_tmp_<random>/   → `regrid_dir`
```

See [`Utilities.setup_output_dirs`](@ref) for how these directories are created and
how to customise the root path via the `coupler_output_dir` configuration option.

## Callbacks (`cs.callbacks`)

`cs.callbacks` is a `NamedTuple` of `TimeManager.Callback` objects, which each execute
some operation at predetermined times throughout a simulation. Each callback
wraps a schedule and a function, and is triggered at the appropriate time during
`step!`. The default callbacks include a walltime progress reporter and, if enabled,
a checkpointing callback. See [TimeManager](@ref) for the `Callback` API and
[Checkpointer](@ref) for how to add a checkpointing callback.

## Conservation checks (`cs.conservation_checks`)

When energy and water conservation checking is enabled (via the `energy_check`
configuration option), `cs.conservation_checks` is a `NamedTuple` with keys
`:energy` and `:water`, holding an `EnergyConservationCheck` and a
`WaterConservationCheck` respectively. Otherwise it is `nothing`. Conservation
checks are only supported for slabplanet simulations. See
[ConservationChecker](@ref) for details.

## Diagnostics handler (`cs.diags_handler`)

`cs.diags_handler` is a `ClimaDiagnostics.DiagnosticsHandler` that schedules and
writes diagnostic output variables throughout the simulation. It is `nothing` when
`use_coupler_diagnostics` is disabled. See [Diagnostics](@ref) for how to
configure which variables are saved and at what frequency.

## Thermodynamic parameters (`cs.thermo_params`)

`cs.thermo_params` is a `Thermodynamics.Parameters.ThermodynamicsParameters` object
shared across all component models. It is used by `FluxCalculator` when computing
turbulent surface fluxes. In testing contexts only it may be `nothing`.

## Checkpoint settings (`cs.save_cache`, `cs.prev_checkpoint_t`)

`cs.save_cache` controls whether model caches (in addition to states) are written
when a checkpoint is saved. Setting this to `false` produces smaller checkpoint
files but prevents exact cache restoration on restart. `cs.prev_checkpoint_t` is a
`Ref` tracking the simulation time of the last checkpoint, used internally to remove
intermediate checkpoint files. See [Checkpointer](@ref) for the full checkpointing
documentation.
