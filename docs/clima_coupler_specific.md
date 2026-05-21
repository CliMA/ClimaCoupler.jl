# ClimaCoupler.jl — Repo-Specific Guide

<!-- This file documents things specific to ClimaCoupler.jl that are NOT
     covered by the shared DeveloperGuides at docs/dev-guides/. -->

## Package overview

ClimaCoupler.jl provides coupled-system time-stepping control and boundary-field
exchange between component models (atmosphere, ocean, land, sea ice). It is
agnostic to component internals: each component exposes a small interface defined
in `src/Interfacer.jl`, and the coupler orchestrates time stepping and field
mapping through `src/SimCoordinator.jl`. The primary target use-cases are
[AMIP](https://pcmdi.llnl.gov/mips/amip/home/overview.html) and
[CMIP](https://wcrp-cmip.org)-style global climate simulations, though the
design also supports lighter slabplanet configurations useful for conservation
testing.

## Directory map

Mapped to the architectural layers in
[architectural_boundaries.md](dev-guides/architecture/architectural_boundaries.md).

| Layer | Directory / File | Description |
|:---|:---|:---|
| Coupler Core | `src/Interfacer.jl` | Component-model interface definitions (`SurfaceModelSimulation`, `AtmosModelSimulation`, `CoupledSimulation`) |
| Coupler Core | `src/FieldExchanger.jl` | Import/export field mapping between component models |
| Coupler Core | `src/FluxCalculator.jl` | Turbulent surface-flux computation at the atmosphere–surface boundary |
| Coupler Core | `src/SimCoordinator.jl` | Top-level time-stepping loop and component orchestration |
| Coupler Core | `src/ConservationChecker.jl` | Energy and water conservation diagnostics |
| Coupler Core | `src/TimeManager.jl` | Simulation calendar and time-stepping helpers |
| Coupler Core | `src/Checkpointer.jl` | JLD2-based checkpoint save/load |
| Coupler Core | `src/Input.jl` | YAML configuration loading and parameter dispatch |
| Coupler Core | `src/Utilities.jl` | Shared utility functions (regridding, masking, etc.) |
| Coupler Core | `src/CalibrationTools.jl` | EKI/EnsembleKalmanProcesses helpers |
| Simple Models | `src/Models/` | Built-in stub component models (slab ocean, prescribed ocean/sea ice) |
| Simple Models | `src/surface_stub.jl` | Minimal surface stub for testing |
| Output | `src/SimOutput/` | Diagnostics output helpers |
| Output | `src/Plotting.jl` | Makie-based post-processing plots (loaded via weak dep) |
| Component Extensions | `ext/ClimaCouplerClimaAtmosExt.jl` | ClimaAtmos component-model integration |
| Component Extensions | `ext/ClimaCouplerClimaLandExt/` | ClimaLand component-model integration |
| Component Extensions | `ext/ClimaCouplerCMIPExt/` | CMIP-mode (Oceananigans + ClimaSeaIce) integration |
| Component Extensions | `ext/ClimaCouplerMakieExt/` | Visualization extensions |
| Driver Scripts | `experiments/AMIP/` | AMIP run scripts and configuration |
| Driver Scripts | `experiments/CMIP/` | CMIP run scripts and configuration |
| Driver Scripts | `experiments/calibration/` | EKI calibration experiment drivers |
| Configuration | `config/` | YAML configuration files (CI configs, experiment presets) |
| Tests | `test/` | Unit and integration tests (see Test groups below) |
| Documentation | `docs/` | Documenter.jl setup; `docs/dev-guides/` is the DeveloperGuides subtree |

## Key abstractions

1. **`AbstractSurfaceSimulation` / `AbstractAtmosSimulation`** (`src/Interfacer.jl`) — abstract types that component models subtype to participate in coupling. Implement the interface methods (e.g., `get_field`, `update_field!`, `step!`) to plug a new component in. Surface types further specialize into `AbstractOceanSimulation`, `AbstractLandSimulation`, and `AbstractSeaIceSimulation`.

2. **`CoupledSimulation`** (`src/Interfacer.jl`) — the top-level container holding all component simulations, shared fields, and the ClimaComms context.

3. **`FieldExchanger` module** (`src/FieldExchanger.jl`) — maps fields from each component into the coupler's shared field space and back, handling regridding and masking via `ClimaCore` and `ClimaUtilities`. Key function: `update_sim!`.

4. **`FluxCalculator` module** (`src/FluxCalculator.jl`) — computes turbulent surface fluxes (latent, sensible, momentum) using `SurfaceFluxes.jl` at the atmosphere–surface boundary. Key function: `turbulent_fluxes!`.

5. **`SimCoordinator` module** (`src/SimCoordinator.jl`) — drives the outer coupling loop: calls each component's `step!`, triggers field exchanges, writes checkpoints, and runs conservation checks. Key functions: `run!`, `step!`, `setup_and_run`.

## Running experiments

### AMIP (standard)
```bash
julia --project=experiments/AMIP -E "using Pkg; Pkg.instantiate()"
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl \
    --config_file config/ci_configs/amip_default.yml --job_id amip_default
```
Output lands in `output/<job_id>/`. Default run takes ~10 min on a single CPU.

### Slabplanet (conservation testing)
```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl \
    --config_file config/ci_configs/slabplanet_default.yml --job_id slabplanet_default
```
Set `energy_check: true` in the config to enable conservation tracking.

### Environment variables
| Variable | Values | Effect |
|:---|:---|:---|
| `CLIMACOMMS_DEVICE` | `CPU` (default), `CUDA` | Switch between CPU and GPU |
| `CLIMACOMMS_CONTEXT` | `SINGLETON` (default), `MPI` | Disable or enable MPI |

### Caltech HPC (Central / clima)
```bash
# Central
export MODULEPATH="/resnick/groups/esm/modules:$MODULEPATH"
module load climacommon

# clima
module load common
```
See the [slurm-buildkite wiki](https://github.com/CliMA/slurm-buildkite/wiki) for cluster access details.

## Test groups

Run the full suite with:
```bash
julia --project=test test/runtests.jl
```

| Test file | What it covers |
|:---|:---|
| `aqua.jl` | Aqua.jl package-quality checks (ambiguities, unbound type params, stale deps) |
| `interfacer_tests.jl` | Component interface contracts, `get_field`/`update_field!` dispatch |
| `conservation_checker_tests.jl` | Energy and water conservation bookkeeping |
| `utilities_tests.jl` | Regridding, masking, and shared utility functions |
| `field_exchanger_tests.jl` | Field mapping between atmosphere and surface |
| `flux_calculator_tests.jl` | Turbulent surface flux computations |
| `input_tests.jl` | YAML config loading and parameter parsing |
| `sim_output_tests.jl` | Diagnostics output helpers |
| `models/slab_ocean_tests.jl` | Slab ocean built-in model |
| `models/prescr_ocean_tests.jl` | Prescribed ocean built-in model |
| `models/prescr_seaice_tests.jl` | Prescribed sea-ice built-in model |
| `calibration_tools_tests.jl` | EKI/calibration helpers |
| `mpi_tests/` | MPI-parallel coupling tests (run separately via Buildkite) |

CI is managed by Buildkite. The GitHub Actions matrix (`.github/workflows/ci.yml`) runs Julia 1.10 and 1.12 on CPU.

## Repo-specific conventions

- **Configuration-driven experiments**: All experiment parameters live in YAML files under `config/`. Prefer adding a new config file over hard-coding values in driver scripts.
- **Weak dependencies for heavy visualisation**: `Makie`, `CairoMakie`, `GeoMakie`, `Oceananigans`, etc. are `[weakdeps]`. Extensions in `ext/` are loaded only when the user explicitly imports those packages. Keep the core coupler free of heavy visualization or ocean-model imports.
- **No component model state inside the coupler core**: `src/` must remain agnostic to the internals of ClimaAtmos, ClimaLand, etc. Component-specific logic belongs in the corresponding `ext/` extension.
- **ClimaComms device/context agnosticism**: All coupler code must run on CPU (SINGLETON), CPU+MPI, and GPU (CUDA). Use `ClimaComms.context`, `ClimaComms.device`, and `ClimaComms.iamroot` rather than device-specific intrinsics.
- **News entries**: Breaking changes and notable features must be recorded in `NEWS.md` following the format in [changelogs_and_versions.md](dev-guides/code-quality/changelogs_and_versions.md).

## Self-correction

If this guide is discovered to be stale or missing a pattern, update it.
