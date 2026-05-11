# ML Flux Correction for ClimaAtmos EDMF

## What this does

ClimaAtmos uses an EDMF (Eddy Diffusivity / Mass Flux) turbulence scheme to represent sub-grid scale vertical transport. This scheme has two components: eddy diffusivity (ED) and mass flux (MF). Both are imperfect. This code trains a neural network to learn a **third flux correction term** that closes the gap between the model and ERA5 reanalysis.

The key idea: instead of predicting tendencies directly, the NN predicts a **vertical flux profile** F(z) at cell faces. The tendency correction is then the flux divergence -dF/dz. This guarantees the correction is **conservative** (the column integral of the tendency equals the net flux through the boundaries, which we set to zero).

```
Total SGS tendency = ED term + MF term + NN flux correction
```
1. **Compute target**: For each hourly timestep, the correction tendency is the residual between ERA5 and model tendencies:
  ```
   dT_correction/dt = (T_era5(t+1) - T_era5(t))/dt - (T_model(t+1) - T_model(t))/dt
  ```
   Same for specific humidity q.
2. **Train the NN**: The network takes a column of atmospheric features and outputs flux profiles at cell faces for both T and q. The loss minimises:
  ```
   || -(F[k+1] - F[k]) / dz[k]  -  target_tendency[k] ||²
  ```
   with boundary conditions F[top] = F[bottom] = 0 enforced by construction.
3. **Deploy online**: Once trained, the flux correction enters ClimaAtmos as an additional SGS flux term alongside ED and MF.

## Data

- **Model output**: Hourly instantaneous snapshots from ClimaAtmos on a 192x96x70 (lon/lat/z) grid. NetCDF files at `{MODEL_DIR}/{varname}_1h_inst.nc`.
- **ERA5**: Hourly pressure-level data (37 levels, 0.25° global). One file per hour at `{ERA5_DIR}/era5_pressure_levels_YYYYMMDD_HH00.nc`. Fields: T, q, geopotential, u, v, w.
- **Regridding**: ERA5 is bilinearly interpolated to the model's horizontal grid, then vertically interpolated from pressure levels to model z-levels using geopotential height as the common coordinate.

## Architecture

Two options, selected via `ARCHITECTURE`:

- `**:unet`** (default) — 1D UNet over the vertical column. Two encoder stages (Conv + MaxPool) downsample the vertical dimension, a bottleneck processes at coarse resolution, then two decoder stages (ConvTranspose + skip connections) upsample back. A 1x1 conv produces T and q values at cell centres, which are linearly interpolated to interior face fluxes. This captures both local and nonlocal vertical structure.
- `**:mlp**` — Plain MLP on the flattened column (3 hidden layers + output). Simpler and faster at inference but has no structural inductive bias for the vertical dimension.

## Configuration

All configuration is via constants at the top of `train_flux_correction.jl`:


| Parameter                   | Default                       | Description                                     |
| --------------------------- | ----------------------------- | ----------------------------------------------- |
| `MODEL_DIR`                 | `.../output_0004/clima_atmos` | Path to ClimaAtmos hourly output                |
| `ERA5_DIR`                  | `.../ml_correct_v1`           | Path to ERA5 hourly pressure-level files        |
| `Z_MAX`                     | `20000.0` m                   | Only train on levels below this height          |
| `ARCHITECTURE`              | `:unet`                       | `:unet` or `:mlp`                               |
| `HIDDEN_DIM`                | `128`                         | Hidden dim (MLP) or base channels×4 (UNet)      |
| `N_EPOCHS`                  | `100`                         | Training epochs                                 |
| `BATCH_SIZE`                | `256`                         | Columns per batch                               |
| `LEARNING_RATE`             | `1e-3`                        | Adam learning rate                              |
| `LOSS_WEIGHT_Q`             | `1e6`                         | Upweights q loss (q is O(1e-3) vs T O(100))     |
| `TIME_DAY_START`            | `2`                           | First day of data to use (skip spinup)          |
| `TIME_DAY_END`              | `4`                           | Last day of data to use                         |
| `COLUMN_SUBSAMPLE_FRACTION` | `0.1`                         | Fraction of columns to randomly sample per step |
| `INPUT_FEATURES`            | see below                     | List of (varname, transform) pairs              |


### Input features

Configured as a list of `(netcdf_varname, transform)` tuples. To add a variable, append a line — everything else adapts automatically.

Current defaults (5 features):

```julia
("ta",    identity)                    # temperature
("hus",   identity)                    # specific humidity
("pfull", x -> log(max(x, 1f0)))      # log-pressure
("hur",   identity)                    # relative humidity
("tke",   identity)                    # turbulent kinetic energy
```

Other available variables to add: `edt` (eddy diffusivity), `rhoa` (air density), `arup` (updraft area fraction), `waup` (updraft velocity), `cl` (cloud fraction), `lmix` (mixing length), `clw`/`cli` (cloud liquid/ice).

## Running

**CPU** (for debugging / small subsets):

```bash
cd /home/cchristo/clima/ClimaCoupler.jl
julia --project=ml_flux -e 'using Pkg; Pkg.instantiate()'
julia --project=ml_flux ml_flux/train_flux_correction.jl
```

**GPU via SLURM**:

```bash
sbatch ml_flux/run_training.sbatch
```

GPU is auto-detected — the same script runs on CPU if no GPU is available.

## Output

The trained model is saved to `ml_flux/flux_correction_model.bson` containing:

- `model` — the Lux model architecture
- `ps` — trained parameters (on CPU)
- `st` — model state
- `X_norm_stats` — input normalisation (mean/std per feature)
- `nz` — number of vertical levels used

## Files

```
ml_flux/
├── README.md                     # this file
├── Project.toml                  # Julia dependencies
├── train_flux_correction.jl      # training script
├── run_training.sbatch           # SLURM job script (GPU)
└── flux_correction_model.bson    # saved model (after training)
```

