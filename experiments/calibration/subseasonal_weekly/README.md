# Subseasonal Calibration

Calibrate ClimaCoupler parameters using EnsembleKalmanProcesses against ERA5 observations.

## Quick Start


### Option A: TransformUnscented (recommended, simpler)

```bash
cd /glade/u/home/cchristo/clima/copies2/ClimaCoupler.jl

# 1. Generate observations (only needed once, or when changing obs settings)
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal_weekly/generate_observations.jl

# 2. Run calibration (from tmux on login node)
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal_weekly/run_calibration.jl
```

### Option B: TransformInversion (more robust, needs precompute)

```bash
cd /glade/u/home/cchristo/clima/copies2/ClimaCoupler.jl

# 1. Generate observations (only needed once)
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal_weekly/generate_observations.jl

# 2. Run full calibration (from tmux - handles precompute automatically)
./experiments/calibration/subseasonal_weekly/run_full_calibration.sh
```

## Switching Between Process Types

Edit `run_calibration.jl` and comment/uncomment the appropriate OPTION block (~line 116-156).

## Changing Parameters

Edit `calibration_setup.jl` - this is the single source of truth for:
- `CALIBRATION_PRIORS` - list of parameters to calibrate
- `CALIBRATION_ENSEMBLE_SIZE` - ensemble size (for TransformInversion only)
- `NORMALIZE_VARIABLES` - if true, normalize variables to zero mean and unit variance
- `CALIBRATION_NOISE_SCALAR` - noise variance for observation covariance

### Variable Normalization

When `NORMALIZE_VARIABLES = true` (default), each variable is normalized to zero mean and unit variance:
1. `generate_observations.jl` computes latitude-weighted mean and std from ERA5 data
2. Stats are saved to `norm_stats.jld2` 
3. Observations are normalized before creating the observation vector
4. `observation_map.jl` loads the same stats and applies them to model output

This ensures variables with different magnitudes (e.g., temperature in K vs radiation in W/m²) contribute equally to the calibration objective.

## Data Source

Uses weekly-averaged ERA5 data from:
`/glade/campaign/univ/ucit0011/cchristo/wxquest_data/daily_weekly_stats/weekly`

# Analyzing Calibration 
    `julia --project=experiments/ClimaEarth experiments/calibration/subseasonal_weekly/analyze_calibration.jl`

## Key Files

| File | Purpose |
|------|---------|
| `calibration_setup.jl` | defines parameters and settings to calibrate |
| `run_calibration.jl` | Main calibration script |
| `generate_observations.jl` | Loads weekly ERA5 data, creates observation vector |
| `observation_map.jl` | Maps model output to observation space |
| `run_full_calibration.sh` | Orchestrator for TransformInversion |

## Output

Results go to the `output_dir` specified in `run_calibration.jl` (currently `/glade/derecho/scratch/cchristo/calibration/expXX`).
Logs are saved to `calibration_YYYYMMDD_HHMMSS.log` in the output directory.

Plot results:
```bash
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal_weekly/quick_plot_cal.jl
```

## Troubleshooting

**"Killed" on login node**: Use TransformUnscented, or run `./run_full_calibration.sh` which uses cpudev for heavy computation.

**NetCDF errors**: Retry logic is built in. If persistent, check that model output files exist.

**Bus error**: JLD2 mmap issue on Lustre - already mitigated by using IOStream mode.
