# Subseasonal Calibration

Calibrate ClimaCoupler parameters using EnsembleKalmanProcesses against ERA5 observations.

## Quick Start


### Option A: TransformUnscented (recommended, simpler)

```bash
cd /glade/u/home/cchristo/clima/copies3/ClimaCoupler.jl

# 1. Generate observations (only needed once, or when changing obs settings)
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal/generate_observations.jl

# 2. Run calibration (from tmux on login node)
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal/run_calibration.jl
```

### Option B: TransformInversion (more robust, needs precompute)

```bash
cd /glade/u/home/cchristo/clima/copies3/ClimaCoupler.jl

# 1. Generate observations (only needed once)
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal/generate_observations.jl

# 2. Run full calibration (from tmux - handles precompute automatically)
./experiments/calibration/subseasonal/run_full_calibration.sh
```

## Switching Between Process Types

Edit `run_calibration.jl` and comment/uncomment the appropriate OPTION block (~line 116-156).

## Changing Parameters

Edit `calibration_priors.jl` - this is the single source of truth for:
- `CALIBRATION_PRIORS` - list of parameters to calibrate
- `CALIBRATION_ENSEMBLE_SIZE` - ensemble size (for TransformInversion only)

# Analizing Calibraiton 
    `julia --project=experiments/ClimaEarth experiments/calibration/subseasonal/analyze_calibration.jl`

## Key Files

| File | Purpose |
|------|---------|
| `calibration_priors.jl` | **Edit this** - defines parameters to calibrate |
| `run_calibration.jl` | Main calibration script |
| `generate_observations.jl` | Loads ERA5 data, creates observation vector |
| `observation_map.jl` | Maps model output to observation space |
| `run_full_calibration.sh` | Orchestrator for TransformInversion |

## Output

Results go to the `output_dir` specified in `run_calibration.jl` (default: `/glade/derecho/scratch/cchristo/calibration/exp9`).

Plot results:
```bash
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal/quick_plot_cal.jl
```

## Troubleshooting

**"Killed" on login node**: Use TransformUnscented, or run `./run_full_calibration.sh` which uses cpudev for heavy computation.

**NetCDF errors**: Retry logic is built in. If persistent, check that model output files exist.

**Bus error**: JLD2 mmap issue on Lustre - already mitigated by using IOStream mode.
