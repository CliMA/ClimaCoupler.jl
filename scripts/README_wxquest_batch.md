# WxQuest / WeatherBench batch submission (Derecho)

Submit one PBS job per initialization date. Default base config is subseasonal
`wxquest_progedmf`; pass `--base-config` for WeatherBench.

## Runbook

Run from repo root on Derecho (`/glade/u/home/cchristo/clima/copies2/ClimaCoupler.jl`).

### Submit — WeatherBench (example)

```bash
# Dry-run
python scripts/submit_wxquest_batch.py \
  --base-config config/weather/weatherbench_progedmf.yml \
  --dates 20101001 \
  --ic-dir /glade/campaign/univ/ucit0011/cchristo/wxquest_data/initial_conditions/initial_conditions_0p25deg_dev \
  --gpus 4 \
  --skip-ic-check \
  --dry-run

# Submit
python scripts/submit_wxquest_batch.py \
  --base-config config/weather/weatherbench_progedmf.yml \
  --dates 20101001 \
  --ic-dir /glade/campaign/univ/ucit0011/cchristo/wxquest_data/initial_conditions/initial_conditions_0p25deg_dev \
  --gpus 4 \
  --skip-ic-check \
  --yes
```

Notes for WeatherBench:
- Set `era5_initial_condition_dir` in the base YAML to the IC tree you want.
  The batcher does **not** overwrite it (`--ic-dir` is validation only).
- `--skip-ic-check` is needed when the IC tree has processed files only
  (missing raw intermediates like `era5_raw`, `era5_surface`, …).
- Follow the log: `tail -f generated/batch_wxquest_progedmf/logs/ccwx_progedmf_20101001_0000.out`
- Run output: `/glade/derecho/scratch/cchristo/coupler_runs_v2/weatherbench/20101001`

### Submit — subseasonal wxquest_progedmf

```bash
# Preview planned jobs (no writes, no qsub)
python scripts/submit_wxquest_batch.py --dry-run

# Test one date (nightly env, default)
python scripts/submit_wxquest_batch.py --dates 19970709 --yes

# Full default batch — 10 dates in DEFAULT_START_DATES
python scripts/submit_wxquest_batch.py --yes

# Remaining dates after a successful test (skip the test date)
python scripts/submit_wxquest_batch.py \
  --dates 19820615 20150627 20230711 19870731 20020617 20150514 20140622 20020731 19930601 \
  --yes

# Switch env profile
python scripts/submit_wxquest_batch.py --env nightly --yes   # main branches (default)
python scripts/submit_wxquest_batch.py --env amip --yes      # Manifest.toml pins
```

### Monitor

```bash
# Queue status (your jobs)
qstat -u $USER

# Follow a job log
tail -f generated/batch_wxquest_progedmf/logs/ccwx_progedmf_20101001_0000.out

# Cancel a job
qdel <jobid>
```

### Interactive (login node)

```bash
module load climacommon/2025_02_25
source scripts/setup_amip_env.sh nightly   # or: amip

julia --project=experiments/AMIP experiments/AMIP/run_simulation.jl \
  --config_file config/weather/weatherbench_progedmf.yml
# or: config/subseasonal_configs/wxquest_progedmf.yml
```

### Paths

| What | Where |
|------|-------|
| Generated configs | `generated/batch_wxquest_progedmf/configs/` |
| PBS scripts | `generated/batch_wxquest_progedmf/pbs/` |
| Job logs | `generated/batch_wxquest_progedmf/logs/ccwx_progedmf_<YYYYMMDD>_<HHMM>.out` |
| Run output | `coupler_output_dir/<YYYYMMDD>/` from the base YAML (scratch) |
| IC files | path in `era5_initial_condition_dir` (validate with `--ic-dir`) |

## Scripts

| File | Purpose |
|------|---------|
| `submit_wxquest_batch.py` | Generate per-date configs + PBS scripts and `qsub` them |
| `setup_amip_env.sh` | Set up `experiments/AMIP` Julia env (used by PBS jobs and interactively) |

## Environment profiles

Both profiles use `climacommon/2025_02_25`. The difference is Julia package resolution:

| `--env` | Packages |
|---------|----------|
| `nightly` (default) | `instantiate` + main on same set as Buildkite nightly: ClimaAtmos, ClimaCore, ClimaCoreMakie, ClimaTimeSteppers, Thermodynamics, ClimaLand, SurfaceFluxes, RRTMGP |
| `amip` | `instantiate` + `resolve` only — pins from `experiments/AMIP/Manifest.toml` |

## What each job does

1. `cd` to the repo workspace
2. Load `climacommon/2025_02_25`
3. Run `setup_amip_env.sh` for the chosen profile
4. Launch `experiments/AMIP/run_simulation.jl` on `--gpus` GPUs via MPI (default 4)

Per date, the script writes:

- `generated/batch_wxquest_progedmf/configs/wxquest_progedmf_<YYYYMMDD>_<HHMM>.yml`
- `generated/batch_wxquest_progedmf/pbs/runner_wxquest_progedmf_<YYYYMMDD>_<HHMM>.pbs`
- Logs: `generated/batch_wxquest_progedmf/logs/ccwx_progedmf_<YYYYMMDD>_<HHMM>.out`

Each generated config updates only `start_date` and appends the date to
`coupler_output_dir`. It does **not** change `era5_initial_condition_dir`.

## Common options

```bash
--dates 19970709 19820615     # explicit dates (overrides DEFAULT_START_DATES)
--dates-file dates.txt        # one date per line (# comments OK)
--base-config PATH            # default: config/subseasonal_configs/wxquest_progedmf.yml
--ic-dir PATH                 # IC validation dir only (does not rewrite YAML)
--skip-ic-check               # skip IC file existence check
--walltime 06:00:00
--nodes 1                     # number of GPU nodes (default 1); ranks = nodes*gpus
--gpus 4                      # GPUs/MPI ranks per node (1-4, default 4)
--climaatmos-path PATH        # Pkg.develop local ClimaAtmos instead of #main
--select "1:ncpus=48:ngpus=4" # optional override; else derived from --gpus
--module climacommon/2025_02_25
--dry-run
--yes                         # skip confirmation prompt
```

## IC directory

The script validates that expected NetCDF prefixes exist per date under `--ic-dir`.
Processed-only trees (e.g. `initial_conditions_0p25deg_dev`) may need `--skip-ic-check`.
Always set `era5_initial_condition_dir` in the base YAML to the same IC tree the run will use.
