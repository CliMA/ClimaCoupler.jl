#!/usr/bin/env python3
"""Submit many short ClimaCoupler forecasts as a pool of PBS workers (ML training data).

``submit_wxquest_batch.py`` submits one PBS job per start date, so every date
pays the full Julia package load + GPU compilation cost. For short forecasts
(e.g. 12 h) that setup dominates the run time. This script instead:

  1. clones the base coupler YAML once per start date (same edits as the batch
     script: ``start_date`` and a per-date ``coupler_output_dir``),
  2. writes one shared config list, and
  3. submits ``--workers`` independent PBS jobs of ``--nodes`` nodes each. Every
     worker runs ``experiments/AMIP/run_simulations_sequential.jl``, which
     claims the next unfinished date from the shared list, runs it, and repeats
     in the same Julia session until nothing is left.

Workers start whenever the scheduler frees enough nodes for one of them, so a
busy machine only delays part of the pool, and a late worker simply finds less
(or no) work left.

Completed runs leave ``artifacts/run_complete.txt`` in their output directory.
Rerunning the same command submits a fresh pool for whatever is left: claims
held by jobs that are no longer queued or running are released first. Failed
runs are not retried unless ``--retry-failed`` is given.

By default every date with a complete set of IC files in ``--ic-dir`` is run.

Example:

  python3 scripts/submit_wxquest_ml_batch.py --dry-run
  python3 scripts/submit_wxquest_ml_batch.py --dates-file dates.txt --yes
"""

import argparse
import datetime as dt
import getpass
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from submit_wxquest_batch import (  # noqa: E402
    ENV_PROFILES,
    check_ic_files,
    confirm,
    format_start_date,
    make_pbs_script,
    parse_start_date,
    read_text,
    update_yaml_config_text,
    write_text,
)

DRIVER = "experiments/AMIP/run_simulations_sequential.jl"
COMPLETE_MARKER = "run_complete.txt"
IC_DATE_RE = re.compile(r"^era5_init_processed_internal_(\d{8})_(\d{4})\.nc$")


def discover_ic_dates(ic_dir: Path) -> list[str]:
    """Start dates (YYYYMMDDHHMM) with a complete set of IC files in ``ic_dir``."""
    dates = []
    for path in sorted(ic_dir.iterdir()):
        m = IC_DATE_RE.match(path.name)
        if m and not check_ic_files(ic_dir, parse_start_date(m.group(1) + m.group(2))):
            dates.append(m.group(1) + m.group(2))
    return dates


def walltime_seconds(walltime: str) -> int:
    h, m, s = (int(x) for x in walltime.split(":"))
    return 3600 * h + 60 * m + s


def active_job_ids() -> set[str]:
    """Numeric IDs of this user's queued or running PBS jobs."""
    res = subprocess.run(
        ["qstat", "-u", getpass.getuser()], check=True, capture_output=True, text=True
    )
    return {
        line.split(".")[0]
        for line in res.stdout.splitlines()
        if re.match(r"^\d+\.", line)
    }


def claim_state(claim: Path, active: set[str]) -> str:
    """One of 'none', 'running', 'failed', 'stale'."""
    if not claim.is_dir():
        return "none"
    owner_file = claim / "owner.txt"
    owner = owner_file.read_text().split()[0].split(".")[0] if owner_file.is_file() else ""
    if owner in active:
        return "running"
    return "failed" if (claim / "failed.txt").is_file() else "stale"


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description=(
            "Run many ClimaCoupler forecasts through a pool of independent PBS "
            "workers that share one work queue, compiling once per worker."
        )
    )
    parser.add_argument(
        "--dates",
        nargs="*",
        default=[],
        help=(
            "Start dates YYYYMMDD or YYYYMMDDHHMM. Default: every complete date in "
            "--ic-dir. Use --dates-file for YYYYMMDD-HHMM (argparse reads a leading "
            "'-' as an option)."
        ),
    )
    parser.add_argument(
        "--dates-file",
        type=str,
        default=None,
        help="File with one start date per line (YYYYMMDD[HH[MM]] or YYYYMMDD-HHMM; # comments OK).",
    )
    parser.add_argument(
        "--base-config",
        type=str,
        default=str(repo_root / "config/weather_configs/weatherbench_progedmf_0m.yml"),
        help="Base coupler YAML cloned per date.",
    )
    parser.add_argument(
        "--ic-dir",
        type=str,
        default=(
            "/glade/campaign/univ/ucit0011/cchristo/wxquest_data/"
            "initial_conditions/initial_conditions_0p25deg_dev/no_hz_stretch"
        ),
        help=(
            "ERA5 IC directory, used to discover dates and validate files. It does "
            "not rewrite era5_initial_condition_dir in the YAML."
        ),
    )
    parser.add_argument(
        "--name",
        type=str,
        default="ic_ml",
        help="Batch label for generated filenames and PBS job names.",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=None,
        help="Where to write configs, list, claims, PBS scripts, and logs. Default: generated/ml_batch_<name>.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of independent PBS worker jobs (default 4; capped at the number of pending dates).",
    )
    parser.add_argument(
        "--nodes",
        type=int,
        default=1,
        help="GPU nodes per worker (default 1; enough for h_elem 30).",
    )
    parser.add_argument("--gpus", type=int, default=4, help="GPUs (MPI ranks) per node (1-4).")
    parser.add_argument("--select", type=str, default=None, help="PBS select override.")
    parser.add_argument(
        "--walltime",
        type=str,
        default="06:00:00",
        help=(
            "PBS walltime per worker (default 6 h; shorter jobs backfill more easily). "
            "Workers stop claiming dates that would not finish in time."
        ),
    )
    parser.add_argument(
        "--est-run-minutes",
        type=int,
        default=40,
        help="Assumed minutes per date until a worker has timed its own runs (default 40).",
    )
    parser.add_argument(
        "--retry-failed",
        action="store_true",
        help="Release claims of failed runs so they are attempted again.",
    )
    parser.add_argument(
        "--env",
        type=str,
        choices=sorted(ENV_PROFILES),
        default="nightly",
        help="Julia environment profile ('amip' = Manifest pins, 'nightly' = main branches).",
    )
    parser.add_argument("--module", type=str, default=None, help="HPC module override.")
    parser.add_argument(
        "--climaatmos-path",
        type=str,
        default="/glade/u/home/cchristo/clima/copies2/ClimaAtmos.jl",
        help="Local ClimaAtmos.jl checkout to use (pass '' to use the Manifest/main version).",
    )
    parser.add_argument("--workspace", type=str, default=str(repo_root), help="ClimaCoupler workspace.")
    parser.add_argument("--skip-ic-check", action="store_true", help="Skip IC file existence check.")
    parser.add_argument("--dry-run", action="store_true", help="Print plan; write/submit nothing.")
    parser.add_argument("--no-submit", action="store_true", help="Write files but do not qsub.")
    parser.add_argument("--yes", action="store_true", help="Skip confirmation prompt.")
    args = parser.parse_args()

    env_profile = ENV_PROFILES[args.env]
    climacommon_module = args.module or env_profile["module"]
    if not 1 <= args.gpus <= 4:
        raise SystemExit(f"--gpus must be between 1 and 4; got {args.gpus}")
    if args.nodes < 1 or args.workers < 1:
        raise SystemExit("--nodes and --workers must be >= 1")
    nranks = args.nodes * args.gpus
    select = args.select or (
        f"{args.nodes}:ncpus={16 * args.gpus}:mpiprocs={args.gpus}:ngpus={args.gpus}"
    )

    climaatmos_path = None
    if args.climaatmos_path:
        climaatmos_path = str(Path(args.climaatmos_path).resolve())
        if not Path(climaatmos_path).is_dir():
            raise SystemExit(f"--climaatmos-path is not a directory: {climaatmos_path}")

    ic_dir = Path(args.ic_dir).resolve()
    workspace = Path(args.workspace).resolve()
    out_base = Path(args.out_dir or repo_root / f"generated/ml_batch_{args.name}").resolve()
    claim_dir = out_base / "claims"
    list_path = out_base / "lists" / f"configs_{args.name}.txt"
    base_config_path = Path(args.base_config).resolve()
    if not base_config_path.is_file():
        raise SystemExit(f"Base config not found: {base_config_path}")
    base_yaml = read_text(base_config_path)
    m_out = re.search(r'^\s*coupler_output_dir:+\s*"([^"]*)"', base_yaml, flags=re.MULTILINE)
    if not m_out:
        raise SystemExit(f"Base config has no quoted coupler_output_dir: {base_config_path}")
    base_out_dir = m_out.group(1).rstrip("/")

    raw_dates: list[str] = []
    if args.dates_file:
        with open(args.dates_file, "r", encoding="utf-8") as f:
            raw_dates = [l.strip() for l in f if l.strip() and not l.strip().startswith("#")]
    elif args.dates:
        raw_dates = list(args.dates)
    else:
        raw_dates = discover_ic_dates(ic_dir)
        if not raw_dates:
            raise SystemExit(f"No complete IC date sets found in {ic_dir}")
    starts = sorted({parse_start_date(d) for d in raw_dates})

    if not args.skip_ic_check:
        missing_any = False
        for sd in starts:
            missing = check_ic_files(ic_dir, sd)
            if missing:
                missing_any = True
                print(f"ERROR: Missing IC files for {format_start_date(sd)}: {', '.join(missing)}")
        if missing_any:
            raise SystemExit("Use --skip-ic-check to proceed anyway.")

    active = active_job_ids()
    runs = []
    for sd in starts:
        start_date_str = format_start_date(sd)
        stem = f"wxquest_progedmf_{sd.strftime('%Y%m%d_%H%M')}"
        output_root = Path(f"{base_out_dir}/{start_date_str}/{stem}")
        claim = claim_dir / stem
        if (output_root / "artifacts" / COMPLETE_MARKER).is_file():
            state = "complete"
        else:
            state = claim_state(claim, active)
        runs.append({
            "start_date_str": start_date_str,
            "config_path": out_base / "configs" / f"{stem}.yml",
            "claim": claim,
            "state": state,
        })

    releasable = {"stale", "failed"} if args.retry_failed else {"stale"}
    to_release = [r for r in runs if r["state"] in releasable]
    n_pending = sum(r["state"] in {"none"} | releasable for r in runs)
    n_workers = min(args.workers, n_pending)

    def dates_in(*states: str) -> str:
        return " ".join(r["start_date_str"] for r in runs if r["state"] in states) or "-"

    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    workers = []
    for w in range(n_workers):
        tag = f"{args.name}_w{w:02d}"
        workers.append({
            "job_name": f"ccml_{tag}"[:60],
            "pbs_path": out_base / "pbs" / f"runner_{tag}.pbs",
            "log_path": out_base / "logs" / f"ccml_{tag}_{stamp}.out",
        })

    print(f"Base config:  {base_config_path}")
    print(f"IC directory: {ic_dir}")
    print(f"Output root:  {base_out_dir}/<start_date>/")
    print(f"Generated in: {out_base}")
    print(f"Env profile:  {args.env} ({env_profile['description']}), module {climacommon_module}")
    print(f"ClimaAtmos:   {climaatmos_path or '(Manifest / #main)'}")
    print(f"Workers:      {n_workers} x ({args.nodes} node(s) x {args.gpus} GPU/node = {nranks} ranks, "
          f"select={select}), walltime {args.walltime}")
    print(f"Dates:        {len(runs)} total")
    print(f"  pending:    {dates_in('none')}")
    print(f"  complete:   {dates_in('complete')}")
    print(f"  running:    {dates_in('running')}")
    print(f"  failed:     {dates_in('failed')}{' (will retry)' if args.retry_failed else ''}")
    print(f"  stale:      {dates_in('stale')} (claim released, will rerun)")
    print()

    if args.dry_run:
        return
    if n_workers == 0:
        print("Nothing left to run.")
        return
    if not args.yes:
        action = "write files" if args.no_submit else "write files and submit"
        if not confirm(f"Proceed to {action} {n_workers} worker(s)?"):
            print("Aborted by user.")
            return

    for run in runs:
        write_text(run["config_path"], update_yaml_config_text(base_yaml, run["start_date_str"]))
    write_text(list_path, "".join(f"{r['config_path']}\n" for r in runs))
    claim_dir.mkdir(parents=True, exist_ok=True)
    for run in to_release:
        shutil.rmtree(run["claim"])

    driver_args = (
        f"--claim_dir {claim_dir} "
        f"--walltime_seconds {walltime_seconds(args.walltime)} "
        f"--est_run_seconds {60 * args.est_run_minutes}"
    )
    for worker in workers:
        pbs_text = make_pbs_script(
            workspace=workspace,
            climacommon_module=climacommon_module,
            env_mode=args.env,
            config_path=list_path,
            job_name=worker["job_name"],
            log_path=worker["log_path"],
            walltime=args.walltime,
            select=select,
            nodes=args.nodes,
            gpus_per_node=args.gpus,
            climaatmos_path=climaatmos_path,
            driver=DRIVER,
            config_flag="--config_list",
            driver_args=driver_args,
        )
        write_text(worker["pbs_path"], pbs_text)
        worker["log_path"].parent.mkdir(parents=True, exist_ok=True)

        if args.no_submit:
            print(f"Wrote {worker['job_name']} (not submitted): {worker['pbs_path']}")
            continue
        try:
            res = subprocess.run(
                ["qsub", str(worker["pbs_path"])], check=True, capture_output=True, text=True
            )
            print(f"Submitted {worker['job_name']} -> {res.stdout.strip()}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to submit {worker['pbs_path']}: {e.stderr.strip()}")


if __name__ == "__main__":
    main()
