#!/usr/bin/env python3

import argparse
import datetime as dt
import re
import subprocess
from pathlib import Path
from textwrap import dedent
from typing import Optional


DEFAULT_START_DATES = [
    "19970709",
    "19820615",
    "20150627",
    "20230711",
    "19870731",
    "20020617",
    "20150514",
    "20140622",
    "20020731",
    "19930601",
]

# Files expected per initialization date in the IC directory
EXPECTED_IC_PREFIXES = [
    "albedo_processed",
    "era5_bucket_processed",
    "era5_init_processed_internal",
    "era5_land",
    "era5_land_processed",
    "era5_model_levels",
    "era5_raw",
    "era5_surface",
    "sic_processed",
    "sst_processed",
    "surf_processed",
]


def parse_start_date(raw: str) -> dt.datetime:
    raw = raw.strip()
    formats = ["%Y%m%d-%H%M", "%Y%m%d%H%M", "%Y%m%d%H", "%Y%m%d"]
    for fmt in formats:
        try:
            return dt.datetime.strptime(raw, fmt)
        except ValueError:
            continue
    raise ValueError(
        f"Unrecognized date format: {raw}. Expected YYYYMMDD, YYYYMMDDHHMM, or YYYYMMDD-HHMM"
    )


def format_start_date(start: dt.datetime) -> str:
    """Coupler/Atmos config string: YYYYMMDD at 00Z, else YYYYMMDD-HHMM."""
    if start.hour == 0 and start.minute == 0:
        return start.strftime("%Y%m%d")
    return start.strftime("%Y%m%d-%H%M")


def read_text(path: Path) -> str:
    with path.open("r", encoding="utf-8") as f:
        return f.read()


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write(content)


def check_ic_files(ic_dir: Path, start: dt.datetime) -> list[str]:
    """Return list of missing IC files for a given date."""
    ymd = start.strftime("%Y%m%d")
    hhmm = start.strftime("%H%M")
    missing = []
    for prefix in EXPECTED_IC_PREFIXES:
        fname = f"{prefix}_{ymd}_{hhmm}.nc"
        if not (ic_dir / fname).is_file():
            missing.append(fname)
    return missing


def update_yaml_config_text(base_text: str, start_date_str: str) -> str:
    """Update start_date and append date suffix to coupler_output_dir."""
    text = base_text
    # Match coupler_output_dir (single or double colon for robustness)
    m = re.search(r'(^\s*coupler_output_dir:+\s*)"([^"]*)"', text, flags=re.MULTILINE)
    if m:
        base_out_dir = m.group(2)
        new_out_dir = f"{base_out_dir.rstrip('/')}/{start_date_str}"
        text = re.sub(
            r'(^\s*coupler_output_dir:+\s*).*$',
            rf'\1"{new_out_dir}"',
            text,
            flags=re.MULTILINE,
        )
    text = re.sub(r'(^\s*start_date:\s*).*$', rf'\1"{start_date_str}"', text, flags=re.MULTILINE)
    return text


ENV_PROFILES = {
    "amip": {
        "module": "climacommon/2025_02_25",
        "description": "Manifest-pinned experiments/AMIP environment",
    },
    "nightly": {
        "module": "climacommon/2025_02_25",
        "description": "Main branches of key upstream packages (Buildkite nightly)",
    },
}


def make_pbs_script(
    workspace: Path,
    climacommon_module: str,
    env_mode: str,
    config_path: Path,
    job_name: str,
    log_path: Path,
    walltime: str,
    select: str,
    nodes: int,
    gpus_per_node: int,
    climaatmos_path: Optional[str] = None,
) -> str:
    nranks = nodes * gpus_per_node
    setup_script = workspace / "scripts/setup_amip_env.sh"
    climaatmos_export = (
        f'export CLIMAATMOS_PATH="{climaatmos_path}"\n        '
        if climaatmos_path
        else ""
    )
    return dedent(
        f"""\
        #!/bin/bash
        #PBS -N {job_name}
        #PBS -j oe
        #PBS -A UCIT0012
        #PBS -q main
        # #PBS -l job_priority=premium
        #PBS -o {log_path}
        #PBS -l walltime={walltime}
        #PBS -l select={select}
        
        cd {workspace}
        
        AMIP_PATH='experiments/AMIP/'
        DRIVER='experiments/AMIP/run_simulation.jl'
        CONFIG_FILE='{config_path}'
        ENV_MODE='{env_mode}'
        
        export MODULEPATH="/glade/campaign/univ/ucit0011/ClimaModules-Derecho:$MODULEPATH"
        unset LD_LIBRARY_PATH
        module load {climacommon_module}
        module list
        
        export CLIMACOMMS_DEVICE="CUDA"
        export CLIMACOMMS_CONTEXT="MPI"

        # Cray MPICH: the SMP-aware MPI_Reduce algorithm (MPIR_Reduce_intra_smp)
        # fails with MPI_ERR_TRUNCATE on large GPU-buffer reduces from the
        # ClimaDiagnostics NetCDF remap (all ranks on one node). Force the
        # non-SMP reduce algorithm and disable the optimized reduce collective.
        export MPICH_REDUCE_NO_SMP=1
        export MPICH_COLL_OPT_OFF=mpi_reduce
        export MPICH_GPU_IPC_ENABLED=0
        export MPICH_SMP_SINGLE_COPY_MODE=NONE
        
        if [ -n "$SCRATCH" ]; then
            export TMPDIR=${{SCRATCH}}/temp
        else
            export TMPDIR=/tmp/clima_${{USER}}
        fi
        mkdir -p ${{TMPDIR}}
        
        export AMIP_PATH
        {climaatmos_export}source {setup_script} "$ENV_MODE"
        
        CMDRN="julia --project=$AMIP_PATH $DRIVER --config_file $CONFIG_FILE"
        "$MPITRAMPOLINE_MPIEXEC" -np {nranks} -ppn {gpus_per_node} set_gpu_rank "$CMDRN"
        """
    ).strip() + "\n"


def confirm(prompt: str) -> bool:
    try:
        ans = input(f"{prompt} [y/N]: ").strip().lower()
    except EOFError:
        return False
    return ans in {"y", "yes"}


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Submit multiple ClimaCoupler wxquest_progedmf jobs on Derecho, "
            "one per start date, by generating per-date configs and PBS scripts."
        )
    )
    repo_root = Path(__file__).resolve().parent.parent
    parser.add_argument(
        "--dates",
        nargs="*",
        default=[],
        help=(
            "Start dates as YYYYMMDD[HH[MM]] (space-separated). "
            "If omitted, uses the built-in default list in the script."
        ),
    )
    parser.add_argument(
        "--dates-file",
        type=str,
        default=None,
        help="File with one start date per line (YYYYMMDD[HH[MM]]).",
    )
    parser.add_argument(
        "--base-config",
        type=str,
        default=str(repo_root / "config/subseasonal_configs/wxquest_progedmf.yml"),
        help="Path to the base YAML config to clone and modify per job",
    )
    parser.add_argument(
        "--ic-dir",
        type=str,
        default="/glade/campaign/univ/ucit0011/cchristo/wxquest_data/initial_conditions/initial_conditions_v1_dev",
        help="Directory containing the ERA5 IC files (used for validation)",
    )
    parser.add_argument(
        "--workspace",
        type=str,
        default=str(repo_root),
        help="ClimaCoupler workspace directory (where runner scripts expect to cd)",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=str(repo_root / "generated/batch_wxquest_progedmf"),
        help="Output directory to write generated configs and PBS scripts",
    )
    parser.add_argument(
        "--walltime",
        type=str,
        default="06:00:00",
        help="PBS walltime per job",
    )
    parser.add_argument(
        "--nodes",
        type=int,
        default=1,
        help=(
            "Number of GPU nodes to request (default: 1). Total MPI ranks = "
            "nodes * gpus-per-node. Use >1 for jobs too large for a single "
            "node (e.g. high h_elem that OOMs on 4 GPUs)."
        ),
    )
    parser.add_argument(
        "--gpus",
        type=int,
        default=4,
        help=(
            "GPUs (MPI ranks) per node (default: 4). Derecho GPU nodes have "
            "4 A100s, so max is 4. Combined with --nodes to derive the PBS "
            "select string and mpiexec -np/-ppn."
        ),
    )
    parser.add_argument(
        "--select",
        type=str,
        default=None,
        help=(
            "PBS select resource string override. If omitted, derived from "
            "--nodes/--gpus as "
            "'<nodes>:ncpus=<16*gpus>:mpiprocs=<gpus>:ngpus=<gpus>'."
        ),
    )
    parser.add_argument(
        "--env",
        type=str,
        choices=sorted(ENV_PROFILES),
        default="nightly",
        help=(
            "Julia environment profile: 'nightly' tracks main on key upstream packages; "
            "'amip' uses Manifest.toml pins only"
        ),
    )
    parser.add_argument(
        "--module",
        type=str,
        default=None,
        help=(
            "Module to load for the job environment. "
            "Defaults to climacommon/2025_02_25 for both profiles."
        ),
    )
    parser.add_argument(
        "--climaatmos-path",
        type=str,
        default=None,
        help=(
            "Local ClimaAtmos.jl checkout to Pkg.develop instead of ClimaAtmos#main "
            "(exported as CLIMAATMOS_PATH into the PBS job)."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Only print planned actions; do not write files or submit jobs",
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Do not prompt for confirmation; proceed to write/submit",
    )
    parser.add_argument(
        "--skip-ic-check",
        action="store_true",
        help="Skip existence check for IC files in the IC directory",
    )

    args = parser.parse_args()

    env_profile = ENV_PROFILES[args.env]
    climacommon_module = args.module or env_profile["module"]

    if args.gpus < 1 or args.gpus > 4:
        raise SystemExit(
            f"--gpus must be between 1 and 4 (Derecho GPU nodes have 4 A100s); got {args.gpus}"
        )
    if args.nodes < 1:
        raise SystemExit(f"--nodes must be >= 1; got {args.nodes}")
    nranks = args.nodes * args.gpus
    select = args.select or (
        f"{args.nodes}:ncpus={16 * args.gpus}:mpiprocs={args.gpus}:ngpus={args.gpus}"
    )

    climaatmos_path = None
    if args.climaatmos_path:
        climaatmos_path = str(Path(args.climaatmos_path).resolve())
        if not Path(climaatmos_path).is_dir():
            raise SystemExit(f"--climaatmos-path is not a directory: {climaatmos_path}")

    dates: list[str] = []
    if args.dates_file:
        file_path = Path(args.dates_file)
        if not file_path.is_file():
            raise SystemExit(f"Dates file not found: {file_path}")
        with file_path.open("r", encoding="utf-8") as f:
            for line in f:
                raw = line.strip()
                if raw and not raw.startswith("#"):
                    dates.append(raw)
    elif args.dates:
        dates.extend(args.dates)
    else:
        dates = list(DEFAULT_START_DATES)

    parsed_dates: list[dt.datetime] = []
    for raw in dates:
        parsed_dates.append(parse_start_date(raw))

    base_config_path = Path(args.base_config).resolve()
    if not base_config_path.is_file():
        raise SystemExit(f"Base config not found: {base_config_path}")

    ic_dir = Path(args.ic_dir).resolve()
    workspace = Path(args.workspace).resolve()
    out_base = Path(args.out_dir).resolve()

    base_yaml = read_text(base_config_path)
    m_out = re.search(r'(^\s*coupler_output_dir:+\s*)"([^"]*)"', base_yaml, flags=re.MULTILINE)
    base_out_dir = m_out.group(2) if m_out else None

    planned = []
    for sd in parsed_dates:
        ymd = sd.strftime("%Y%m%d")
        hhmm = sd.strftime("%H%M")
        start_date_str = format_start_date(sd)
        out_tag = start_date_str.replace("-", "_")
        output_dir = f"{base_out_dir.rstrip('/')}/{out_tag}" if base_out_dir else None

        config_out = out_base / "configs" / f"wxquest_progedmf_{ymd}_{hhmm}.yml"
        pbs_out = out_base / "pbs" / f"runner_wxquest_progedmf_{ymd}_{hhmm}.pbs"
        log_out = out_base / "logs" / f"ccwx_progedmf_{ymd}_{hhmm}.out"
        job_name = f"cc_wx_prog_{ymd}" if hhmm == "0000" else f"cc_wx_prog_{ymd}_{hhmm}"

        planned.append({
            "start": sd,
            "start_date_str": start_date_str,
            "coupler_output_dir": output_dir,
            "config_path": str(config_out),
            "pbs_path": str(pbs_out),
            "log_path": str(log_out),
            "job_name": job_name,
        })

    print(f"IC directory: {ic_dir}")
    print(f"Base config:  {base_config_path}")
    print(f"Output base:  {out_base}")
    print(f"Env profile:  {args.env} ({env_profile['description']})")
    print(f"HPC module:   {climacommon_module}")
    print(f"Resources:    {args.nodes} node(s) x {args.gpus} GPU/node = {nranks} ranks (select={select})")
    print(f"ClimaAtmos:   {climaatmos_path or '(nightly #main / Manifest)'}")
    print(f"\nPlanned submissions ({len(planned)} jobs):")
    for item in planned:
        print(
            f"  {item['job_name']:20s}  start_date={item['start_date_str']}  "
            f"output_dir={item['coupler_output_dir'] or 'N/A'}"
        )
    print()

    if args.dry_run:
        return

    # Validate IC files
    if not args.skip_ic_check:
        all_ok = True
        for item in planned:
            missing = check_ic_files(ic_dir, item["start"])
            if missing:
                all_ok = False
                print(f"ERROR: Missing IC files for {item['start_date_str']}:")
                for f in missing:
                    print(f"  - {ic_dir / f}")
        if not all_ok:
            print("\nUse --skip-ic-check to proceed anyway.")
            raise SystemExit(1)
        print(f"IC check passed: all {len(EXPECTED_IC_PREFIXES)} file types present for all {len(planned)} dates.\n")

    if not args.yes:
        if not confirm("Proceed to write files and submit these jobs?"):
            print("Aborted by user.")
            return

    # Materialize files and submit
    for item in planned:
        cfg_text = update_yaml_config_text(base_yaml, item["start_date_str"])
        write_text(Path(item["config_path"]), cfg_text)

        pbs_text = make_pbs_script(
            workspace=workspace,
            climacommon_module=climacommon_module,
            env_mode=args.env,
            config_path=Path(item["config_path"]).resolve(),
            job_name=item["job_name"],
            log_path=Path(item["log_path"]).resolve(),
            walltime=args.walltime,
            select=select,
            nodes=args.nodes,
            gpus_per_node=args.gpus,
            climaatmos_path=climaatmos_path,
        )
        write_text(Path(item["pbs_path"]), pbs_text)

        Path(item["log_path"]).parent.mkdir(parents=True, exist_ok=True)

        try:
            res = subprocess.run(["qsub", item["pbs_path"]], check=True, capture_output=True, text=True)
            job_id = res.stdout.strip()
            print(f"Submitted {item['job_name']} -> {job_id}")
        except subprocess.CalledProcessError as e:
            print(f"Failed to submit {item['pbs_path']}: {e.stderr.strip()}")


if __name__ == "__main__":
    main()
