#!/usr/bin/env bash
# Set up the experiments/AMIP Julia environment for batch or interactive runs.
#
# Usage (from repo root):
#   source scripts/setup_amip_env.sh amip      # manifest-pinned deps
#   source scripts/setup_amip_env.sh nightly   # main branches (Buildkite nightly style)
#
# Optional:
#   export CLIMAATMOS_PATH=/path/to/ClimaAtmos.jl
#   # when set, Pkg.develop that checkout instead of ClimaAtmos#main
#
# Requires AMIP_PATH (defaults to experiments/AMIP/).

set -euo pipefail

MODE="${1:-amip}"
AMIP_PATH="${AMIP_PATH:-experiments/AMIP/}"

julia --project="$AMIP_PATH" -e 'using Pkg; Pkg.instantiate(;verbose=true)'

if [[ "$MODE" == "nightly" ]]; then
    # Match .buildkite/nightly/pipeline.yml UPSTREAM_PACKAGES default.
    # `Name@version` pins a release; otherwise track main.
    echo "--- nightly: tracking main on Buildkite nightly upstream packages"
    NIGHTLY_PKGS=(ClimaAtmos ClimaCore ClimaTimeSteppers Thermodynamics ClimaLand SurfaceFluxes RRTMGP CloudMicrophysics)

    if [[ -n "${CLIMAATMOS_PATH:-}" ]]; then
        if [[ ! -d "$CLIMAATMOS_PATH" ]]; then
            echo "CLIMAATMOS_PATH is set but not a directory: $CLIMAATMOS_PATH" >&2
            return 1 2>/dev/null || exit 1
        fi
        echo "--- developing local ClimaAtmos from: $CLIMAATMOS_PATH"
        # Drop ClimaAtmos from the #main list (do not leave an empty element).
        filtered=()
        for pkg in "${NIGHTLY_PKGS[@]}"; do
            [[ "$pkg" == "ClimaAtmos" ]] && continue
            filtered+=("$pkg")
        done
        NIGHTLY_PKGS=("${filtered[@]}")
    fi

    pkgs_julia=""
    for entry in "${NIGHTLY_PKGS[@]}"; do
        if [[ "$entry" == *@* ]]; then
            pkgs_julia+="Pkg.PackageSpec(; name=\"${entry%@*}\", version=\"${entry#*@}\"), "
        else
            pkgs_julia+="Pkg.PackageSpec(; name=\"$entry\", rev=\"main\"), "
        fi
    done
    julia --project="$AMIP_PATH" -e "using Pkg; Pkg.add([${pkgs_julia}])"

    if [[ -n "${CLIMAATMOS_PATH:-}" ]]; then
        julia --project="$AMIP_PATH" -e \
            "using Pkg; Pkg.develop(Pkg.PackageSpec(; name=\"ClimaAtmos\", path=\"$CLIMAATMOS_PATH\"))"
    fi
elif [[ "$MODE" == "amip" ]]; then
    echo "--- amip: using Manifest.toml pins (instantiate only)"
    if [[ -n "${CLIMAATMOS_PATH:-}" ]]; then
        if [[ ! -d "$CLIMAATMOS_PATH" ]]; then
            echo "CLIMAATMOS_PATH is set but not a directory: $CLIMAATMOS_PATH" >&2
            return 1 2>/dev/null || exit 1
        fi
        # Avoid Pkg.develop/resolve: Manifest pins (ClimaLand 1.12, GeometryOpsCore
        # 0.1.12, ...) may be newer than what local registries currently list.
        # Point ClimaAtmos at the local checkout by rewriting Manifest (+ Project sources).
        CLIMAATMOS_PATH="$(cd "$CLIMAATMOS_PATH" && pwd)"
        echo "--- pinning local ClimaAtmos in Manifest from: $CLIMAATMOS_PATH"
        AMIP_PATH_ABS="$(cd "$AMIP_PATH" && pwd)"
        python3 - "$AMIP_PATH_ABS" "$CLIMAATMOS_PATH" <<'PY'
import re
import sys
from pathlib import Path

amip = Path(sys.argv[1])
atmos = Path(sys.argv[2]).resolve()
project = amip / "Project.toml"

# Julia 1.11+ prefers Manifest-v1.11.toml; also patch Manifest.toml if present.
manifests = [p for p in (amip / "Manifest-v1.11.toml", amip / "Manifest.toml") if p.is_file()]
if not manifests:
    raise SystemExit(f"No Manifest.toml found under {amip}")

def pin_manifest(manifest: Path) -> None:
    text = manifest.read_text()
    out, in_atmos, replaced = [], False, False
    for line in text.splitlines(keepends=True):
        if line.startswith("[[deps.ClimaAtmos]]"):
            in_atmos = True
            out.append(line)
            continue
        if in_atmos and line.startswith("[["):
            in_atmos = False
        if in_atmos and (
            line.startswith("git-tree-sha1")
            or line.startswith("repo-url")
            or line.startswith("repo-rev")
            or line.startswith("path =")
        ):
            if not replaced:
                out.append(f'path = "{atmos}"\n')
                replaced = True
            continue
        out.append(line)
    if not replaced:
        raise SystemExit(f"Could not find ClimaAtmos entry to pin in {manifest}")
    manifest.write_text("".join(out))
    print(f"Pinned ClimaAtmos in {manifest.name} -> {atmos}")

for m in manifests:
    pin_manifest(m)

# Keep [deps] ClimaAtmos as a UUID string; only record the path under [sources].
CLIMAATMOS_UUID = "b2c96348-7fb7-4fe0-8da9-78d88439e717"
proj = project.read_text()
proj = re.sub(
    r'^ClimaAtmos\s*=\s*\{path\s*=\s*"[^"]*"\}\s*$',
    f'ClimaAtmos = "{CLIMAATMOS_UUID}"',
    proj,
    count=1,
    flags=re.M,
)
src_line = f'ClimaAtmos = {{path = "{atmos}"}}'
if re.search(r"^\[sources\]", proj, flags=re.M):
    # Replace only inside the [sources] section.
    def _sources_repl(match: re.Match[str]) -> str:
        block = match.group(0)
        if re.search(r"^ClimaAtmos\s*=", block, flags=re.M):
            block = re.sub(
                r"^ClimaAtmos\s*=\s*.*$",
                src_line,
                block,
                count=1,
                flags=re.M,
            )
        else:
            block = block.rstrip() + "\n" + src_line + "\n"
        return block

    proj = re.sub(
        r"^\[sources\][^\[]*",
        _sources_repl,
        proj,
        count=1,
        flags=re.M,
    )
else:
    if not proj.endswith("\n"):
        proj += "\n"
    proj += "\n[sources]\n" + src_line + "\n"
project.write_text(proj)
print(f"Project [sources] ClimaAtmos -> {atmos}")
PY
    fi
else
    echo "Unknown mode: $MODE (expected 'amip' or 'nightly')" >&2
    return 1 2>/dev/null || exit 1
fi

# Nightly may have changed many pins; amip keeps Manifest pins and should not
# re-resolve (ClimaLand 1.12+ may be Manifest-only vs current registries).
if [[ "$MODE" == "nightly" ]]; then
    julia --project="$AMIP_PATH" -e 'using Pkg; Pkg.resolve()'
fi
julia --project="$AMIP_PATH" -e 'using Pkg; Pkg.precompile()'
julia --project="$AMIP_PATH" -e 'using Pkg; Pkg.status()'

echo "--- AMIP env ready (mode=$MODE, project=$AMIP_PATH)"
