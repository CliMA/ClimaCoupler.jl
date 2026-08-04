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
    echo "--- nightly: tracking main on Buildkite nightly upstream packages"
    NIGHTLY_PKGS=(ClimaAtmos ClimaCore ClimaCoreMakie ClimaTimeSteppers Thermodynamics ClimaLand SurfaceFluxes RRTMGP)

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

    pkgs_julia=$(printf 'Pkg.PackageSpec(; name="%s", rev="main"), ' "${NIGHTLY_PKGS[@]}")
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
        echo "--- developing local ClimaAtmos from: $CLIMAATMOS_PATH"
        julia --project="$AMIP_PATH" -e \
            "using Pkg; Pkg.develop(Pkg.PackageSpec(; name=\"ClimaAtmos\", path=\"$CLIMAATMOS_PATH\"))"
    fi
else
    echo "Unknown mode: $MODE (expected 'amip' or 'nightly')" >&2
    return 1 2>/dev/null || exit 1
fi

julia --project="$AMIP_PATH" -e 'using Pkg; Pkg.resolve()'
julia --project="$AMIP_PATH" -e 'using Pkg; Pkg.precompile()'
julia --project="$AMIP_PATH" -e 'using Pkg; Pkg.status()'

echo "--- AMIP env ready (mode=$MODE, project=$AMIP_PATH)"
