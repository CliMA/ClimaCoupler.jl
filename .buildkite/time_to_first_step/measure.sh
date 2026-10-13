#!/bin/bash
# Measure time to first step of the flagship AMIP simulation.
# Usage: measure.sh JOB_ID [--main]
# Extra arguments are passed to .buildkite/add_upstream_packages.jl.
# Each job uses its own copy of experiments/AMIP and its own depot, since
# jobs in a build share a checkout and depot.
set -euo pipefail

job_id=$1
shift

project="experiments/AMIP_${job_id}"
rm -rf "$project"
cp -r experiments/AMIP "$project"
export JULIA_DEPOT_PATH="${CI_BUILD_DIR:-$PWD}/depot/${job_id}:${JULIA_DEPOT_PATH:-}"

echo "--- Add upstream packages"
julia --project="$project" .buildkite/add_upstream_packages.jl "$@"

echo "--- Instantiate AMIP env"
julia --project="$project" -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project="$project" -e 'using Pkg; Pkg.precompile()'
julia --project="$project" -e 'using Pkg; Pkg.status()'

echo "--- Measure time to first step"
srun --cpu-bind=threads --cpus-per-task=4 julia --threads=3 --color=yes \
    --project="$project" experiments/AMIP/time_to_first_step.jl \
    --config_file config/amip_configs/amip.yml --job_id "$job_id"
