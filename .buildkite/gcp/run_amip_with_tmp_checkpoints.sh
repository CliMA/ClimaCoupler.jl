#!/bin/bash
# Wrapper script to run AMIP simulation with output written to /tmp first,
# then copied back to avoid I/O issues on GCP.

set -e

CONFIG_FILE="$1"
JOB_ID="$2"
COUPLER_OUTPUT_DIR="experiments/ClimaEarth/output/"
mkdir -p $COUPLER_OUTPUT_DIR
TMP_OUTPUT_DIR="/tmp/clima_coupler_output_${SLURM_JOB_ID}"
echo "coupler_output_dir: $TMP_OUTPUT_DIR" >> "$CONFIG_FILE"

# Function to copy output directory from /tmp to final location
copy_output_dir() {
    if [ -d "$TMP_OUTPUT_DIR" ] && [ -n "$(ls -A "$TMP_OUTPUT_DIR" 2>/dev/null)" ]; then
        echo "--- Copying output from $TMP_OUTPUT_DIR to $COUPLER_OUTPUT_DIR"
        cp -r "$TMP_OUTPUT_DIR"/* "$COUPLER_OUTPUT_DIR" 2>/dev/null || true
        echo "Successfully copied all files."
    fi
}

# Set up cleanup on exit
trap copy_output_dir EXIT INT TERM

echo "--- Running AMIP simulation with temporary output directory: $TMP_OUTPUT_DIR"
julia --threads=3 --color=yes --project=experiments/ClimaEarth/ \
    experiments/ClimaEarth/run_amip.jl \
    --config_file "$CONFIG_FILE" \
    --job_id "$JOB_ID"
