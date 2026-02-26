#!/bin/bash
#
# Calibration workflow for subseasonal_weekly pipeline.
# Run this from tmux on the login node.
#
# Usage: ./run_full_calibration.sh
#
# This script:
# 1. Ensures package dependencies are installed
# 2. Runs the calibration (ClimaGPUBackend handles GPU worker submission)
#

set -e  # Exit on error

cd /glade/u/home/zhaoyi/weekly_calibration/ClimaCoupler.jl

# Log to file AND terminal using tee
LOGFILE="calibration_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1
echo "Logging to: $LOGFILE"

# Load modules
export MODULEPATH="/glade/campaign/univ/ucit0011/ClimaModules-Derecho:$MODULEPATH"
module purge
module load climacommon/2025_02_25

echo "=============================================="
echo "  Subseasonal Weekly Calibration"
echo "=============================================="
echo ""

# Step 1: Ensure packages are instantiated
echo "[1/3] Ensuring package dependencies are installed..."
julia --project=experiments/ClimaEarth -e 'using Pkg; Pkg.instantiate()'

# Step 2: Generate observations (creates obs_vec.jld2 and norm_stats.jld2)
# Uses CERES for radiation variables (rsut, rlut, etc.), ERA5 for others (tas, etc.)
echo "[2/3] Generating observations..."
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal_weekly/generate_observations.jl

# Step 3: Run calibration (DerechoBackend handles GPU worker submission)
echo "[3/3] Starting calibration..."
echo ""

julia --project=experiments/ClimaEarth experiments/calibration/subseasonal_weekly/run_calibration.jl

echo ""
echo "=============================================="
echo "  Calibration Complete!"
echo "=============================================="
