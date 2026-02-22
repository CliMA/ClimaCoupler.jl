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

cd /glade/u/home/cchristo/clima/copies3/ClimaCoupler.jl

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
echo "[1/2] Ensuring package dependencies are installed..."
julia --project=experiments/ClimaEarth -e 'using Pkg; Pkg.instantiate()'

# Step 2: Run calibration (ClimaGPUBackend handles worker submission)
echo "[2/2] Starting calibration..."
echo ""

julia --project=experiments/ClimaEarth experiments/calibration/subseasonal_weekly/run_calibration.jl

echo ""
echo "=============================================="
echo "  Calibration Complete!"
echo "=============================================="
