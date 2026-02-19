#!/bin/bash
#
# Full calibration workflow for TransformInversion
# Run this from tmux on the login node!
#
# Usage: ./run_full_calibration.sh
#
# This script:
# 1. Submits precompute job to cpudev queue (fast, ~2 min)
# 2. Waits for it to complete
# 3. Runs the main calibration (which spawns GPU workers)
#

set -e  # Exit on error

cd /glade/u/home/zhaoyi/ClimaCoupler.jl

# Log to file AND terminal using tee
LOGFILE="calibration_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOGFILE") 2>&1
echo "Logging to: $LOGFILE"

# Load modules
export MODULEPATH="/glade/campaign/univ/ucit0011/ClimaModules-Derecho:$MODULEPATH"
module purge
module load climacommon/2025_02_25

echo "=============================================="
echo "  Full Calibration Workflow"
echo "=============================================="
echo ""

# Step 1: Submit precompute job to develop queue
echo "[1/4] Submitting precompute job to develop queue..."
JOB_ID=$(qsub experiments/calibration/subseasonal/precompute.pbs)
echo "      Job submitted: $JOB_ID"

# Extract just the job number
JOB_NUM=$(echo $JOB_ID | cut -d. -f1)

# Step 2: Wait for precompute to complete
echo "[2/4] Waiting for precompute job to complete..."
echo "      (develop queue is usually fast - should take ~2-5 minutes)"

while true; do
    # Check if job is still in queue (qstat returns non-zero if job doesn't exist)
    if ! qstat $JOB_NUM &>/dev/null; then
        echo "      Precompute job finished!"
        break
    fi
    
    # Get job status from qstat output
    STATUS=$(qstat $JOB_NUM 2>/dev/null | grep $JOB_NUM | awk '{print $5}')
    
    case "$STATUS" in
        R) echo "      Job is running..." ;;
        Q) echo "      Job is queued..." ;;
        F|E|"") 
            echo "      Precompute job finished!"
            break 
            ;;
        *) echo "      Job status: $STATUS" ;;
    esac
    
    sleep 10
done

# Wait for Lustre filesystem to sync
sleep 10

# Check if precompute succeeded
if [ -f "experiments/calibration/subseasonal/ekp_inputs.jld2" ]; then
    echo "      ✓ ekp_inputs.jld2 created successfully"
else
    echo "      ✗ ERROR: ekp_inputs.jld2 not found!"
    echo "      Check precompute.log for errors"
    exit 1
fi

# Step 3: Ensure packages are instantiated
echo "[3/3] Ensuring package dependencies are installed..."
julia --project=experiments/ClimaEarth -e 'using Pkg; Pkg.instantiate()'

# Step 4: Run calibration
echo "[4/4] Starting main calibration..."
echo "      (This will spawn GPU worker jobs)"
echo ""

julia --project=experiments/ClimaEarth experiments/calibration/subseasonal/run_calibration.jl

echo ""
echo "=============================================="
echo "  Calibration Complete!"
echo "=============================================="
