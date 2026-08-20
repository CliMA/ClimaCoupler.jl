#!/bin/bash
#PBS -N calibration_smoke
#PBS -A UCIT0011
#PBS -q develop
#PBS -l select=1:ncpus=8:ngpus=1:mem=40GB
#PBS -l walltime=01:00:00
#PBS -j oe

# PBS wrapper for smoke_test.jl - one ensemble member at the prior centre,
# short window, develop queue. Run this AFTER a green preflight and BEFORE
# submitting any ensemble driver. See smoke_test.jl for why the gate exists.
#
# Submit with the config and output paths in the environment, so this file
# needs no editing between calibrations:
#
#   qsub -v CALIBRATION_CONFIG=$REPO/experiments/calibration/amip/config/rlut_pigroups.jl,\
#SMOKE_DIR=/glade/derecho/scratch/$USER/pigroups_smoke,SMOKE_DAYS=3 \
#        -o /glade/derecho/scratch/$USER/pigroups_smoke/smoke.log \
#        $REPO/experiments/calibration/amip/smoke_test.sh
#
# PASS is the line "SMOKE TEST PASSED" in the log AND exit code 0. Anything
# else - "Found NaN", CUDA.KernelException, a nonzero exit - is a stop.

source /etc/profile.d/modules.sh 2>/dev/null
module purge 2>/dev/null
module load climacommon/2025_02_25

REPO=${REPO:-/glade/u/home/$USER/worktree/ClimaCoupler/ne/amip-calibration}
JULIA=${CALIBRATION_JULIA:-/glade/campaign/univ/ucit0011/software/julia/julia-1.11.3/bin/julia}

: "${CALIBRATION_CONFIG:?set CALIBRATION_CONFIG}"
: "${SMOKE_DIR:?set SMOKE_DIR}"
export CALIBRATION_CONFIG SMOKE_DIR
export SMOKE_DAYS=${SMOKE_DAYS:-3}
export SMOKE_SIM_SECONDS=${SMOKE_SIM_SECONDS:-0}   # >0: construction-only window of this many seconds
export CLIMACOMMS_CONTEXT=SINGLETON
export CLIMACOMMS_DEVICE=CUDA

mkdir -p "$SMOKE_DIR"
cd "$REPO"

echo "=== $(date) smoke test"
echo "    config    : $CALIBRATION_CONFIG"
echo "    output    : $SMOKE_DIR"
echo "    extra days: $SMOKE_DAYS (on top of the config's spinup)"
echo "    sim seconds: $SMOKE_SIM_SECONDS (0 = full days mode)"
echo "    julia     : $JULIA"

# --startup-file=no: the user's ~/.julia/config/startup.jl errors under this
# project and would take the job down before any model code runs.
$JULIA --startup-file=no --project=experiments/AMIP \
    experiments/calibration/amip/smoke_test.jl
rc=$?
echo "=== $(date) smoke test exited $rc"
exit $rc
