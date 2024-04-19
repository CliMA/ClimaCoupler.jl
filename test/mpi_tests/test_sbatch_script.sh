#!/bin/bash

# This script submits a job to the Slurm scheduler and waits for it to finish. It
# reports the job status every 30 seconds until the job completes. If the job
# fails or is terminated, the script prints an error message and exits with a
# non-zero status code. This is used by Buildkite to determine whether the job
# truly succeeded or failed.

# Submit the sbatch script and capture its job ID
JOB_ID=$(sbatch  test/mpi_tests/local_checks.sh | awk '{print $4}')
echo "Submitted job with ID: $JOB_ID, output log: slurm-$JOB_ID.out"
START_TIME=$(date +%s)
# Loop until the job finishes
while true; do
    # Check the status of the job
    STATUS=$(scontrol show job $JOB_ID | grep -oP 'JobState=\K\S+')
    sleep 30
    ELAPSED_TIME=$(( $(date +%s) - $START_TIME ))
    # If the job status is 'PD' (pending) or 'R' (running), wait and continue checking
    if [ "$STATUS" == "" ] || [ "$STATUS" == "PENDING" ] || [ "$STATUS" == "RUNNING" ]; then
        echo "Job is still running... Elapsed time: $ELAPSED_TIME seconds."
    # If the job status is 'CF' (completed successfully), print success message and exit
    elif [ "$STATUS" == "COMPLETED" ]; then
        echo "Job completed successfully."
        exit 0
    # If the job status is anything else, print error message and exit
    else
        echo "Error: Job failed or terminated. See slurm-$JOB_ID.out for more information."
        cat "slurm-$JOB_ID.out"
        exit 1
    fi
done
