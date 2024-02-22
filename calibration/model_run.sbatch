#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G

# Extract command-line arguments
experiment_id=$1
iteration=$2

# Find output directory
format_i=$(printf "iteration_%03d" "$iteration")
member=$(printf "member_%03d" "$SLURM_ARRAY_TASK_ID")
output=output/$experiment_id/$format_i/$member/model_log.out

# Run the forward model
srun --output=$output julia --color=no --project=experiments -e "
    using ClimaComms
    using ClimaCoupler
    ClimaComms.init(ClimaComms.context())
    import CalibrateAtmos, ClimaCoupler
    atmos_config = CalibrateAtmos.get_atmos_config($SLURM_ARRAY_TASK_ID, $iteration, \"$experiment_id\")
    CalibrateAtmos.run_forward_model(atmos_config)
"
