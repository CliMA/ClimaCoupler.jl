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
output="output/$experiment_id/$format_i/$member/model_log.out"

experiment_dir = joinpath("experiments", "amip_coupled")
COUPLER_OUTPUT_DIR = joinpath("experiments","$format_i","$member")

# Run the forward model
srun --output=$output julia --color=no --project=experiments -e "
    using ClimaComms
    using ClimaCoupler
    ClimaComms.init(ClimaComms.context())
    include(\"coupler_interface.jl\")
    import CalibrateAtmos
    coupled_simulation = get_coupler_sim($SLURM_ARRAY_TASK_ID, $iteration, \"$experiment_id\");
    run_forward_model(coupled_simulation);
"
