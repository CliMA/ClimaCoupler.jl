#!/bin/bash
#SBATCH --job-name=best_cf_model_EVER
#SBATCH --gpus=1
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

# Match the alias environment
export CLIMACOMMS_DEVICE=CUDA

module load climacommon

# Ensure login-style environment
source ~/.bashrc

julia --project=experiments/ClimaEarth experiments/ClimaEarth/run_amip.jl --config_file="config/subseasonal_configs/wxquest_progedmf.yml"
