#!/bin/bash
#SBATCH --gpus=1
#SBATCH --job-name=ml_flux
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=train_%j.out

module purge
module load climacommon/2025_12_22

unset LD_LIBRARY_PATH

julia --project=. -e '
    using Pkg
    Pkg.instantiate(; verbose=true)
    Pkg.add(name="WeightInitializers", version="1.0.5")
    Pkg.status()
'

# Train (writes latest_run.txt with the run directory path)
USE_GPU=true julia -t 8 --project=. train_flux_correction.jl

echo ""
echo "============================================================"
echo "Training complete. Running evaluation..."
echo "============================================================"

# Evaluate (reads latest_run.txt automatically, saves plots into the run dir)
USE_GPU=true julia -t 32 --project=. evaluate_flux_model.jl
