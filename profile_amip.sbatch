#!/bin/bash

#SBATCH --partition=a3
#SBATCH --output="profile_amip_1proc.txt"
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-gpu=160GB
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=4
set -euo pipefail

julia --project=experiments/ClimaEarth/ -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=experiments/ClimaEarth/ -e 'using Pkg; Pkg.add("CUDA"); using CUDA; CUDA.set_runtime_version!(local_toolkit=true)'
julia --project=experiments/ClimaEarth/ -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; use_system_binary(library_names="/sw/openmpi-5.0.5/lib/libmpi", mpiexec="/sw/openmpi-5.0.5/bin/mpiexec", force=true)'

julia --project=experiments/ClimaEarth/ -e 'using Pkg; Pkg.precompile()'
julia --project=experiments/ClimaEarth/ -e 'using Pkg; Pkg.status()'

# export LD_PRELOAD=/usr/local/cuda-12.4/nsight-systems-2023.4.4/target-linux-x64/libToolsInjection64.so
unset LD_PRELOAD
export CLIMACOMMS_DEVICE="CUDA"
export CLIMACOMMS_CONTEXT="SINGLETON"

export TMPDIR=/home/ext_nefrathe_caltech_edu/nsight_tmp/${SLURM_JOB_ID}
mkdir -p $TMPDIR
rm -f /home/ext_nefrathe_caltech_edu/clima/ClimaCoupler.jl/experiments/ClimaEarth/output/amip/report.nsys-rep
nsys profile --delay 200 --trace=nvtx,mpi,cuda,osrt --output=experiments/ClimaEarth/output/amip/report julia --threads=4 --project=experiments/ClimaEarth/ experiments/ClimaEarth/run_amip.jl --config_file config/amip_configs/amip.yml --job_id amip