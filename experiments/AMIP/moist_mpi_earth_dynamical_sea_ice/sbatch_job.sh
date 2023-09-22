#!/bin/bash


#SBATCH --ntasks=2
#SBATCH --time=10:00:00     # walltime

set -euo pipefail # kill the job if anything fails
set -x # echo script

module purge
module load julia/1.9.3 openmpi/4.1.1 hdf5/1.12.1-ompi411 #netcdf-c/4.6.1

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export JULIA_MPI_BINARY=system
export JULIA_CUDA_USE_BINARYBUILDER=false

# run instantiate/precompile serial
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.build()'

export CLIMACORE_DISTRIBUTED="MPI"
mpiexec julia --project coupler_driver.jl
