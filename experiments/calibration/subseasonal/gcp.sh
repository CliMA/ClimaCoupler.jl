# Script to run `experiments/calibration/subseasonal/run_calibration.jl` on GCP.

# First, setup environment
export PATH="/usr/local/bin:$PATH"
export NVHPC=/sw/nvhpc/Linux_x86_64/24.5
export HPCX_PATH=$NVHPC/comm_libs/12.4/hpcx/hpcx-2.19
# CUDA environment
export CUDA_HOME=$NVHPC/cuda/12.4
export CUDA_PATH=$CUDA_HOME
export CUDA_ROOT=$CUDA_HOME
# MPI via MPIwrapper
export MPITRAMPOLINE_LIB="/sw/mpiwrapper/lib/libmpiwrapper.so"
export OPAL_PREFIX=$HPCX_PATH/ompi
# Library paths - CUDA first, then HPC-X
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$HPCX_PATH/ompi/lib:$LD_LIBRARY_PATH
# Executable paths
export PATH=/sw/mpiwrapper/bin:$CUDA_HOME/bin:$PATH
# Julia
export PATH="/sw/julia/julia-1.11.5/bin:$PATH"

# Instantiate and generate observations
julia --project=experiments/ClimaEarth -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal/generate_observations.jl

# Run
julia --project=experiments/ClimaEarth experiments/calibration/subseasonal/run_calibration.jl
