#!/bin/bash
#SBATCH --gpus=1
#SBATCH --job-name=amip
#SBATCH --time=24:00:00

module purge
module load climacommon

unset LD_LIBRARY_PATH

# CONFIG_FILE='/net/sampo/data1/cchristo/clima/copies/ClimaCoupler.jl/config/amip_configs/amip_land.yml'
CONFIG_FILE='/home/chew/git_repos/ClimaCoupler.jl/config/nightly_configs/amip_coarse.yml'

# rm -rf experiments/ClimaEarth/output/wxquest_diagedmf/checkpoints

# Assuming the script is run from the top level of ClimaCoupler.jl
CLIMAEARTH_PATH='experiments/ClimaEarth/'
DRIVER=${CLIMAEARTH_PATH}'run_amip.jl'


export CLIMACOMMS_DEVICE=CUDA
# export CLIMACOMMS_CONTEXT="SINGLETON"
export CLIMACOMMS_CONTEXT="MPI"

# Use the relative CLIMAEARTH_PATH for the project
julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.instantiate(;verbose=true)'

# Keep the absolute path for the local ClimaAtmos, as specified
julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.develop(Pkg.PackageSpec(;name="ClimaAtmos", path="/home/chew/git_repos/ClimaAtmos.jl"))'
# julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.develop(Pkg.PackageSpec(;name="ClimaParams", path="/net/sampo/data1/cchristo/clima/ClimaParams.jl"))'
# julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="ClimaAtmos", rev="rc/ogwd_gpu"))'
julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="ClimaParams", rev="rc/gw"))'

# julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.update("ClimaAtmos")'
# julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="ClimaLand", rev="main"))'
# julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="ClimaLand", rev="36a35577de3c3d1f67133c70c6c4269781f046b7"))'

# julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="Thermodynamics", rev="main"))'

# julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="ClimaCore", rev="main"))'
# julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="SurfaceFluxes", rev="main"))'
# julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="RRTMGP", rev="main"))'
# julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add(Pkg.PackageSpec(;name="ClimaTimeSteppers", rev="main"))'



julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.add("MPI"); Pkg.add("CUDA")'
julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.resolve()'
julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.precompile()'
# Add CUDA runtime precompilation for GPU jobs
julia --project=$CLIMAEARTH_PATH -e 'using CUDA; CUDA.precompile_runtime()'
julia --project=$CLIMAEARTH_PATH -e 'using Pkg; Pkg.status()'


srun julia --project=$CLIMAEARTH_PATH $DRIVER --config_file $CONFIG_FILE
# srun -n 2 --gpus-per-task=1 julia --project=$CLIMAEARTH_PATH $DRIVER --config_file $CONFIG_FILE
