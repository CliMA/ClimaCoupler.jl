#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --job-name=mpi_amip
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=expansion

export MODULEPATH="/groups/esm/modules:$MODULEPATH"
module purge
module load climacommon/2024_04_05

export CC_PATH=$(pwd)/ # adjust this to the path of your ClimaCoupler.jl directory
export RUN_NAME=coarse_single_ft64_hourly_checkpoints_restart
export CONFIG_FILE=${CC_PATH}config/model_configs/${RUN_NAME}.yml
export RESTART_DIR=experiments/AMIP/output/amip/${RUN_NAME}_artifacts/

export OPENBLAS_NUM_THREADS=1
export JULIA_NVTX_CALLBACKS=gc
export OMPI_MCA_opal_warn_on_missing_libcuda=0
export JULIA_MAX_NUM_PRECOMPILE_FILES=100
export SLURM_KILL_BAD_EXIT=1

julia --project=experiments/AMIP/ -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=experiments/AMIP/ -e 'using Pkg; Pkg.precompile()'
julia --project=experiments/AMIP/ -e 'using Pkg; Pkg.status()'

julia --project=artifacts -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=artifacts -e 'using Pkg; Pkg.precompile()'
julia --project=artifacts -e 'using Pkg; Pkg.status()'
julia --project=artifacts artifacts/download_artifacts.jl

srun -K julia --project=experiments/AMIP/ experiments/AMIP/coupler_driver.jl --config_file $CONFIG_FILE

# restart from simulation time of 400 seconds
export RESTART_T=400

# setup the new config file with ammened checkpointing frequency
export RESTART_CONFIG_FILE=${CONFIG_FILE::-4}_tmp.yml
cp $CONFIG_FILE $RESTART_CONFIG_FILE
sed -i 's/t_end: \"800secs\"/t_end: \"3600secs\"/g' $RESTART_CONFIG_FILE

# rerun the model
srun -K julia --project=experiments/AMIP/ experiments/AMIP/coupler_driver.jl --config_file $RESTART_CONFIG_FILE --restart_dir $RESTART_DIR --restart_t $RESTART_T

# throw an error if no restart checkpoint files are found
if [ $(ls -1 $RESTART_DIR/checkpoint | wc -l) -lt 5 ]; then
  echo "Error: RESTART_DIR does not contain enough files"
  exit 1
else
  echo "Successful: RESTART_DIR contains $(ls -1 $RESTART_DIR/checkpoint | wc -l) files"
  exit 0
fi

# Trouble shooting?
# - ensure you're using the latest module file of climacommon and set MODULEPATH to the correct location
# - ensure you're using the latest version of ClimaCoupler.jl
# - did you cd to your version of ClimaCoupler.jl?
