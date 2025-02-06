#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --job-name=mpi_amip
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=16G
#SBATCH --partition=expansion

export MODULEPATH="/groups/esm/modules:$MODULEPATH"
module purge
module load climacommon/2024_10_09

export CC_PATH=$(pwd)/ # adjust this to the path of your ClimaCoupler.jl directory
export JOB_ID=amip_coarse_ft64_hourly_checkpoints_restart
export CONFIG_FILE=${CC_PATH}config/ci_configs/${JOB_ID}.yml
export RESTART_DIR=experiments/ClimaEarth/output/${JOB_ID}/checkpoints/

export OPENBLAS_NUM_THREADS=1
export JULIA_MAX_NUM_PRECOMPILE_FILES=100
export SLURM_KILL_BAD_EXIT=1

julia --project=experiments/ClimaEarth/ -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=experiments/ClimaEarth/ -e 'using Pkg; Pkg.precompile()'
julia --project=experiments/ClimaEarth/ -e 'using Pkg; Pkg.status()'

srun -K julia --project=experiments/ClimaEarth/ experiments/ClimaEarth/run_amip.jl --config_file $CONFIG_FILE --job_id $JOB_ID

# restart from simulation time of 400 seconds
export RESTART_T=400

# setup the new config file with ammened checkpointing frequency
export RESTART_CONFIG_FILE=${CONFIG_FILE::-4}_tmp.yml
cp $CONFIG_FILE $RESTART_CONFIG_FILE
sed -i 's/t_end: \"800secs\"/t_end: \"3600secs\"/g' $RESTART_CONFIG_FILE

# rerun the model
srun -K julia --project=experiments/ClimaEarth/ experiments/ClimaEarth/run_amip.jl --config_file $RESTART_CONFIG_FILE --job_id $JOB_ID --restart_dir $RESTART_DIR --restart_t $RESTART_T

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
