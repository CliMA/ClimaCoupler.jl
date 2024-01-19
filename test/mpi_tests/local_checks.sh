#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --job-name=mpi_restart_test
#SBATCH --reservation=clima
#SBATCH --mem=32GB
#SBATCH --ntasks=2

module purge
module load julia/1.10.0
export JULIA_MPI_BINARY=system
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export CLIMACORE_DISTRIBUTED="MPI"
export JULIA_HDF5_PATH=""

export RUN_NAME=amip_restart_mpi_test
export RESTART_DIR=experiments/AMIP/output/amip/${RUN_NAME}_artifacts/
export RESTART_T=200

julia -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; use_system_binary()'
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project -e 'using Pkg; Pkg.build("MPI")'
julia --project -e 'using Pkg; Pkg.build("HDF5")'
julia --project -e 'using Pkg; Pkg.API.precompile()'

julia --project=experiments/AMIP/ -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=experiments/AMIP/ -e 'using Pkg; Pkg.precompile()'
julia --project=experiments/AMIP/ -e 'using Pkg; Pkg.status()'

julia --project=artifacts -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=artifacts -e 'using Pkg; Pkg.precompile()'
julia --project=artifacts -e 'using Pkg; Pkg.status()'
julia --project=artifacts artifacts/download_artifacts.jl

# run spin up
# - specify `--hourly_checkpoint true` to save monthly checkpoints of all model prognostic states
mpiexec julia --color=yes --project=experiments/AMIP/ experiments/AMIP/coupler_driver.jl --run_name $RUN_NAME --coupled true   --start_date 19790101 --hourly_checkpoint true  --anim true --surface_setup PrescribedSurface --dt_cpl 200 --energy_check false --mode_name amip --mono_surface false --vert_diff true --moist equil --rad clearsky --precip_model 0M --z_elem 35 --dz_bottom 50 --h_elem 12 --kappa_4 3e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --dt 200secs --t_end 0.1days --job_id $RUN_NAME --dt_save_to_sol 1000days --dt_save_state_to_disk 10days --apply_limiter false --FLOAT_TYPE Float64

# init using a restart
# - specify the directory of the `checkpoint/` folder (i.e.,  `--restart_dir`) and time (in secs; `--restart_t`) of the restart file
mpiexec julia --color=yes --project=experiments/AMIP/ experiments/AMIP/coupler_driver.jl --run_name $RUN_NAME --coupled true --restart_dir $RESTART_DIR --restart_t $RESTART_T --start_date 19790102 --anim true --surface_setup PrescribedSurface --dt_cpl 200 --energy_check false --mode_name amip --mono_surface false --vert_diff true --moist equil --rad clearsky --precip_model 0M --z_elem 35 --dz_bottom 50 --h_elem 12 --kappa_4 3e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --dt 200secs --t_end 0.1days --job_id $RUN_NAME --dt_save_to_sol 1000days --dt_save_state_to_disk 10days --apply_limiter false --FLOAT_TYPE Float64
