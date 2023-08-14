#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --job-name=coupler_driver_modular_test
#SBATCH --reservation=clima
#SBATCH --mem=20GB
#SBATCH --ntasks=8

#SBATCH -J "coupler_driver_modular_test"   # job name
#SBATCH --mail-user=gdecker@caltech.edu   # email address

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

set -euo pipefail # kill the job if anything fails
set -x # echo script

module purge
module load julia/1.8.5 openmpi/4.1.1 hdf5/1.12.1-ompi411

export JULIA_MPI_BINARY=system
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export CLIMACORE_DISTRIBUTED="MPI"
export JULIA_HDF5_PATH=""

export RUN_NAME=AMIP_modular_Float64+monthly_checkpoint
export RESTART_DIR=experiments/AMIP/modular/output/amip/${RUN_NAME}_artifacts/
export RESTART_T=200

julia -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; use_system_binary()'
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project -e 'using Pkg; Pkg.build("MPI")'
julia --project -e 'using Pkg; Pkg.build("HDF5")'
julia --project -e 'using Pkg; Pkg.API.precompile()'

julia --project=experiments/AMIP/modular/ -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=experiments/AMIP/modular/ -e 'using Pkg; Pkg.precompile()'
julia --project=experiments/AMIP/modular/ -e 'using Pkg; Pkg.status()'

julia --project=artifacts -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=artifacts -e 'using Pkg; Pkg.precompile()'
julia --project=artifacts -e 'using Pkg; Pkg.status()'
julia --project=artifacts artifacts/download_artifacts.jl

# run spin up
# - specify `--monthly_checkpoint true` to save monthly checkpoints of all model prognostic states

# BUSINGER (Float64)
mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Businger --run_name coarse_single_modular --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# BUSINGER (Float32)
mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Businger --run_name coarse_single_modular --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# GRYANIK (Float64)
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Gryanik --run_name coarse_single_modular --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# GRYANIK (Float32)
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Gryanik --run_name coarse_single_modular --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# Grachev (Float64)
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Grachev --run_name coarse_single_modular --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# Grachev (Float32)
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Grachev --run_name coarse_single_modular --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# Cheng (Float64)
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Cheng --run_name coarse_single_modular --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# Cheng (Float32)
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Cheng --run_name coarse_single_modular --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# Beljaars (Float64)
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Beljaars --run_name coarse_single_modular --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# Beljaars (Float64)
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Beljaars --run_name coarse_single_modular --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# Holtslag (Float64)
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Holtslag --run_name coarse_single_modular --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular

# Holtslag (Float32)
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Holtslag --run_name coarse_single_modular --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular