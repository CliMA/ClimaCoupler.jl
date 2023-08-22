#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --job-name=test_cheng
#SBATCH --reservation=clima
#SBATCH --mem-per-cpu=5GB
#SBATCH --ntasks=32

#SBATCH --mail-user=gdecker@caltech.edu   # email address

# Notify at the beginning, end of job and on failure.
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

set -euo pipefail # kill the job if anything fails
set -x # echo script

module purge
module load julia/1.8.1 openmpi/4.1.1 hdf5/1.12.1-ompi411 #netcdf-c/4.6.1


export BUILD_HISTORY_HANDLE=""
export CI="true"
export CLIMATEMACHINE_SETTINGS_FIX_RNG_SEED="true"
export CUDA_VERSION="11.2"
export FLAME_PLOT=""
export GKSwstype="100"
export JULIA_MAX_NUM_PRECOMPILE_FILES="100"
export JULIA_VERSION="1.8.5"
export MPI_IMPL="openmpi"
export OPENBLAS_NUM_THREADS="1"
export OPENMPI_VERSION="4.1.1"

export JULIA_MPI_BINARY=system
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export CLIMACORE_DISTRIBUTED=""
export JULIA_HDF5_PATH=""

export RUN_NAME=test_cheng


echo "--- Configure MPI"
julia -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; use_system_binary()'

echo "--- Instantiate package env"
julia --project -e 'using Pkg; Pkg.instantiate(;verbose=true)'
# julia --project -e 'using Pkg; Pkg.rm("SurfaceFluxes")'
julia --project -e 'using Pkg; Pkg.add(url="https://github.com/CliMA/SurfaceFluxes.jl", rev = "gd/test_ufs")'
julia --project -e 'using Pkg; Pkg.precompile()'
julia --project -e 'using Pkg; Pkg.status()'

# echo "--- Instantiate sea breeze env"
# julia --project=experiments/ClimaCore/sea_breeze -e 'using Pkg; Pkg.instantiate(;verbose=true)'
# julia --project=experiments/ClimaCore/sea_breeze -e 'using Pkg; Pkg.precompile()'
# julia --project=experiments/ClimaCore/sea_breeze -e 'using Pkg; Pkg.status()'

echo "--- Instantiate amip modular"
julia --project=experiments/AMIP/modular/ -e 'using Pkg; Pkg.instantiate(;verbose=true)'
# julia --project=experiments/AMIP/modular/ -e 'using Pkg; Pkg.rm("SurfaceFluxes")'
julia --project=experiments/AMIP/modular/ -e 'using Pkg; Pkg.add(url="https://github.com/CliMA/SurfaceFluxes.jl", rev = "gd/test_ufs")'
julia --project=experiments/AMIP/modular/ -e 'using Pkg; Pkg.precompile()'
julia --project=experiments/AMIP/modular/ -e 'using Pkg; Pkg.status()'
julia --project=experiments/AMIP/modular/ -e 'using Pkg; Pkg.build("MPI"); Pkg.build("HDF5")'
julia --project=experiments/AMIP/modular/ -e 'using Pkg; Pkg.build()'


# echo "--- Instantiate perf env"
# julia --project=perf/ -e 'using Pkg; Pkg.instantiate(;verbose=true)'
# julia --project=perf/ -e 'using Pkg; Pkg.precompile()'
# julia --project=perf/ -e 'using Pkg; Pkg.status()'

# echo "--- Instantiate test env"
# julia --project=test/ -e 'using Pkg; Pkg.instantiate(;verbose=true)'
# julia --project=test/ -e 'using Pkg; Pkg.precompile()'
# julia --project=test/ -e 'using Pkg; Pkg.status()'

echo "--- Download artifacts"
julia --project=artifacts -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia --project=artifacts -e 'using Pkg; Pkg.precompile()'
julia --project=artifacts -e 'using Pkg; Pkg.status()'
julia --project=artifacts artifacts/download_artifacts.jl


# run spin up
# - specify `--monthly_checkpoint true` to save monthly checkpoints of all model prognostic states

# BUSINGER (Float64)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Businger --run_name coarse_single_modular_businger_ft64 --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --precip_model 0M --job_id coarse_single_modular_businger_ft64

# BUSINGER (Float32)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Businger --run_name coarse_single_modular_businger_ft32 --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --precip_model 0M --job_id coarse_single_modular_businger_ft32

# GRYANIK (Float64)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Gryanik --run_name coarse_single_modular_gryanik_ft64 --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6  --precip_model 0M --job_id coarse_single_modular_gryanik_ft64

# GRYANIK (Float32)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Gryanik --run_name coarse_single_modular_gryanik_ft32 --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --precip_model 0M --job_id coarse_single_modular_gryanik_ft32

# Grachev (Float64)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Grachev --run_name coarse_single_modular_grachev_ft64 --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular_grachev_ft64

# Grachev (Float32)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Grachev --run_name coarse_single_modular_grachev_ft32 --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular_grachev_ft32

# Cheng (Float64)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Cheng --run_name coarse_single_modular_cheng_ft64 --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular_cheng_ft64

# Cheng (Float32)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Cheng --run_name coarse_single_modular_cheng_ft32 --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular_cheng_ft32

# Beljaars (Float64)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Beljaars --run_name coarse_single_modular_beljaars_ft64 --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular_beljaars_ft64

# Beljaars (Float34)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Beljaars --run_name coarse_single_modular_beljaars_ft32 --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular_beljaars_ft32

# Holtslag (Float64)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Holtslag --run_name coarse_single_modular_holtslag_ft64 --FLOAT_TYPE Float64 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular_holtslag_ft64

# Holtslag (Float32)
#julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --uft Holtslag --run_name coarse_single_modular_holtslag_ft32 --FLOAT_TYPE Float32 --coupled true --monthly_checkpoint true --surface_setup PrescribedSurface --moist equil --vert_diff true --rad gray --energy_check false --mode_name amip --anim true --t_end 32days --dt_save_to_sol 1days --dt_cpl 400 --dt 400secs --mono_surface false --h_elem 6 --dt_save_restart 10days --precip_model 0M --job_id coarse_single_modular_holtslag_ft32