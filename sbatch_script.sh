#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=mpi_restart_test
#SBATCH --mem=32GB
#SBATCH --ntasks=32
module purge
module load julia/1.8.5 openmpi/4.1.1 hdf5/1.12.1-ompi411
export JULIA_MPI_BINARY=system
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:=1}
export CLIMACORE_DISTRIBUTED="MPI"
export JULIA_HDF5_PATH=""

export RESTART_T=18316800
export RUN_NAME=n32_amip_restart_from_${RESTART_T}
export RESTART_DIR=/central/scratch/esm/slurm-buildkite/climacoupler-reports/31/climacoupler-reports/experiments/AMIP/modular/output/amip/amip_longrun_target_topo_artifacts/
export NEW_RESTART_DIR=experiments/AMIP/modular/output/amip/${RUN_NAME}_artifacts/edmf_cache_hack/


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
# mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --run_name $RUN_NAME --coupled true   --start_date 19790101 --monthly_checkpoint true  --anim true --surface_setup PrescribedSurface --dt_cpl 200 --energy_check false --mode_name amip --mono_surface false --vert_diff true --moist equil --rad clearsky --precip_model 0M --z_elem 35 --dz_bottom 50 --h_elem 12 --kappa_4 3e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --dt 200secs --t_end 1days --job_id $RUN_NAME --dt_save_to_sol 1000days --dt_save_to_disk 10days --apply_limiter false --FLOAT_TYPE Float64
mkdir $NEW_RESTART_DIR
cp -rf $RESTART_DIR/checkpoint/ $NEW_RESTART_DIR/

# run edmf hack to add ρatke to overwrite the atmos restart at time RESTART_T from RESTART_DIR
# append variable with one process
julia --project edmf_append_variable_hack.jl

wait

# init using a restart
# - specify the directoryof the `checkpoint/` folder (i.e.,  `--restart_dir`) and time (in secs; `--restart_t`) of the restart file
mpiexec julia --color=yes --project=experiments/AMIP/modular/ experiments/AMIP/modular/coupler_driver_modular.jl --run_name amip_longrun_target_offline_edmf_topo --coupled true --run_name $RUN_NAME --coupled true --restart_dir $NEW_RESTART_DIR --restart_t $RESTART_T --start_date 19790801 --monthly_checkpoint true --surface_setup PrescribedSurface --dt_cpl 100 --energy_check false --mode_name amip --mono_surface false --vert_diff true --moist equil --rad clearsky --precip_model 0M --turbconv diagnostic_edmfx --edmfx_entr_detr true --topography Earth --topo_smoothing true --use_reference_state false --z_elem 35 --dz_bottom 50 --h_elem 12 --kappa_4 4e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --dt 100secs --t_end 400days --job_id amip_longrun_target_offline_edmf_topo --dt_save_to_sol 1000days --dt_save_to_disk 1days --apply_limiter false --FLOAT_TYPE Float64  --post_process false
