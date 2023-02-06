#!/bin/bash
#SBATCH --ntasks=64          
#SBATCH --time=30:00:00     # walltime

set -euo pipefail # kill the job if anything fails
set -x # echo script

module purge
module load julia/1.8.1 openmpi/4.1.1 hdf5/1.12.1-ompi411 #netcdf-c/4.6.1

export CLIMACORE_DISTRIBUTED="MPI"
export JUlIA_MPI_BINARY="system"
export JULIA_HDF5_PATH=""

# run instantiate/precompile serial
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.build()'
julia --project -e 'using Pkg; Pkg.build("MPI"); Pkg.build("HDF5")'

#mpiexec julia --project coupler_driver.jl --surface_scheme monin_obukhov --moist equil --rad clearsky --microphy 0M --z_elem 50 --dz_top 3000 --dz_bottom 30 --h_elem 16 --kappa_4 1e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --alpha_rayleigh_w 10 --dt_cpl 150 --dt_rad 1hours --idealized_insolation true --FLOAT_TYPE Float64 --vert_diff true --energy_check false --mode_name amip --t_end 360days --dt_save_to_sol 10days --mono_surface false --dt_save_to_disk 1days --run_name monin_spg10_fixmsk_csky_iinsol_vereq3_dztop_Wf0p5
#mpiexec julia --project coupler_driver.jl --surface_scheme monin_obukhov --moist equil --rad clearsky --microphy 0M --z_elem 50 --dz_top 3000 --dz_bottom 30 --h_elem 16 --kappa_4 1e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --alpha_rayleigh_w 10 --dt_cpl 150 --dt_rad 1hours --idealized_insolation true --FLOAT_TYPE Float64 --vert_diff true --energy_check false --mode_name amip --t_end 360days --dt_save_to_sol 10days --mono_surface false --dt_save_to_disk 1days --run_name test_full_
#mpiexec julia --project coupler_driver.jl --surface_scheme monin_obukhov --moist equil --rad clearsky --microphy 0M --z_elem 50 --dz_top 3000 --dz_bottom 30 --h_elem 16 --kappa_4 1e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --alpha_rayleigh_w 10 --dt_cpl 150 --dt_rad 1hours --idealized_insolation true --FLOAT_TYPE Float64 --vert_diff true --energy_check false --mode_name amip --t_end 360days --dt_save_to_sol 10days --mono_surface false --dt_save_to_disk 1days --run_name test_full_maxmin
# mpiexec julia --project coupler_driver.jl --surface_scheme monin_obukhov --moist equil --rad clearsky --microphy 0M --z_elem 50 --dz_top 3000 --dz_bottom 30 --h_elem 16 --kappa_4 1e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --alpha_rayleigh_w 10 --dt_cpl 50 --dt_rad 1hours --idealized_insolation true --FLOAT_TYPE Float64 --vert_diff true --energy_check false --mode_name amip --t_end 360days --dt_save_to_sol 10days --mono_surface false --dt_save_to_disk 1days --run_name test_full_maxmin50
#mpiexec julia --project coupler_driver.jl --surface_scheme monin_obukhov --moist equil --rad clearsky --microphy 0M --z_elem 45 --dz_bottom 30 --h_elem 16 --kappa_4 1e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --dt_cpl 150 --dt_rad 6hours --idealized_insolation true --FLOAT_TYPE Float64 --vert_diff true --energy_check false --mode_name amip --t_end 360days --dt_save_to_sol 10days --mono_surface false --dt_save_to_disk 1days --run_name test_full_maxmin150_lessdz
#mpiexec julia --project coupler_driver.jl --surface_scheme monin_obukhov --moist equil --rad allskywithclear --microphy 0M --z_elem 45 --dz_bottom 30 --h_elem 16 --kappa_4 1e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --dt_cpl 150 --dt_rad 1hours --idealized_insolation false --FLOAT_TYPE Float64 --vert_diff true --energy_check false --mode_name amip --t_end 360days --dt_save_to_sol 10days --mono_surface false --dt_save_to_disk 1days --run_name test_full_maxminq150_lessdz_albedoadjusted04allsky_rinsol
mpiexec julia --project coupler_driver.jl --surface_scheme monin_obukhov --moist equil --rad allskywithclear --microphy 0M --z_elem 45 --dz_bottom 30 --h_elem 16 --kappa_4 1e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --dt_cpl 150 --dt_rad 1hours --idealized_insolation false --FLOAT_TYPE Float64 --vert_diff true --energy_check false --mode_name amip --t_end 360days --dt_save_to_sol 10days --mono_surface false --dt_save_to_disk 1days --run_name test_rinsol_ntasks64



#mpiexec julia --project coupler_driver.jl --surface_scheme monin_obukhov --moist equil --rad clearsky --microphy 0M --z_elem 50 --dz_top 3000 --dz_bottom 30 --h_elem 16 --kappa_4 1e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --alpha_rayleigh_w 10 --dt_cpl 50 --dt_rad 1hours --idealized_insolation true --FLOAT_TYPE Float64 --vert_diff true --energy_check false --mode_name amip --t_end 120days --dt_save_to_sol 10days --mono_surface false --dt_save_to_disk 1days --run_name test_full50

# goal: mpiexec julia --project=examples examples/hybrid/driver.jl --surface_scheme monin_obukhov --moist equil --rad allskywithclear --microphy 0M --z_elem 45 --dz_bottom 30 --h_elem 16 --kappa_4 1e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --dt 150secs --dt_rad 1hours --idealized_insolation false --t_end 365days --job_id longrun_aquaplanet_rhoe_equil_highres_allsky --dt_save_to_disk 10days --FLOAT_TYPE Float64