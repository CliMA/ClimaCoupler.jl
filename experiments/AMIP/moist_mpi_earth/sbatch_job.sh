#!/bin/bash
#SBATCH --ntasks=2          
#SBATCH --time=10:00:00     # walltime

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

mpiexec julia --project coupler_driver.jl --energy_check=false --mode_name=amip --anim=true --t_end=32days --dt_save_to_sol=1days --dt_cpl=200 --mono_surface=false --h_elem=6 --dt_save_restart=5days

# goal: mpiexec julia --project=examples examples/hybrid/driver.jl --surface_scheme monin_obukhov --moist equil --rad allskywithclear --microphy 0M --z_elem 45 --dz_bottom 30 --h_elem 16 --kappa_4 1e16 --rayleigh_sponge true --alpha_rayleigh_uh 0 --dt 150secs --dt_rad 1hours --idealized_insolation false --t_end 365days --job_id longrun_aquaplanet_rhoe_equil_highres_allsky --dt_save_to_disk 10days --FLOAT_TYPE Float64