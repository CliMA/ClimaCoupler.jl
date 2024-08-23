import MPI
#=
# if running locally:
module purge
module load julia/1.10.1 cuda/12.2 ucx/1.14.1_cuda-12.2 nsight-systems/2023.4.1
export OPENBLAS_NUM_THREADS=1
export JULIA_NVTX_CALLBACKS=gc
export OMPI_MCA_opal_warn_on_missing_libcuda=0
export JULIA_MAX_NUM_PRECOMPILE_FILES=100
export JULIA_LOAD_PATH=<directory with your JuliaProject.toml>

export CLIMACORE_DISTRIBUTED="MPI"

julia -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; use_system_binary()'

julia -e 'using Pkg; Pkg.instantiate(;verbose=true)'
julia -e 'using Pkg; Pkg.precompile()'
julia -e 'using Pkg; Pkg.status()'

# then precompile all packages and run with `julia` and the default env
=#

function runmpi(file; ntasks = 1)
    MPI.mpiexec() do cmd
        Base.run(
            `$cmd -n $ntasks $(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $file`;
            wait = true,
        )
        true
    end
end

if !Sys.iswindows()
    @info "tests started"
    runmpi(joinpath(@__DIR__, "regridder_mpi_tests.jl"), ntasks = 2)
    runmpi(joinpath(@__DIR__, "checkpointer_mpi_tests.jl"), ntasks = 1)
end
