include("mpi/mpi_init.jl") # setup MPI context for distributed runs #hide

include(joindir(pkgdir(CA), "examples/hybrid/driver.jl"))
