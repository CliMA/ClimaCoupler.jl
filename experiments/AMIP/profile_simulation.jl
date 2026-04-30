# # AMIP Profiling Driver

#=
## Overview

This script runs a short AMIP simulation under the Nsight Systems profiler.
It is intended to be launched via:

```bash
nsys profile --start-later=true --capture-range=cudaProfilerApi \
    --trace=nvtx,mpi,cuda,osrt --cuda-memory-usage=true --kill=none \
    --output=path/to/report \
    julia --project=experiments/AMIP experiments/AMIP/profile_simulation.jl \
    --config_file config/benchmark_configs/amip_diagedmf_profile.yml
```

Two coupling steps are run first to trigger precompilation. The Nsight profiler
then starts capturing at the `CUDA.@profile` boundary, so that compilation and
initialization are excluded from the profile.
=#

# Load the necessary modules to run the coupled simulation
include("code_loading.jl")

import CUDA

# Get the configuration file from the command line (or manually set it here)
config_file = Input.parse_commandline(Input.argparse_settings())["config_file"]

# Set up the coupled simulation
cs = CoupledSimulation(config_file)

# Run two steps to trigger precompilation, then profile three steps
step!(cs);
step!(cs);
GC.gc()
CUDA.@profile external = true begin
    for _ in 1:3
        step!(cs)
    end
end
