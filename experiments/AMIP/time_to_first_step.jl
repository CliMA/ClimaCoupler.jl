# # AMIP Time-to-First-Step Driver

#=
## Overview

This script measures the wall time to load packages, set up a coupled
simulation, and take its first coupling step. Julia startup is not included.
The timings are written to `time_to_first_step.md` and
`time_to_first_step.csv` in the artifacts directory.

```bash
julia --project=experiments/AMIP experiments/AMIP/time_to_first_step.jl \
    --config_file config/amip_configs/amip.yml
```
=#

load_time = @elapsed include("code_loading.jl")

import ClimaCoupler: CoupledSimulation, Input, step!
import ClimaCoupler.Utilities: @timed_log, format_duration

# Get the configuration file from the command line (or manually set it here)
config_file = Input.parse_commandline(Input.argparse_settings())["config_file"]

# `@eval` keeps compilation of the timed call inside the timed region
setup_time = @elapsed @eval cs = CoupledSimulation(config_file)
step_time = @elapsed @eval @timed_log "First coupling step" step!(cs)

timings = [
    "Load packages" => load_time,
    "Set up CoupledSimulation" => setup_time,
    "First coupling step" => step_time,
    "Total" => load_time + setup_time + step_time,
]

open(joinpath(cs.dir_paths.artifacts_dir, "time_to_first_step.md"), "w") do io
    println(io, "### Time to first step: `$config_file`\n")
    upstream_packages = get(ENV, "UPSTREAM_PACKAGES", "")
    isempty(upstream_packages) || println(io, "Upstream packages: `$upstream_packages`\n")
    println(io, "| Phase | Wall time |")
    println(io, "|---|---|")
    for (phase, seconds) in timings
        println(io, "| $phase | $(format_duration(seconds)) |")
    end
end

open(joinpath(cs.dir_paths.artifacts_dir, "time_to_first_step.csv"), "w") do io
    println(io, "phase,seconds")
    for (phase, seconds) in timings
        println(io, "$phase,$seconds")
    end
end
