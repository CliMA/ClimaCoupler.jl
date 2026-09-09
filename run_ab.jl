# Run several coupled configurations back to back in ONE process.
#
# Compilation is per-process, not per-simulation, and it dominates these short
# runs (~40 min of a ~48 min run). Building each configuration in the same
# process pays it once, and guarantees every arm runs on identical compiled
# code -- which matters for an A/B whose whole point is comparing output.
#
#   julia --project=experiments/CMIP -t 4 run_ab.jl ov_off ov_on
#
# Each argument names a config in config/ci_configs/<name>.yml. Outputs land in
# output/<name>/ as usual, one per job_id.

# Load the full CMIP stack, exactly as experiments/CMIP/run_simulation.jl does.
# `using ClimaCoupler` alone is not enough: the component models live behind
# package extensions, so without them the atmosphere type will not resolve.
include(joinpath(@__DIR__, "experiments", "CMIP", "code_loading.jl"))

import CUDA
import Dates

const Interfacer = ClimaCoupler.Interfacer
const Input = ClimaCoupler.Input

# Capture the run list before anything touches ARGS: `get_coupler_config_dict`
# re-parses ARGS through argparse, so it must hold `--config_file <path>` (what
# run_simulation.jl passes) and not our bare names.
const RUN_NAMES = copy(ARGS)
config_path(name) = joinpath("config", "ci_configs", "$(name).yml")

function set_args_for!(name)
    empty!(ARGS)
    append!(ARGS, ["--config_file", config_path(name)])
    return nothing
end

stamp() = Dates.format(Dates.now(), "HH:MM:SS")

# --- pre-flight: fail on a bad config in seconds, not after a 40 min compile ---
for name in RUN_NAMES
    path = config_path(name)
    isfile(path) || error("no such config: $path")
    set_args_for!(name)
    cfg = Input.get_coupler_config_dict(path)
    @info "preflight $name" job_id = cfg["job_id"] step_concurrently =
        cfg["step_concurrently"] overlap_slow_surfaces = cfg["overlap_slow_surfaces"] t_end =
        cfg["t_end"] dt_cpl = cfg["dt_cpl"] dt_ocean = cfg["dt_ocean"]
end
@info "preflight OK for $(join(RUN_NAMES, ", "))"

"""
Build and run one configuration, then drop it and reclaim the device memory.

The sim is built and released inside this function so the only reference dies
before the next build; two full coupled states will not fit on one GPU.
"""
function run_one(name)
    set_args_for!(name)
    @info "=== $name: building at $(stamp()) ==="
    t_build = @elapsed cs = Interfacer.CoupledSimulation(config_path(name))
    @info "=== $name: built in $(round(t_build, digits = 1)) s; running ==="
    t_run = @elapsed ClimaCoupler.SimCoordinator.run!(cs)
    @info "=== $name: done at $(stamp()), build $(round(t_build, digits = 1)) s, run $(round(t_run, digits = 1)) s ==="
    return (; name, t_build, t_run)
end

results = NamedTuple[]
for name in RUN_NAMES
    # A run may legitimately blow up (NaNs, instability). Record it and carry on
    # to the remaining configurations rather than losing the whole batch.
    try
        push!(results, run_one(name))
    catch e
        @error "=== $name FAILED at $(stamp()) ===" exception = (e, catch_backtrace())
        push!(results, (; name, t_build = NaN, t_run = NaN))
    end
    GC.gc()
    GC.gc()
    CUDA.reclaim()
    @info "GPU free after reclaim: $(round(CUDA.available_memory() / 2^30, digits = 2)) GiB"
end

println("\n", "="^62)
for r in results
    status = isnan(r.t_run) ? "  <-- FAILED" : ""
    println(
        rpad(r.name, 14),
        " build ",
        lpad(round(r.t_build, digits = 1), 8),
        " s   run ",
        lpad(round(r.t_run, digits = 1), 8),
        " s",
        status,
    )
end
println("="^62)
