# # Work-queue multi-run driver

#=
Run coupled simulations from a shared list of configs, several back to back in a
single Julia session, so that package loading and GPU kernel compilation are paid
only once per job. Intended for generating many short forecasts (e.g. one per
ERA5 initialization date) that share a model configuration and differ only in
`start_date` and `coupler_output_dir`.

Any number of independent jobs (workers) can run this driver on the same config
list and claim directory. Each worker repeatedly claims the next config that is
neither complete nor claimed, runs it, and exits when nothing is left, so workers
can start whenever the scheduler frees nodes.

Usage (from the repository root):
```bash
julia --project=experiments/AMIP experiments/AMIP/run_simulations_sequential.jl \
    --config_list configs.txt --claim_dir claims [--walltime_seconds 21600]
```
`configs.txt` lists one coupler YAML config per line (blank lines and lines
starting with `#` are ignored).

- A config is claimed by atomically creating `<claim_dir>/<config name>/`, which
  records the owning `PBS_JOBID` in `owner.txt`.
- A completed run writes `artifacts/run_complete.txt` in its output directory and
  is never rerun.
- A failed run writes `failed.txt` into its claim and is not retried by this
  batch. A worker stops after `--max_consecutive_failures` failures in a row.
- With `--walltime_seconds`, a worker only claims a new config if its longest run
  so far still fits before the deadline, measured from `JOB_START_EPOCH` if set.
  The first run (which includes compilation) is excluded; until a second run
  finishes, `--est_run_seconds` is used instead.
- Claims left by a job that died mid-run are not reclaimed here; the submit
  script releases claims whose owning job is no longer queued or running.

Errors must be raised on all MPI ranks for a worker to continue; an error raised
on a single rank leaves the others waiting in a collective until the job
walltime expires.
=#

include("code_loading.jl")

import ArgParse
import CUDA
import Dates

const COMPLETE_MARKER = "run_complete.txt"

function parse_worker_args()
    settings = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table! settings begin
        "--config_list"
        help = "Text file with one coupler YAML config path per line"
        arg_type = String
        required = true
        "--claim_dir"
        help = "Directory shared by all workers for claiming configs"
        arg_type = String
        required = true
        "--walltime_seconds"
        help = "Job walltime in seconds; 0 disables the time budget"
        arg_type = Int
        default = 0
        "--est_run_seconds"
        help = "Assumed duration of one run before any run has finished"
        arg_type = Int
        default = 2400
        "--max_consecutive_failures"
        help = "Stop this worker after this many failed runs in a row"
        arg_type = Int
        default = 2
    end
    args = ArgParse.parse_args(ARGS, settings)
    # `Input.get_coupler_config_dict` re-parses `ARGS` with the coupler's own
    # settings, which reject the options above.
    empty!(ARGS)
    lines = strip.(readlines(args["config_list"]))
    args["configs"] = filter(l -> !isempty(l) && !startswith(l, "#"), lines)
    return args
end

function marker_path(config_file)
    config_dict = Input.get_coupler_config_dict(config_file)
    output_dir_root = joinpath(config_dict["coupler_output_dir"], config_dict["job_id"])
    return joinpath(output_dir_root, "artifacts", COMPLETE_MARKER)
end

claim_path(claim_dir, config_file) = joinpath(claim_dir, splitext(basename(config_file))[1])

function try_claim(path)
    try
        mkdir(path)
    catch err
        ispath(path) && return false
        rethrow()
    end
    owner = "$(get(ENV, "PBS_JOBID", "unknown")) $(gethostname()) $(Dates.now())\n"
    write(joinpath(path, "owner.txt"), owner)
    return true
end

"""Index of the next config this worker claimed, or 0 if none is left."""
function claim_next(configs, markers, claim_dir)
    for (i, config_file) in enumerate(configs)
        isfile(markers[i]) && continue
        try_claim(claim_path(claim_dir, config_file)) && return i
    end
    return 0
end

function free_device_memory!()
    GC.gc(true)
    CUDA.functional() && CUDA.reclaim()
    return nothing
end

function run_worker(args)
    configs = args["configs"]
    claim_dir = args["claim_dir"]
    walltime = args["walltime_seconds"]
    job_start = parse(Float64, get(ENV, "JOB_START_EPOCH", string(time())))
    deadline = job_start + walltime

    ctx = ClimaComms.context()
    ClimaComms.init(ctx)
    is_root = ClimaComms.iamroot(ctx)
    is_root && mkpath(claim_dir)
    markers = marker_path.(configs)

    @info "Worker $(get(ENV, "PBS_JOBID", "local")) starting on $(length(configs)) configs"
    # The first run of a worker includes compilation, so it is left out of the estimate.
    run_times = Float64[]
    n_done, n_failed, consecutive_failures = 0, 0, 0
    while true
        next = 0
        if is_root
            longest_run =
                length(run_times) >= 2 ? maximum(run_times[2:end]) : args["est_run_seconds"]
            if walltime > 0 && time() + 1.15 * longest_run > deadline
                @info "Not enough walltime left for another run (longest run $(round(longest_run)) s)"
            else
                next = claim_next(configs, markers, claim_dir)
            end
        end
        next = ClimaComms.bcast(ctx, next)
        next == 0 && break

        config_file = configs[next]
        claim = claim_path(claim_dir, config_file)
        @info "Starting: $config_file"
        t0 = time()
        try
            cs = CoupledSimulation(config_file)
            run!(cs)
            elapsed = time() - t0
            if is_root
                mkpath(dirname(markers[next]))
                write(markers[next], "completed $(Dates.now()) in $(round(elapsed; digits = 1)) s\n")
            end
            push!(run_times, elapsed)
            n_done += 1
            consecutive_failures = 0
            @info "Finished in $(round(elapsed; digits = 1)) s: $config_file"
        catch err
            is_root && write(joinpath(claim, "failed.txt"), sprint(showerror, err))
            n_failed += 1
            consecutive_failures += 1
            @error "Failed: $config_file" exception = (err, catch_backtrace())
        end
        free_device_memory!()
        if consecutive_failures >= args["max_consecutive_failures"]
            @error "Stopping worker after $consecutive_failures consecutive failures"
            break
        end
    end
    @info "Worker finished: $n_done completed, $n_failed failed"
    return nothing
end

run_worker(parse_worker_args())
