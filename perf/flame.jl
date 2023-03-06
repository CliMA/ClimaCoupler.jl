# flame.jl: provides allocation breakdown for individual backtraces for single-process unthredded runs
# and check for overall allocation limits based on previous runs
# copied and modified from `ClimaAtmos/perf`

import Profile
import ProfileCanvas
using Test

cc_dir = joinpath(dirname(@__DIR__));
include(joinpath(cc_dir, "experiments", "AMIP", "modular", "cli_options.jl"));

# assuming a common driver for all tested runs
filename = joinpath(cc_dir, "experiments", "AMIP", "modular", "coupler_driver_modular.jl")

# selected runs for performance analysis and their expected allocations (based on previous runs)
run_name_list =
    ["default_modular_unthreaded", "coarse_single_modular", "target_amip_n32_shortrun", "target_amip_n1_shortrun"]
run_name = run_name_list[parse(Int, ARGS[2])]
allocs_limit = Dict()
allocs_limit["perf_default_modular_unthreaded"] = 2685744
allocs_limit["perf_coarse_single_modular"] = 3864624
allocs_limit["perf_target_amip_n32_shortrun"] = 172134848

# number of time steps used for profiling
const n_samples = 2

# flag to split coupler init from its solve
ENV["CI_PERF_SKIP_COUPLED_RUN"] = true

# pass in the correct arguments, overriding defaults with those specific to each run_name (in `pipeline.yaml`)
dict = parsed_args_per_job_id(; trigger = "--run_name $run_name")
parsed_args_prescribed = parsed_args_from_ARGS(ARGS)
parsed_args_target = dict[run_name]
global parsed_args = merge(parsed_args_target, parsed_args_prescribed) # global scope needed to recognize this definition in the coupler driver
run_name = "perf_" * run_name
parsed_args["job_id"] = run_name
parsed_args["run_name"] = run_name
parsed_args["enable_threading"] = false

@info run_name

# initialize the coupler
try
    include(filename)
catch err
    if err.error !== :exit_profile_init
        rethrow(err.error)
    end
end

#####
##### Profiling
#####

function step_coupler!(cs, n_samples)
    cs.tspan[1] = cs.model_sims.atmos_sim.integrator.t
    cs.tspan[2] = cs.tspan[1] + n_samples * cs.Δt_cpl
    solve_coupler!(cs)
end

# compile coupling loop first
step_coupler!(cs, n_samples)
Profile.clear_malloc_data()
Profile.clear()

# profile the coupling loop
prof = Profile.@profile begin
    step_coupler!(cs, n_samples)
end

(; output_dir) = simulation

# produce flamegraph
if haskey(ENV, "BUILDKITE_COMMIT") || haskey(ENV, "BUILDKITE_BRANCH")
    output_dir = "perf/output/$run_name/"
    mkpath(output_dir)
    ProfileCanvas.html_file(joinpath(output_dir, "flame.html"))
else
    ProfileCanvas.view(Profile.fetch())
end

#####
##### Allocation tests
#####

# We're grouping allocation tests here for convenience.

buffer = 1.4 # increase slightly for (nondeterministic) threaded runs

# profile the coupling loop
allocs = @allocated step_coupler!(cs, n_samples)
@timev step_coupler!(cs, n_samples)

@info "`allocs ($run_name)`: $(allocs)"

if allocs < allocs_limit[run_name] * buffer
    @info "TODO: lower `allocs_limit[$run_name]` to: $(allocs)"
end

Δallocs = allocs / allocs_limit[run_name]
@info "Allocation change (allocs/allocs_limit): $Δallocs"
percent_alloc_change = (1 - Δallocs) * 100
if percent_alloc_change ≥ 0
    @info "Allocations improved by: $percent_alloc_change %"
else
    @info "Allocations worsened by: $percent_alloc_change %"
end

@testset "Allocations limit" begin
    @test allocs ≤ allocs_limit[run_name] * buffer
end
