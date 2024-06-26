# flame.jl: provides allocation breakdown for individual backtraces for single-process unthredded runs
# and check for overall allocation limits based on previous runs
# copied and modified from `ClimaAtmos/perf`

import Test: @test, @testset
import ClimaAtmos as CA
import Profile
import ProfileCanvas
import YAML

cc_dir = joinpath(dirname(@__DIR__));
config_dir = joinpath(cc_dir, "config", "perf_configs");
include(joinpath(cc_dir, "experiments", "ClimaEarth", "cli_options.jl"));

# assuming a common driver for all tested runs
filename = joinpath(cc_dir, "experiments", "ClimaEarth", "run_amip.jl")

# currently tested jobs and their allowed allocation limits
allocs_limit = Dict()
allocs_limit["perf_default_unthreaded"] = 8638304
allocs_limit["perf_coarse_single_ft64"] = 18280800
allocs_limit["perf_target_amip_n32_shortrun"] = 172134848

# number of time steps used for profiling
const n_samples = 2

# import parsed command line arguments
global parsed_args = parse_commandline(argparse_settings())

# select the configuration file and extract the job ID
config_file =
    parsed_args["config_file"] =
        isinteractive() ? "../config/perf_configs/perf_default_unthreaded.yml" : parsed_args["config_file"]
job_id = parsed_args["job_id"]

# import config setup
config_dict = YAML.load_file(config_file)

# global scope needed to recognize this definition in the coupler driver
parsed_args = merge(config_dict, parsed_args)

# disable threading
parsed_args["enable_threading"] = false

# flag to split coupler init from its solve
ENV["CI_PERF_SKIP_COUPLED_RUN"] = true

@info job_id

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

# produce flamegraph
if haskey(ENV, "BUILDKITE_COMMIT") || haskey(ENV, "BUILDKITE_BRANCH")
    output_dir = "perf/output/$job_id/"
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

@info "`allocs ($job_id)`: $(allocs)"

if allocs < allocs_limit[job_id] * buffer
    @info "TODO: lower `allocs_limit[$job_id]` to: $(allocs)"
end

Δallocs = allocs / allocs_limit[job_id]
@info "Allocation change (allocs/allocs_limit): $Δallocs"
percent_alloc_change = (1 - Δallocs) * 100
if percent_alloc_change ≥ 0
    @info "Allocations improved by: $percent_alloc_change %"
else
    @info "Allocations worsened by: $percent_alloc_change %"
end

@testset "Allocations limit" begin
    @test allocs ≤ allocs_limit[job_id] * buffer
end
