# flame_diff.jl: provides allocation breakdown for individual backtraces for single-process unthredded runs
# and check for fractional change in allocation compared to the last staged run
if isinteractive()
    buildkite_branch = buildkite_commit = buildkite_number = "none"
    buildkite_cc_dir = build_path = "."
    perf_run_no = "1"
else
    buildkite_branch = ENV["BUILDKITE_BRANCH"]
    buildkite_commit = ENV["BUILDKITE_COMMIT"]
    buildkite_number = ENV["BUILDKITE_BUILD_NUMBER"]
    buildkite_cc_dir = "/groups/esm/slurm-buildkite/climacoupler-ci/"
    build_path = "/central/scratch/esm/slurm-buildkite/climacoupler-ci/$buildkite_number/climacoupler-ci/perf/"
    perf_run_no = ARGS[2]
end

cwd = pwd()
@info "build_path is: $build_path"

import Profile
using Test
import Base: view
include("ProfileCanvasDiff.jl")
import .ProfileCanvasDiff
using JLD2

cc_dir = joinpath(dirname(@__DIR__));
include(joinpath(cc_dir, "experiments", "AMIP", "moist_mpi_earth", "cli_options.jl"));

# assuming a common driver for all tested runs
filename = joinpath(cc_dir, "experiments", "AMIP", "moist_mpi_earth", "coupler_driver_modular.jl")

# selected runs for performance analysis and their expected allocations (based on previous runs)
run_name_list = ["default_modular", "coarse_single_modular", "target_amip_n32_shortrun", "target_amip_n1_shortrun"]
run_name = run_name_list[parse(Int, perf_run_no)]

# number of time steps used for profiling
n_samples = 2

# flag to split coupler init from its solve
ENV["CI_PERF_SKIP_COUPLED_RUN"] = true

# pass in the correct arguments, overriding defaults with those specific to each run_name (in `pipeline.yaml`)
dict = parsed_args_per_job_id(; trigger = "--run_name $run_name")
parsed_args_prescribed = parsed_args_from_ARGS(ARGS)
parsed_args_target = dict[run_name]
parsed_args = merge(parsed_args_target, parsed_args_prescribed) # global scope needed to recognize this definition in the coupler driver
run_name = "perf_diff_" * run_name
parsed_args["job_id"] = run_name
parsed_args["run_name"] = run_name
parsed_args["enable_threading"] = false

@info run_name

if !isinteractive()
    function step_coupler!(cs, n_samples)
        cs.tspan[1] = cs.model_sims.atmos_sim.integrator.t
        cs.tspan[2] = cs.tspan[1] + n_samples * cs.Δt_cpl
        solve_coupler!(cs)
    end

    try # initialize the coupler
        ENV["CI_PERF_SKIP_COUPLED_RUN"] = true
        include(filename)
    catch err
        if err.error !== :exit_profile_init
            rethrow(err.error)
        end
    end
end
#####
##### Profiling
#####

# obtain the stacktree from the last saved file in `buildkite_cc_dir`
ref_file = joinpath(buildkite_cc_dir, "$run_name.jld2")
tracked_list = isfile(ref_file) ? load(ref_file) : Dict{String, Float64}()

# compile coupling loop first
step_coupler!(cs, n_samples)

if !isinteractive()
    # clear compiler allocs
    Profile.clear_malloc_data()
    Profile.clear()
end

# profile the coupling loop
prof = Profile.@profile begin
    step_coupler!(cs, n_samples + 1000)
end

# produce flamegraph with colors highlighting the allocation differences relative to the last saved run
# profile_data
if haskey(ENV, "BUILDKITE_COMMIT") || haskey(ENV, "BUILDKITE_BRANCH")
    output_dir = "perf/output/$run_name"
    mkpath(output_dir)
    ProfileCanvasDiff.html_file(joinpath(output_dir, "flame_diff.html"), build_path = build_path)
    ProfileCanvasDiff.html_file(
        joinpath(output_dir, "flame_diff_self_count.html"),
        build_path = build_path,
        self_count = true,
    )
end

# save (and reset) the stack tree if this is running on the `staging` branch
@info "This branch is: $buildkite_branch, commit $buildkite_commit"
profile_data, new_tracked_list = ProfileCanvasDiff.view(Profile.fetch(), tracked_list = tracked_list, self_count = true);
#if buildkite_branch == "staging"
isfile(ref_file) ?
mv(ref_file, joinpath(buildkite_cc_dir, "flame_reference_file.$run_name.$buildkite_commit.jld2"), force = true) :
nothing
save(ref_file, new_tracked_list) # reset ref_file upon staging
#end



# tesr


# rr.data["all"].children[1].children[1].children[1].children[1].children[1].children[1].children[1].children[1].children[1].children[1].children[1].children[1].children[2].children

# function step_coupler!(cs, n)
#     m = add_n!(n)
#     times_minus_n!(m)

# end

# function add_n!(n)
#     n + n
# end

# function times_minus_n!(n)
#     f = 300 * zeros(3000)
#     b = (n * n - n) .+ f
#     sqrt.(b)
# end
# cs = []
