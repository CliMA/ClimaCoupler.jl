import Profile
using Test
import Base: view
include("ProfileCanvasDiff.jl")
import .ProfileCanvasDiff
using JLD2

if isinteractive()
    buildkite_cc_dir = build_path = "."
else
    buildkite_number = ENV["BUILDKITE_BUILD_NUMBER"]
    buildkite_cc_dir = "/groups/esm/slurm-buildkite/climacoupler-ci/"
    build_path = "/central/scratch/esm/slurm-buildkite/climacoupler-ci/$buildkite_number/climacoupler-ci/perf/"
end

output_dir = joinpath(buildkite_cc_dir, "test/")
mkpath(output_dir)

# dummy functions to analyze
function cumsum_sqrt(x)
    y = cumsum(x)
    sqrt.(y)
end

function get_y(x)
    f = collect(1:1:10000) .* x
    cumsum(f)
end

function step_coupler!(n, y_all = [])
    for i in collect(1:1:n)
        y = get_y(i)
        push!(y_all, y)
    end
    return y_all
end

# init
step_coupler!(1)

# clear compiler allocs
Profile.clear_malloc_data()
Profile.clear()

# profile the coupling loop
prof = Profile.@profile begin
    step_coupler!(100)
end

# ref file
ref_file = joinpath(output_dir, "reference.jld2")
tracked_list = isfile(ref_file) ? load(ref_file) : Dict{String, Float64}()

# save ref file
profile_data, new_tracked_list = ProfileCanvasDiff.view(Profile.fetch(), tracked_list = tracked_list, self_count = true);

save(ref_file, new_tracked_list) # reset ref_file upon staging


"""
    find_child(flame_tree, target_name; self_count = false)

Helper function to find a particular flame tree node.
"""
function find_child(flame_tree, target_name; self_count = false)

    line = flame_tree.line
    file = flame_tree.file
    func = flame_tree.func

    node_name = "$func.$file.$line"
    node_name = self_count ? "self_count_" * node_name : node_name

    if node_name == target_name
        global out = flame_tree
    else
        if !(isempty(flame_tree.children))
            for sf in flame_tree.children
                find_child(sf, target_name, self_count = self_count)
            end
        end
    end
    return @isdefined(out) ? out : nothing
end

@testset "flame diff tests" begin
    # load the dictionary of tracked counts from the reference file
    tracked_list = isfile(ref_file) ? load(ref_file) : Dict{String, Float64}()

    test_func_name = "get_y.flame_test.jl.26"

    # test flame diff
    tracked_list["$test_func_name"] = 100  # reference value for node count with child sum
    profile_data, new_tracked_list =
        ProfileCanvasDiff.view(Profile.fetch(), tracked_list = tracked_list, self_count = false)
    node = find_child(profile_data.data["all"], test_func_name, self_count = false)
    @test node.count_change == node.count - 100

    # test flame diff with self_count
    tracked_list["self_count_$test_func_name"] = 50 # reference value for node count w/o child sum
    profile_data, new_tracked_list =
        ProfileCanvasDiff.view(Profile.fetch(), tracked_list = tracked_list, self_count = true)

    node = find_child(profile_data.data["all"], test_func_name, self_count = false)
    child_sum = sum(map(x -> x.count, node.children))
    @test node.count_change == (node.count - child_sum) - 50

    # html_file test
    ProfileCanvasDiff.html_file(
        joinpath(output_dir, "flame_diff.html"),
        build_path = build_path,
        tracked_list = tracked_list,
        self_count = false,
    )

    @test isfile(joinpath(output_dir, "flame_diff.html"))
    rm(output_dir; recursive = true, force = true)

end
