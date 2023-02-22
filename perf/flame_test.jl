import Profile
using Test
import Base: view
include("ProfileCanvasDiff.jl")
import .ProfileCanvasDiff
using JLD2

output_dir = "test/"
mkpath(output_dir)

# dummy functions to analyse

function cumsum_sqrt(x)
    y = cumsum(x)
    sqrt.(y)
end

function get_y(x)
    f = collect(1:1:3500) .* x
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
step_coupler!(100)

# clear compiler allocs
# Profile.clear_malloc_data()
# Profile.clear()

# profile the coupling loop
prof = Profile.@profile begin
    step_coupler!(50)
end

# html_file test
ProfileCanvasDiff.html_file(joinpath(output_dir, "flame_diff.html"), build_path = ".");

# ref file
ref_file = joinpath(output_dir, "reference.jld2")
tracked_list = isfile(ref_file) ? load(ref_file) : Dict{String, Float64}()

# save ref file
profile_data, new_tracked_list =
    ProfileCanvasDiff.view(Profile.fetch(), tracked_list = tracked_list, independent_count = true);

save(ref_file, new_tracked_list) # reset ref_file upon staging

# # load new ref file
# tracked_list = isfile(ref_file) ? load(ref_file) : Dict{String, Float64}()

# # run new allocations
# prof = Profile.@profile begin
#     step_coupler!(300)
# end
# profile_data_2, new_tracked_list_2 =
#     ProfileCanvasDiff.view(Profile.fetch(), tracked_list = tracked_list, independent_count = true);
