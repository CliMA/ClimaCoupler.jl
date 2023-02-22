# temporarily copied and modified from https://github.com/pfitzseb/ProfileCanvas.jl

module ProfileCanvasDiff

using Profile, JSON, REPL, Pkg.Artifacts, Base64

export @profview, @profview_allocs

struct ProfileData
    data::Any
    typ::Any
end

mutable struct ProfileFrame
    func::String
    file::String # human readable file name
    path::String # absolute path
    line::Int # 1-based line number
    count::Int # number of samples in this frame
    countLabel::Union{Missing, String} # defaults to `$count samples`
    flags::UInt8 # any or all of ProfileFrameFlag
    taskId::Union{Missing, UInt}
    children::Vector{ProfileFrame}
    count_change::Float64 # fractional change in count: (new - old) / old
end

struct ProfileDisplay <: Base.Multimedia.AbstractDisplay end

function __init__()
    pushdisplay(ProfileDisplay())

    atreplinit(i -> begin
        while ProfileDisplay() in Base.Multimedia.displays
            popdisplay(ProfileDisplay())
        end
        pushdisplay(ProfileDisplay())
    end)
end

function jlprofile_data_uri(build_path)
    path = joinpath(build_path, "ProfileViewerDiff.js")
    @info "the module path is $path"
    str = read(path, String)

    return string("data:text/javascript;base64,", base64encode(str))
end

function Base.show(io::IO, ::MIME"text/html", canvas::ProfileData, build_path = "")
    id = "profiler-container-$(round(Int, rand()*100000))"

    data_uri = (jlprofile_data_uri(build_path))
    println(
        io,
        """
        <div id="$(id)" style="height: 400px; position: relative;"></div>
        <script type="module">
            const ProfileCanvas = await import('$data_uri')
            const viewer = new ProfileCanvas.ProfileViewer("#$(id)", $(JSON.json(canvas.data)), "$(canvas.typ)")
        </script>
        """,
    )
end

function Base.display(_::ProfileDisplay, canvas::ProfileData)

    file = html_file(string(tempname(), ".html"), canvas)
    url = "file://$file"

    if Sys.iswindows()
        run(`cmd /c "start $url"`)
    elseif Sys.isapple()
        run(`open $url`)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $url`)
    end
end

html_file(filename, data = Profile.fetch(); build_path = "", kwargs...) =
    html_file(filename, build_path = build_path, view(data; kwargs...)[1])

function html_file(file::AbstractString, canvas::ProfileData; build_path = "")
    @assert endswith(file, ".html")

    data_uri = jlprofile_data_uri(build_path)
    open(file, "w") do io
        id = "profiler-container-$(round(Int, rand()*100000))"

        println(
            io,
            """
            <html>
            <head>
            <style>
                #$(id) {
                    margin: 0;
                    padding: 0;
                    width: 100vw;
                    height: 100vh;
                    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif, "Apple Color Emoji", "Segoe UI Emoji";
                    overflow: hidden;
                }
                body {
                    margin: 0;
                    padding: 0;
                }
            </style>
            </head>
            <body>
                <div id="$(id)"></div>
                <script type="module">
                    const ProfileCanvas = await import('$data_uri')
                    const viewer = new ProfileCanvas.ProfileViewer("#$(id)", $(JSON.json(canvas.data)), "$(canvas.typ)")
                </script>
            </body>
            </html>
            """,
        )
    end
    return file
end

using Profile

# https://github.com/timholy/FlameGraphs.jl/blob/master/src/graph.jl
const ProfileFrameFlag = (
    RuntimeDispatch = UInt8(2^0),
    GCEvent = UInt8(2^1),
    REPL = UInt8(2^2),
    Compilation = UInt8(2^3),
    TaskEvent = UInt8(2^4),
)

function view(
    data = Profile.fetch();
    C = false,
    tracked_list = Dict{String, Int}(;),
    independent_count = false,
    kwargs...,
)

    d = Dict{String, ProfileFrame}()

    new_tracked_list = Dict{String, Int}(;)

    if VERSION >= v"1.8.0-DEV.460"
        threads = ["all", 1:Threads.nthreads()...]
    else
        threads = ["all"]
    end

    if isempty(data)
        Profile.warning_empty()
        return
    end

    lidict = Profile.getdict(unique(data))
    data_u64 = convert(Vector{UInt64}, data)
    for thread in threads
        graph = stackframetree(data_u64, lidict; thread = thread, kwargs...)
        curr_tree = make_tree(
            ProfileFrame("root", "", "", 0, graph.count, missing, 0x0, missing, ProfileFrame[], 999), #root process
            graph;
            C = C,
            tracked_list = tracked_list,
            independent_count = independent_count,
            kwargs...,
        )

        curr_tree = independent_count ? iterate_independent_change!(curr_tree, tracked_list) : curr_tree

        new_tracked_list = collect_child_counts(curr_tree)

        d[string(thread)] = curr_tree

    end


    return (ProfileData(d, "Thread"), new_tracked_list)
end

function iterate_independent_change!(flame_tree, tracked_list)

    if isempty(flame_tree.children)
        child_sum = 0
        log_independent_change!(flame_tree, child_sum, tracked_list)
    else
        child_sum = sum(map(x -> x.count, flame_tree.children))
        log_independent_change!(flame_tree, child_sum, tracked_list)
        #@show flame_tree.countLabel
        for sf in flame_tree.children
            iterate_independent_change!(sf, tracked_list)

        end
    end
    return flame_tree
end

function log_independent_change!(flame_tree, child_sum, tracked_list)
    name = string(flame_tree.file)
    func = string(flame_tree.func)
    line = string(flame_tree.line)
    file = string(flame_tree.file)
    func_sign = "independent_count_$func.$file.$line"
    current_indeppendent_count = flame_tree.count - child_sum

    old_independent_count = func_sign in keys(tracked_list) ? tracked_list[func_sign] : -999 * flame_tree.count  #TODO delete
    overall_count_change = deepcopy(flame_tree.count_change)  #TDO: do %ages
    flame_tree.count_change = (current_indeppendent_count - old_independent_count)

    ct = flame_tree.count
    flame_tree.countLabel =
        func_sign * "," *
        string(flame_tree.count) *
        " samples \n csum = $child_sum \n ct = $ct \n Δ_ref = $overall_count_change \n Δ_ref_i = $current_indeppendent_count - $old_independent_count =" *
        string(flame_tree.count_change) *
        " counts"
end

"""
    collect_child_counts(flame_tree, ct = 0, dict = Dict{String, Float64}())

Iterate over all children of a stack tree and save their names ("\$func.\$file.\$line") and
corresponding count values in a Dict.
"""
function collect_child_counts(flame_tree, ct = 0, dict = Dict{String, Float64}())
    ct += 1
    line = flame_tree.line
    file = flame_tree.file
    func = flame_tree.func
    push!(dict, "$func.$file.$line" => flame_tree.count)

    if isempty(flame_tree.children)
        # 0 upstream contribution if end child
        push!(dict, "independent_count_$func.$file.$line" => flame_tree.count)

    else
        child_sum = sum(map(x -> x.count, flame_tree.children))
        push!(dict, "independent_count_$func.$file.$line" => flame_tree.count - child_sum)
        for sf in flame_tree.children
            collect_child_counts(sf, ct, dict)
        end
    end
    return dict
end

function stackframetree(data_u64, lidict; thread = nothing, combine = true, recur = :off)
    root = combine ? Profile.StackFrameTree{StackTraces.StackFrame}() : Profile.StackFrameTree{UInt64}()
    if VERSION >= v"1.8.0-DEV.460"
        thread = thread == "all" ? (1:Threads.nthreads()) : thread
        root, _ = Profile.tree!(root, data_u64, lidict, true, recur, thread)
    else
        root = Profile.tree!(root, data_u64, lidict, true, recur)
    end
    if !isempty(root.down)
        root.count = sum(pr -> pr.second.count, root.down)
    end

    return root
end

function status(sf::StackTraces.StackFrame)
    st = UInt8(0)
    if sf.from_c && (sf.func === :jl_invoke || sf.func === :jl_apply_generic || sf.func === :ijl_apply_generic)
        st |= ProfileFrameFlag.RuntimeDispatch
    end
    if sf.from_c && startswith(String(sf.func), "jl_gc_")
        st |= ProfileFrameFlag.GCEvent
    end
    if !sf.from_c && sf.func === :eval_user_input && endswith(String(sf.file), "REPL.jl")
        st |= ProfileFrameFlag.REPL
    end
    if !sf.from_c && occursin("./compiler/", String(sf.file))
        st |= ProfileFrameFlag.Compilation
    end
    if !sf.from_c && occursin("task.jl", String(sf.file))
        st |= ProfileFrameFlag.TaskEvent
    end
    return st
end

function status(node::Profile.StackFrameTree, C::Bool)
    st = status(node.frame)
    C && return st
    # If we're suppressing C frames, check all C-frame children
    for child in values(node.down)
        child.frame.from_c || continue
        st |= status(child, C)
    end
    return st
end

function add_child(graph::ProfileFrame, node, C::Bool; tracked_list = Dict{String, Int}(;), independent_count = false)
    name = string(node.frame.file)
    func = string(node.frame.func)
    line = string(node.frame.line)
    file = basename(string(node.frame.file))
    func_sign = "$func.$file.$line"

    if func == ""
        func = "unknown"
    end

    old_count = func_sign in keys(tracked_list) ? tracked_list[func_sign] : 999
    current_count = Float64(node.count)

    frame = ProfileFrame(
        func,
        basename(name),
        name,
        node.frame.line,
        node.count,
        missing,
        status(node, C),
        missing,
        ProfileFrame[],
        current_count - Float64(old_count),
    )

    push!(graph.children, frame)

    return frame
end

function make_tree(
    graph,
    node::Profile.StackFrameTree;
    C = false,
    tracked_list = Dict{String, Int}(;),
    independent_count = false,
)
    for child_node in sort!(collect(values(node.down)); rev = true, by = node -> node.count)
        # child not a hidden frame
        if C || !child_node.frame.from_c
            child = add_child(graph, child_node, C, tracked_list = tracked_list, independent_count = independent_count)
            make_tree(child, child_node; C = C, tracked_list = tracked_list, independent_count = independent_count)
        else
            make_tree(graph, child_node, tracked_list = tracked_list, independent_count = independent_count)
        end
    end

    return graph
end

"""
    @profview f(args...) [C = false]

Clear the Profile buffer, profile `f(args...)`, and view the result graphically.

The default of `C = false` will only show Julia frames in the profile graph.
"""
macro profview(ex, args...)
    return quote
        Profile.clear()
        Profile.@profile $(esc(ex))
        view(; $(esc.(args)...))
    end
end

## Allocs

"""
    @profview_allocs f(args...) [sample_rate=0.0001] [C=false]

Clear the Profile buffer, profile `f(args...)`, and view the result graphically.
"""
macro profview_allocs(ex, args...)
    sample_rate_expr = :(sample_rate = 0.0001)
    for arg in args
        if Meta.isexpr(arg, :(=)) && length(arg.args) > 0 && arg.args[1] === :sample_rate
            sample_rate_expr = arg
        end
    end
    if isdefined(Profile, :Allocs)
        return quote
            Profile.Allocs.clear()
            Profile.Allocs.@profile $(esc(sample_rate_expr)) $(esc(ex))
            view_allocs()
        end
    else
        return :(@error "This version of Julia does not support the allocation profiler.")
    end
end

function view_allocs(_results = Profile.Allocs.fetch(); C = false)
    results = _results::Profile.Allocs.AllocResults
    allocs = results.allocs

    allocs_root = ProfileFrame("root", "", "", 0, 0, missing, 0x0, missing, ProfileFrame[], 0)
    counts_root = ProfileFrame("root", "", "", 0, 0, missing, 0x0, missing, ProfileFrame[], 0)
    for alloc in allocs
        this_allocs = allocs_root
        this_counts = counts_root

        for sf in Iterators.reverse(alloc.stacktrace)
            if !C && sf.from_c
                continue
            end
            file = string(sf.file)
            this_counts′ = ProfileFrame(
                string(sf.func),
                basename(file),
                file,
                sf.line,
                0,
                missing,
                0x0,
                missing,
                ProfileFrame[],
                0,
            )
            ind = findfirst(
                c -> (c.func == this_counts′.func && c.path == this_counts′.path && c.line == this_counts′.line),
                this_allocs.children,
            )

            this_counts, this_allocs = if ind === nothing
                push!(this_counts.children, this_counts′)
                this_allocs′ = deepcopy(this_counts′)
                push!(this_allocs.children, this_allocs′)

                (this_counts′, this_allocs′)
            else
                (this_counts.children[ind], this_allocs.children[ind])
            end
            this_allocs.count += alloc.size
            this_allocs.countLabel = memory_size(this_allocs.count)
            this_counts.count += 1
            this_allocs.count_change = 0.6
        end

        alloc_type = replace(string(alloc.type), "Profile.Allocs." => "")
        ind = findfirst(c -> (c.func == alloc_type), this_allocs.children)
        if ind === nothing
            push!(
                this_allocs.children,
                ProfileFrame(
                    alloc_type,
                    "",
                    "",
                    0,
                    this_allocs.count,
                    memory_size(this_allocs.count),
                    ProfileFrameFlag.GCEvent,
                    missing,
                    ProfileFrame[],
                    0.6,
                ),
            )
            push!(
                this_counts.children,
                ProfileFrame(alloc_type, "", "", 0, 1, missing, ProfileFrameFlag.GCEvent, missing, ProfileFrame[], 1),
            )
        else
            this_counts.children[ind].count += 1
            this_allocs.children[ind].count += alloc.size
            this_allocs.children[ind].countLabel = memory_size(this_allocs.children[ind].count)
            this_allocs.children[ind].count_change = 0.6
        end

        counts_root.count += 1
        allocs_root.count += alloc.size
        allocs_root.countLabel = memory_size(allocs_root.count)
        allocs_root.count_change = 0.6
    end

    d = Dict{String, ProfileFrame}("size" => allocs_root, "count" => counts_root)

    return ProfileData(d, "Allocation")
end

const prefixes = ["bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"]
function memory_size(size)
    i = 1
    while size > 1000 && i + 1 < length(prefixes)
        size /= 1000
        i += 1
    end
    return string(round(Int, size), " ", prefixes[i])
end


end
