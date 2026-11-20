"""
    PhaseProbe

Phase-resolved timing of a Julia program's startup. Each phase records wall time, exact
Julia compile time (`Base.cumulative_compile_time_ns`) and, when sampled, a split of wall
time into GPU compile, inference, host codegen, GC, I/O and execution.

The split uses the sampling profiler with C frames, which sees both host codegen and GPU
kernel compilation; SnoopCompile sees inference only.
"""
module PhaseProbe

import Profile
import Base.StackTraces

const CATEGORIES = (:gpu_compile, :inference, :codegen, :gc, :io, :execute)

mutable struct Phase
    name::String
    wall_ns::UInt64
    compile_ns::UInt64
    recompile_ns::UInt64
    samples::Dict{Symbol, Int}
    nsamples::Int
end

const PHASES = Phase[]

"When non-empty, a flat profile for each phase is written here."
const DUMPDIR = Ref("")

function dump_phase(name)
    isempty(DUMPDIR[]) && return
    mkpath(DUMPDIR[])
    fn = joinpath(DUMPDIR[], replace(name, r"[^A-Za-z0-9]+" => "_") * ".flat.txt")
    open(fn, "w") do io
        Profile.print(io; C = true, format = :flat, sortedby = :count, mincount = 5)
    end
end

"Category => frame-name and file-path patterns, in precedence order; the first match wins."
const RULES = [
    :gpu_compile => (
        paths = [
            "/GPUCompiler/",
            "/CUDA/src/compiler/",
            "/LLVM/src/",
            "/KernelAbstractions/src/compiler",
        ],
        funcs = [
            "cufunction",
            "ptxas",
            "nvlink",
            "compile_method_instance",
            "emit_llvm",
            "emit_asm",
        ],
    ),
    :inference => (
        paths = [
            "compiler/typeinfer.jl",
            "compiler/abstractinterpretation.jl",
            "compiler/inferencestate.jl",
            "compiler/ssair/",
        ],
        funcs = ["jl_type_infer", "typeinf", "abstract_call", "abstract_eval"],
    ),
    :codegen => (
        paths = String[],
        funcs = [
            "jl_create_native",
            "jl_emit_code",
            "jl_emit_codeinst",
            "_jl_compile_codeinst",
            "jl_compile_method_internal",
            "jl_generate_fptr",
            "emit_function",
            "jl_compile_codeinst",
            "llvm::",
            "ir_to_codeinf",
            "jl_compile_extern_c",
            "jl_jit_",
        ],
    ),
    :gc => (paths = String[], funcs = ["jl_gc_", "gc_mark", "jl_gc_collect", "sweep_"]),
    :io => (
        paths = ["/NCDatasets/", "/HDF5/", "/NetCDF_jll/", "/Downloads/", "/JLD2/"],
        funcs = ["nc_", "H5", "uv_fs_", "jl_fs_"],
    ),
]

_matches(s, pats) = any(p -> occursin(p, s), pats)

"""
    classify(ips, cache)

Return the highest-precedence category in `RULES` matched by any frame of the stack `ips`,
or `:execute`. `cache` memoizes the category of each instruction pointer.
"""
function classify(ips::Vector{UInt64}, cache::Dict{UInt64, Symbol})
    best = length(RULES) + 1          # index into RULES; larger == lower precedence
    for ip in ips
        c = get(cache, ip, nothing)
        if c === nothing
            c = :execute
            frames = StackTraces.lookup(ip)
            for fr in frames
                fname = string(fr.func)
                fpath = string(fr.file)
                for (i, (cat, pats)) in enumerate(RULES)
                    if _matches(fpath, pats.paths) || _matches(fname, pats.funcs)
                        c = cat
                        break
                    end
                end
                c === :execute || break
            end
            cache[ip] = c
        end
        if c !== :execute
            i = findfirst(r -> r.first === c, RULES)
            i < best && (best = i)
        end
    end
    best > length(RULES) ? :execute : RULES[best].first
end

"""
    samples(data)

Split a profile buffer fetched with `include_meta = false` into one instruction-pointer
vector per sample.
"""
function samples(data::Vector{UInt})
    out = Vector{UInt64}[]
    cur = UInt64[]
    for d in data
        if d == 0
            isempty(cur) || (push!(out, copy(cur)); empty!(cur))
        else
            push!(cur, d)
        end
    end
    isempty(cur) || push!(out, cur)
    return out
end

"""
    @phase name expr

Run `expr` and record its wall time, Julia compile time and sampled category breakdown
under `name`. Return the value of `expr`.
"""
macro phase(name, expr)
    :(@phase_impl($(esc(name)), true, $(esc(expr))))
end

"""
    @phase_nosample name expr

Like `@phase`, without the sampling profiler. Use it for package loading.
"""
macro phase_nosample(name, expr)
    :(@phase_impl($(esc(name)), false, $(esc(expr))))
end

macro phase_impl(name, sample, expr)
    quote
        local nm = $(esc(name))
        local smp = $(esc(sample))
        println(stdout, "=== phase: ", nm, " ===")
        flush(stdout)
        GC.gc()
        smp && Profile.clear()
        local (c0, r0) = Base.cumulative_compile_time_ns()
        local t0 = time_ns()
        local val = smp ? (Profile.@profile $(esc(expr))) : $(esc(expr))
        local t1 = time_ns()
        local (c1, r1) = Base.cumulative_compile_time_ns()
        local counts = Dict{Symbol, Int}(c => 0 for c in CATEGORIES)
        local n = 0
        if smp
            local cache = Dict{UInt64, Symbol}()
            for ips in samples(Profile.fetch(include_meta = false))
                counts[classify(ips, cache)] += 1
                n += 1
            end
        end
        push!(PHASES, Phase(nm, t1 - t0, c1 - c0, r1 - r0, counts, n))
        if smp
            dump_phase(nm)
            Profile.clear()
        end
        println(
            stdout,
            "    ",
            nm,
            ": ",
            round((t1 - t0) / 1e9, digits = 1),
            " s wall, ",
            round((c1 - c0) / 1e9, digits = 1),
            " s julia-compile",
        )
        flush(stdout)
        val
    end
end

_s(ns) = round(ns / 1e9, digits = 1)

_row(io, a, b, c, d, e) = println(
    io,
    rpad(string(a), 30),
    lpad(string(b), 11),
    lpad(string(c), 7),
    lpad(string(d), 14),
    lpad(string(e), 13),
)

function report(io::IO = stdout)
    total = sum(p -> p.wall_ns, PHASES; init = UInt64(0))
    total == 0 && return
    println(io, "\n", "="^76)
    println(io, "STARTUP BREAKDOWN   (total ", _s(total), " s)")
    println(io, "="^76)
    _row(io, "phase", "wall(s)", "%", "jlcompile(s)", "recomp(s)")
    println(io, "-"^76)
    for p in PHASES
        _row(
            io,
            p.name,
            _s(p.wall_ns),
            round(100 * p.wall_ns / total, digits = 1),
            _s(p.compile_ns),
            _s(p.recompile_ns),
        )
    end

    println(io, "\n", "="^92)
    println(io, "SAMPLED ATTRIBUTION  (seconds of each phase's wall time)")
    println(io, "="^92)
    print(io, rpad("phase", 30))
    for c in CATEGORIES
        print(io, lpad(string(c), 12))
    end
    println(io, lpad("nsamp", 8))
    println(io, "-"^92)
    for p in PHASES
        print(io, rpad(p.name, 30))
        for c in CATEGORIES
            frac = p.nsamples == 0 ? 0.0 : p.samples[c] / p.nsamples
            print(io, lpad(string(round(frac * _s(p.wall_ns), digits = 1)), 12))
        end
        println(io, lpad(string(p.nsamples), 8))
    end
    println(io)
end

end # module
