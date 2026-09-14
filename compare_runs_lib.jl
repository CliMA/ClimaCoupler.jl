# Shared loaders/diff helpers for run comparison.

"Max absolute and max relative difference between two arrays, ignoring non-finite entries."
function diffstats(a, b)
    a = Array(a);
    b = Array(b)
    size(a) == size(b) || return (nothing, nothing, 0)
    ok = isfinite.(a) .& isfinite.(b)
    any(ok) || return (0.0, 0.0, 0)
    av = Float64.(a[ok]);
    bv = Float64.(b[ok])
    absd = abs.(av .- bv)
    maxabs = maximum(absd)
    denom = max.(abs.(av), abs.(bv))
    rel = ifelse.(denom .> 0, absd ./ max.(denom, eps()), 0.0)
    return (maxabs, maximum(rel), count(absd .> 0))
end

"Collect comparable numeric arrays from a NetCDF file, keyed by variable name."
function load_nc(path)
    out = Dict{String, Any}()
    NCDatasets.NCDataset(path, "r") do ds
        for (name, v) in ds
            arr = try
                Array(v)
            catch
                continue
            end
            eltype(arr) <: Number || continue
            out[name] = arr
        end
    end
    return out
end

"Collect numeric arrays from a JLD2 file, walking nested groups."
function load_jld2(path)
    out = Dict{String, Any}()
    JLD2.jldopen(path, "r") do f
        function walk(g, prefix)
            for k in keys(g)
                key = isempty(prefix) ? k : "$prefix/$k"
                child = try
                    g[k]
                catch
                    continue
                end
                if child isa JLD2.Group
                    walk(child, key)
                elseif child isa AbstractArray && eltype(child) <: Number
                    out[key] = child
                elseif child isa Number
                    out[key] = [child]
                end
            end
        end
        walk(f, "")
    end
    return out
end

load_any(path) = endswith(path, ".nc") ? load_nc(path) : load_jld2(path)


load_any(path) = endswith(path, ".nc") ? load_nc(path) : load_jld2(path)
