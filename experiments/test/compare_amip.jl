# compare.jl provides function to recursively compare complex objects while also
# allowing for some numerical tolerance.

import ClimaComms
import ClimaAtmos as CA
import ClimaCore as CC
import ClimaUtilities.TimeVaryingInputs: AbstractTimeVaryingInput
import NCDatasets

"""
    _error(arr1::AbstractArray, arr2::AbstractArray; ABS_TOL = 100eps(eltype(arr1)))

We compute the error in this way:
- when the absolute value is larger than ABS_TOL, we use the absolute error
- in the other cases, we compare the relative errors

Gather to CPU before comparing so GPU restart tests do not hit scalar-indexing
errors from `maximum` / reductions on `CuArray`.
"""
function _error(arr1::AbstractArray, arr2::AbstractArray; ABS_TOL = 100eps(eltype(arr1)))
    # There are some parameters, e.g. Obukhov length, for which Inf
    # is a reasonable value (implying a stability parameter in the neutral boundary layer
    # regime, for instance). We account for such instances with the `isfinite` function.
    arr1 = Array(arr1) .* isfinite.(Array(arr1))
    arr2 = Array(arr2) .* isfinite.(Array(arr2))
    diff = abs.(arr1 .- arr2)
    denominator = abs.(arr1)
    return ifelse.(denominator .> ABS_TOL, diff ./ denominator, diff)
end

"""
    compare(v1, v2; name = "", ignore = Set([:rc]))

Return whether `v1` and `v2` are the same (up to floating point errors).
`compare` walks through all the properties in `v1` and `v2` until it finds
that there are no more properties. At that point, `compare` tries to match the
resulting objects. When such objects are arrays with floating point, `compare`
defines a notion of `error` that is the following: when the absolute value is
less than `100eps(eltype)`, `error = absolute_error`, otherwise it is relative
error. The `error` is then compared against a tolerance.
Keyword arguments
=================
- `name` is used to collect the name of the property while we go recursively
  over all the properties. You can pass a base name.
- `ignore` is a collection of `Symbol`s that identify properties that are
  ignored when walking through the tree. This is useful for properties that
  are known to be different (e.g., `output_dir`).
`:rc` is some CUDA/CuArray internal object that we don't care about
"""
function compare(
    v1::T1,
    v2::T2;
    name = "",
    ignore = Set([:rc]),
) where {
    T1 <: Union{CC.Fields.FieldVector, CC.Spaces.AbstractSpace, NamedTuple, CA.AtmosCache},
    T2 <: Union{CC.Fields.FieldVector, CC.Spaces.AbstractSpace, NamedTuple, CA.AtmosCache},
}
    pass = true
    return _compare(pass, v1, v2; name, ignore)
end

function _compare(pass, v1::T, v2::T; name, ignore) where {T}
    properties = filter(x -> !(x in ignore), propertynames(v1))
    # Dictionaries.jl Dictionaries store entries in :slots/:keys/:vals arrays with
    # #undef. Walking those throws UndefRefError.
    if (:keys in properties || :slots in properties) &&
       (:vals in properties || :values in properties)
        return pass && print_maybe(v1 == v2, "$name differs")
    end
    if isempty(properties)
        pass &= _compare(v1, v2; name, ignore)
    else
        # Recursive case. Catch exceptions so failures report the field path
        # instead of a multi-megabyte AtmosCache type stacktrace.
        for p in properties
            child_name = "$(name).$(p)"
            try
                pass &= _compare(
                    pass,
                    getproperty(v1, p),
                    getproperty(v2, p);
                    name = child_name,
                    ignore,
                )
            catch e
                # Avoid printing full AtmosCache type parameters (can be megabytes).
                println("$child_name threw $(nameof(typeof(e)))")
                pass = false
            end
        end
    end
    return pass
end

function _compare(v1::T, v2::T; name, ignore) where {T}
    return print_maybe(v1 == v2, "$name differs")
end

function _compare(v1::T, v2::T; name, ignore) where {T <: Union{AbstractString, Symbol}}
    # What we can safely print without filling STDOUT
    return print_maybe(v1 == v2, "$name differs: $v1 vs $v2")
end

function _compare(v1::T, v2::T; name, ignore) where {T <: Number}
    # We check with triple equal so that we also catch NaNs being equal
    return print_maybe(v1 === v2, "$name differs: $v1 vs $v2")
end

# We ignore NCDatasets. They contain a lot of state-ful information
function _compare(pass, v1::T, v2::T; name, ignore) where {T <: NCDatasets.NCDataset}
    return pass
end

# TimeVaryingInputs are reinitialized on restart (not checkpointed); they embed
# NCDatasets / DataHandlers whose pointers and caches differ across runs.
function _compare(
    pass,
    ::AbstractTimeVaryingInput,
    ::AbstractTimeVaryingInput;
    name,
    ignore,
)
    return pass
end

# FieldVector is AbstractArray, but `Array(fv)` iterates with scalar getindex — forbidden
# on GPU. After restart, concrete FieldVector types often differ (space type params), so
# the same-type property walker does not apply and the AbstractArray path would fire.
# Compare named Field components instead.
function _compare_fieldvector(pass, v1, v2; name, ignore)
    keys1, keys2 = propertynames(v1), propertynames(v2)
    if keys1 != keys2
        return pass && print_maybe(false, "$name differs: keys $keys1 vs $keys2")
    end
    for k in keys1
        k in ignore && continue
        pass &= _compare(
            pass,
            getproperty(v1, k),
            getproperty(v2, k);
            name = "$(name).$(k)",
            ignore,
        )
    end
    return pass
end

function _compare(pass, v1::CC.Fields.FieldVector, v2::CC.Fields.FieldVector; name, ignore)
    return _compare_fieldvector(pass, v1, v2; name, ignore)
end

# Same-type FieldVector must beat the generic `where {T}` walker.
function _compare(pass, v1::T, v2::T; name, ignore) where {T <: CC.Fields.FieldVector}
    return _compare_fieldvector(pass, v1, v2; name, ignore)
end

# ClimaCore Fields: compare parent arrays (gathered in AbstractArray leaf), never walk
# Field/grid properties. Same-type method beats generic property walker.
function _compare_climacore_field(pass, v1, v2; name, ignore)
    return pass && _compare(parent(v1), parent(v2); name, ignore)
end

function _compare(pass, v1::T, v2::T; name, ignore) where {T <: CC.Fields.Field}
    return _compare_climacore_field(pass, v1, v2; name, ignore)
end

function _compare(pass, v1::CC.Fields.Field, v2::CC.Fields.Field; name, ignore)
    return _compare_climacore_field(pass, v1, v2; name, ignore)
end

function _compare(
    v1::T,
    v2::T;
    name,
    ignore,
) where {T <: CC.Fields.Field{<:CC.DataLayouts.AbstractData{<:Real}}}
    return _compare(parent(v1), parent(v2); name, ignore)
end

function _compare(pass, v1::T, v2::T; name, ignore) where {T <: CC.DataLayouts.AbstractData}
    return pass && _compare(parent(v1), parent(v2); name, ignore)
end

# Handle views
function _compare(
    pass,
    v1::SubArray{FT},
    v2::SubArray{FT};
    name,
    ignore,
) where {FT <: AbstractFloat}
    return pass && _compare(collect(v1), collect(v2); name, ignore)
end

# Route any AbstractArray through the leaf compare (avoiding the generic property
# walk). This is essential for GPU arrays: without this, `_compare(pass, v1::T, v2::T)`
# recurses into CuArray internals (`.ptr`, `.rc`, …) and either hits scalar-indexing
# errors or compares allocation identity rather than values.
# StaticArrays (e.g. SVector) also go through here; `_error` gathers via `Array(...)`.
# FieldVector is handled above — do not `Array(fv)` on GPU.
function _compare(pass, v1::AbstractArray, v2::AbstractArray; name, ignore)
    return pass && _compare(v1, v2; name, ignore)
end

# Floating-point arrays: tolerance-based comparison on CPU after gather.
function _compare(
    v1::AbstractArray{FT},
    v2::AbstractArray{FT};
    name,
    ignore,
) where {FT <: AbstractFloat}
    a1, a2 = Array(v1), Array(v2)
    if size(a1) != size(a2)
        return print_maybe(false, "$name differs: size $(size(a1)) vs $(size(a2))")
    end
    isempty(a1) && return true
    error = maximum(_error(a1, a2))
    return print_maybe(error <= 100eps(eltype(v1)), "$name error: $error")
end

# Non-float arrays (integer lookup tables, Bool masks, etc.): exact equality on CPU.
function _compare(v1::AbstractArray, v2::AbstractArray; name, ignore)
    a1, a2 = Array(v1), Array(v2)
    if size(a1) != size(a2)
        return print_maybe(false, "$name differs: size $(size(a1)) vs $(size(a2))")
    end
    return print_maybe(a1 == a2, "$name differs")
end

# Device pointers are allocation identity — they differ after every restart by design.
# Array *values* are checked by the AbstractArray dispatches above; this just
# prevents the generic property walker from treating a new-allocation pointer as a fail.
function _compare(v1::Ptr, v2::Ptr; name, ignore)
    return true
end

function _compare(pass, v1::T1, v2::T2; name, ignore) where {T1, T2}
    println("$name differs: type mismatch $(nameof(T1)) vs $(nameof(T2))")
    return false
end

function print_maybe(exp, what)
    exp || println(what)
    return exp
end
