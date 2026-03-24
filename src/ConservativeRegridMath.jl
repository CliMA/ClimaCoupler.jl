"""
    ConservativeRegridMath

Standalone linear-algebra routines for **conservative (finite-volume) regridding** once the
intersection-area matrix `A` and cell areas are known.

This matches the convention used by `ConservativeRegridding.jl`: `A[i, j]` is the area of
intersection between destination cell `i` and source cell `j`, with `size(A) == (n_dest, n_src)`.

See also: [`conservative_regrid_forward!`](@ref), [`conservative_regrid_transpose!`](@ref).
"""
module ConservativeRegridMath

using LinearAlgebra: LinearAlgebra, mul!

"""
    conservative_regrid_forward!(dest, A, src, dest_areas)

In-place conservative remap **source → destination**:

```math
d_i = \\frac{\\sum_j A_{ij}\\, s_j}{a^d_i}
```

- `A`: ``n_{dest} \\times n_{src}`` (e.g. sparse `SparseMatrixCSC` or dense)
- `dest_areas`: length ``n_{dest}`` (physical area of each destination cell; same units as entries of `A`)

Mutates `dest` (length ``n_{dest}``). `src` has length ``n_{src}``.

This is algebraically the same as `ConservativeRegridding.regrid!(dest, regridder, src)` when
`regridder.intersections == A` and `regridder.dst_areas == dest_areas`.

All arguments should live on the same device when using GPU-backed arrays (e.g. CUDA); `./=` follows
standard Julia array semantics.

!!! note "CUDA sparse matrices"
    For `CUDA.CUSPARSE.CuSparseMatrixCSC`, `LinearAlgebra.mul!` may scalar-index and error; the CMIP
    coupler path uses `CUDA.CUSPARSE.mv!` instead. If you call these functions directly with GPU sparse
    `A`, use CUSPARSE or copy `A` to CPU.
"""
function conservative_regrid_forward!(dest, A, src, dest_areas)
    mul!(dest, A, src)
    dest ./= dest_areas
    return dest
end

"""
    conservative_regrid_transpose!(src_buf, A, dest, src_areas)

In-place conservative remap **destination → source** using the transpose of the same intersection matrix:

```math
s_j = \\frac{\\sum_i A_{ij}\\, d_i}{a^s_j}
```

i.e. `src_buf .= (transpose(A) * dest) ./ src_areas`.

- `src_buf`: length ``n_{src}`` (output)
- `dest`: length ``n_{dest}``
- `src_areas`: length ``n_{src}`` (physical area of each source cell)

Equivalent to `ConservativeRegridding.regrid!(src_buf, transpose(regridder), dest)` for a
`Regridder` built from `A` with matching `src_areas` / `dest_areas`.
"""
function conservative_regrid_transpose!(src_buf, A, dest, src_areas)
    mul!(src_buf, transpose(A), dest)
    src_buf ./= src_areas
    return src_buf
end

"""
    precompute_conservative_weights(A, dst_areas, src_areas)

Bundle precomputed conservative regridding data for documentation and optional dispatch.

Returns a `NamedTuple` `(; A, dst_areas, src_areas)` with **copies** of the arrays so later
mutations (e.g. normalization) on a `ConservativeRegridding.Regridder` do not alias this bundle.

`A` should be `n_dest × n_src` intersection areas; `dst_areas` and `src_areas` match
`ConservativeRegridding.Regridder` fields of the same names.
"""
function precompute_conservative_weights(A, dst_areas, src_areas)
    return (;
        A = copy(A),
        dst_areas = copy(dst_areas),
        src_areas = copy(src_areas),
    )
end

"""
    precompute_from_regridder(remapper; copy_arrays::Bool = true)

Convenience wrapper for initialization: read `intersections`, `dst_areas`, and `src_areas`
from any regridder object with those field names (e.g. `ConservativeRegridding.Regridder`).

- `copy_arrays = true` (default): return [`precompute_conservative_weights`](@ref) (CPU/GPU-safe
  `copy` on the same device as the regridder data).
- `copy_arrays = false`: return a named tuple **aliasing** the regridder fields (no extra memory;
  do not mutate the regridder independently). Prefer this on GPU after
  `Oceananigans.on_architecture` so `A` and area vectors stay on-device with no duplicate sparse matrix.

Use once at setup; timestepping then uses [`conservative_regrid_forward!`](@ref) /
[`conservative_regrid_transpose!`](@ref) with arrays on a **single** device consistent with `A`.
"""
function precompute_from_regridder(remapper; copy_arrays::Bool = true)
    if copy_arrays
        return precompute_conservative_weights(
            remapper.intersections,
            remapper.dst_areas,
            remapper.src_areas,
        )
    else
        return (;
            A = remapper.intersections,
            dst_areas = remapper.dst_areas,
            src_areas = remapper.src_areas,
        )
    end
end

end # module
