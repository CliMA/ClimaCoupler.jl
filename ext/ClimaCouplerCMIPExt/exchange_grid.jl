#=
# Exchange (intersection) grid

The polygons formed where the elements of the ClimaCore spectral-element (SE)
boundary space overlap the finite-volume (FV) cells of an Oceananigans grid.
Ocean and sea-ice turbulent fluxes are computed per polygon and aggregated
conservatively to both sides, so coastlines are resolved at intersection
resolution and the same areas define both flux weights and surface fractions
(issue #1838). Geometry and weights are built once on the CPU in `Float64`
(ConservativeRegridding.jl polygon clipping), cast to the simulation float
type, and moved to the device, where every per-step operation is a sparse
gather/scatter.
=#

"""
    ExchangeGrid{FT, VI, VF}

Sparse coupling between the SE boundary space, the FV (Oceananigans) surface
cells, and the intersection polygons that tile their overlap, stored in
compressed-sparse-row (CSR) form for the direction each coupling is consumed
in, so every per-step gather/scatter is a race-free segmented reduction.

SE-side weights come in two flavors over each polygon `Ω_k`:

  - gather (polygon-major): `f̄_k = Σ_n gweight[k,n] f_n` with the SEM basis
    integrals `B_{kn} = ∫_{Ω_k} ϕ_n \\, dA`
    (`ConservativeRegriddingClimaCoreExt.accumulate_principled_b`) and
    `gweight = B_{kn} / Σ_n B_{kn}` (rows sum to 1; constants preserved;
    high-order interpolation of the atmospheric state);
  - scatter (node-major): `F_n = Σ_k sweight[n,k] F_k` with
    `sweight = B̃_{kn} / (Jw)_n`, where `B̃_{kn} ≥ 0` is the area of `Ω_k`
    inside node `n`'s control volume, whose area is `(Jw)_n`; follow with
    `weighted_dss!`;
  - `node_cov[n] = Σ_{k wet} sweight[n,k]`: the fraction of node `n`'s area
    over wet ocean.

Flat index conventions: SE node `n = (e-1) Nq² + (j-1) Nq + i` (i fastest
within each element, then elements — the order of
`ConservativeRegriddingClimaCoreExt.flat_nodal_data` /
`Fields.field2array` on ClimaCore 0.15 `(v, i, j, h)` / `VIJFH` storage);
FV cell `c = (j-1) Nx + i` (`vec(OC.interior(field, :, :, Nz))` and the CR
column index).

# Fields
- `n_poly`, `n_nodes`, `n_elem`, `n_oc`: entity counts
- `gpoly_ptr`, `gnode`, `gweight`: polygon-major CSR of the node→polygon gather
- `snode_ptr`, `spoly`, `sweight`: node-major CSR of the polygon→node scatter
- `node_cov`: nodal wet coverage fractions
- `elem_of_poly`, `oc_of_poly`: SE element / FV cell owning each polygon
- `soc_ptr`, `soc_poly`: FV-cell-major CSR of the polygon→cell scatter
- `area`: geometric polygon areas [m²]
- `b_area`: SE-side scatter areas `Σ_n sweight[n,k] (Jw)_n` [m²]; equals
  `area` up to quadrature error and makes SE-side conservation statements
  exact
- `oc_wet_area`: retained (wet) polygon area per FV cell [m²]; zero for dry
  (immersed) cells and tripolar fold shadow cells
- `oc_area`: the FV grid's own cell area `Az` [m²] where `oc_wet_area > 0`;
  zero elsewhere
"""
struct ExchangeGrid{FT, VI <: AbstractVector{Int32}, VF <: AbstractVector{FT}}
    n_poly::Int
    n_nodes::Int
    n_elem::Int
    n_oc::Int
    gpoly_ptr::VI
    gnode::VI
    gweight::VF
    snode_ptr::VI
    spoly::VI
    sweight::VF
    node_cov::VF
    elem_of_poly::VI
    oc_of_poly::VI
    soc_ptr::VI
    soc_poly::VI
    area::VF
    b_area::VF
    oc_wet_area::VF
    oc_area::VF
end

Adapt.@adapt_structure ExchangeGrid

function Base.show(io::IO, eg::ExchangeGrid{FT}) where {FT}
    print(
        io,
        "ExchangeGrid{$FT}: $(eg.n_poly) polygons from $(eg.n_elem) SE elements × $(eg.n_oc) FV cells",
    )
end

"""
    on_device(arch, eg::ExchangeGrid)

Move an `ExchangeGrid` (or any Adapt-able exchange-grid object) to the memory
of the Oceananigans architecture `arch`.
"""
on_device(arch, x) = Adapt.adapt(OC.Architectures.array_type(arch), x)

# Build a CSR (ptr, perm) pair grouping `1:n_entries` by the row index
# `rows[k]` with `n_rows` rows. `perm` lists entry indices row by row;
# `ptr[r]:(ptr[r+1]-1)` are the positions of row `r` in `perm`.
function _build_csr(rows::Vector{Int}, n_rows::Int)
    perm = sortperm(rows)
    ptr = zeros(Int32, n_rows + 1)
    ptr[1] = 1
    for r in rows
        ptr[r + 1] += 1
    end
    return cumsum!(ptr, ptr), Int32.(perm)
end

"""
    build_exchange_grid(boundary_space, grid_oc)

Construct an [`ExchangeGrid`](@ref) between a ClimaCore cubed-sphere
`boundary_space` and an Oceananigans `grid_oc`, on the CPU in `Float64`:
intersect SE elements with FV cells (fold-aware on a `TripolarGrid`, so fold
shadow cells never produce polygons), drop polygons over dry (immersed)
cells, then integrate the SEM basis (`accumulate_principled_b`, gather
weights) and the control volumes (`control_volume_weights`, scatter weights
and nodal coverages) over each polygon. Move the result to the
device with [`on_device`](@ref).
"""
function build_exchange_grid(boundary_space, grid_oc)
    CRExt = get_ConservativeRegriddingCCExt()
    @assert !isnothing(CRExt) "ConservativeRegriddingClimaCoreExt must be loaded"
    GO = CRExt.GO

    boundary_space_cpu = CC.Adapt.adapt(Array, boundary_space)
    grid_oc_underlying_cpu = OC.on_architecture(OC.CPU(), underlying_grid(grid_oc))

    FT = CC.Spaces.undertype(boundary_space_cpu)
    R = Float64(CC.Spaces.topology(boundary_space_cpu).mesh.domain.radius)
    manifold = CR.Spherical(; radius = R)

    # 1. SE-element × FV-cell intersection polygons.
    dst_tree = CR.Trees.treeify(manifold, boundary_space_cpu)
    src_tree = CR.Trees.treeify(manifold, grid_oc_underlying_cpu)
    intersections = CR.intersection_areas(
        manifold,
        CR.False(),
        dst_tree,
        src_tree;
        intersection_operator = CR.IntersectionGridOperator(manifold),
    )
    elem_of_poly, oc_of_poly, polys = SparseArrays.findnz(intersections)
    n_elem, n_oc = size(intersections)
    area = [GO.area(manifold, poly) for poly in polys]

    # 2. Mark polygons over dry (immersed) surface cells, which are dropped.
    keep = trues(length(area))
    if grid_oc isa OC.ImmersedBoundaryGrid
        grid_with_mask_cpu = OC.on_architecture(OC.CPU(), grid_oc)
        Nx_oc, _, Nz_oc = size(grid_with_mask_cpu)
        for k in eachindex(keep)
            c = oc_of_poly[k]
            i, j = mod1(c, Nx_oc), (c - 1) ÷ Nx_oc + 1
            keep[k] = !OC.ImmersedBoundaries.immersed_cell(i, j, Nz_oc, grid_with_mask_cpu)
        end
    end

    # 3. Per-polygon weights: the SEM basis integrals B_kn for the gather and
    #    the control-volume weights B̃_kn for the scatter and the coverages,
    #    indexed by each kept (wet) polygon's position in the kept numbering.
    qs = CC.Spaces.quadrature_style(boundary_space_cpu)
    Nq = CC.Quadratures.degrees_of_freedom(qs)
    triangle_quad_degree = 2 * (Nq - 1)
    Jw = CRExt.se_node_weights(boundary_space_cpu) # Float64, flat nodal
    n_nodes = length(Jw)
    volumes = ControlVolumes(boundary_space_cpu)

    g_poly, g_node, g_b = Int[], Int[], Float64[] # gather COO: B_kn
    s_poly, s_node, s_w = Int[], Int[], Float64[] # scatter COO: B̃_kn / A_n
    b_sum = zeros(Float64, count(keep))
    kept_id = 0
    for k in eachindex(area)
        keep[k] || continue
        kept_id += 1
        elem = elem_of_poly[k]
        B = CRExt.accumulate_principled_b(
            manifold,
            boundary_space_cpu,
            elem,
            polys[k];
            triangle_quad_degree,
        )
        B̃ = control_volume_weights(volumes, boundary_space_cpu, elem, polys[k])
        node_offset = (elem - 1) * Nq^2
        for j in 1:Nq, i in 1:Nq
            n = node_offset + (j - 1) * Nq + i
            w = B̃[i, j] / Jw[n]
            if w != 0
                push!(s_poly, kept_id)
                push!(s_node, n)
                push!(s_w, w)
            end
            Bij = B[i, j]
            if Bij != 0
                push!(g_poly, kept_id)
                push!(g_node, n)
                push!(g_b, Bij)
                b_sum[kept_id] += Bij
            end
        end
    end
    elem_of_poly, oc_of_poly, area = elem_of_poly[keep], oc_of_poly[keep], area[keep]

    # Degenerate polygons (b_sum ≤ 0 can only come from quadrature round-off
    # on slivers) cannot be normalized for the gather; drop their COO entries
    # in both directions and zero their area so they are inert in every
    # reduction.
    degenerate = findall(<=(0), b_sum)
    if !isempty(degenerate)
        is_degenerate = falses(length(area))
        is_degenerate[degenerate] .= true
        g_keep = .!is_degenerate[g_poly]
        g_poly, g_node, g_b = g_poly[g_keep], g_node[g_keep], g_b[g_keep]
        s_keep = .!is_degenerate[s_poly]
        s_poly, s_node, s_w = s_poly[s_keep], s_node[s_keep], s_w[s_keep]
        area[degenerate] .= 0
        b_sum[degenerate] .= 0
    end

    # 4. CSR assembly for both directions.
    n_poly = length(area)
    gpoly_ptr, gperm = _build_csr(g_poly, n_poly)
    gnode = Int32.(g_node[gperm])
    gweight = [g_b[p] / b_sum[g_poly[p]] for p in gperm]

    snode_ptr, sperm = _build_csr(s_node, n_nodes)
    spoly = Int32.(s_poly[sperm])
    sweight = s_w[sperm]

    node_cov = zeros(Float64, n_nodes)
    b_area = zeros(Float64, n_poly)
    for p in eachindex(s_node)
        node_cov[s_node[p]] += s_w[p]
        b_area[s_poly[p]] += s_w[p] * Jw[s_node[p]]
    end

    soc_ptr, soc_poly = _build_csr(Vector{Int}(oc_of_poly), n_oc)
    oc_wet_area = zeros(Float64, n_oc)
    for k in eachindex(area)
        oc_wet_area[oc_of_poly[k]] += area[k]
    end
    Nx_oc, _, Nz_oc = size(grid_oc_underlying_cpu)
    oc_area = zeros(Float64, n_oc)
    for c in findall(>(0), oc_wet_area)
        i, j = mod1(c, Nx_oc), (c - 1) ÷ Nx_oc + 1
        oc_area[c] = OC.Operators.Azᶜᶜᶜ(i, j, Nz_oc, grid_oc_underlying_cpu)
    end

    return ExchangeGrid{FT, Vector{Int32}, Vector{FT}}(
        n_poly,
        n_nodes,
        n_elem,
        n_oc,
        gpoly_ptr,
        gnode,
        FT.(gweight),
        snode_ptr,
        spoly,
        FT.(sweight),
        FT.(node_cov),
        Int32.(elem_of_poly),
        Int32.(oc_of_poly),
        soc_ptr,
        soc_poly,
        FT.(area),
        FT.(b_area),
        FT.(oc_wet_area),
        FT.(oc_area),
    )
end

#=
# Gather/scatter operations

Segmented reductions over the CSR structures above: race-free and
deterministic. Each wrapper launches a single KernelAbstractions kernel on
the destination's backend, so the same code path runs on CPU and GPU.
=#

import KernelAbstractions

# Launch `kernel!` over `worksize` linear work-items on `backend`, following the
# Oceananigans convention: bake a static workgroup/worksize into the kernel
# instance (workgroup capped at 256, as in `Oceananigans.Utils.heuristic_workgroup`)
# instead of passing a dynamic `ndrange`.
launch_kernel!(kernel!, backend, worksize, args...) =
    kernel!(backend, min(worksize, 256), worksize)(args...)

# dst[r] = Σ_p w[p] src[col[p]] over CSR row r
@kernel function _csr_matvec_kernel!(dst, ptr, col, w, src)
    r = @index(Global)
    acc = zero(eltype(dst))
    @inbounds begin
        for p in ptr[r]:(ptr[r + 1] - 1)
            acc += w[p] * src[col[p]]
        end
        dst[r] = acc
    end
end

function _csr_matvec!(dst, ptr, col, w, src)
    backend = KernelAbstractions.get_backend(dst)
    launch_kernel!(_csr_matvec_kernel!, backend, length(dst), dst, ptr, col, w, src)
    return dst
end

"""
    gather_nodes_to_polys!(poly_values, eg::ExchangeGrid, nodal_values)

Average a flat SE nodal vector onto each polygon: `f̄_k = Σ_n gweight f_n`.
Rows sum to 1, so constants are preserved exactly.
"""
gather_nodes_to_polys!(poly_values, eg::ExchangeGrid, nodal_values) =
    _csr_matvec!(poly_values, eg.gpoly_ptr, eg.gnode, eg.gweight, nodal_values)

@kernel function _gather_cells_kernel!(dst, oc_of_poly, src)
    k = @index(Global)
    @inbounds dst[k] = src[oc_of_poly[k]]
end

"""
    gather_cells_to_polys!(poly_values, eg::ExchangeGrid, cell_values)

Copy the owning FV cell's value onto each polygon (direct indexing; each
polygon lies inside exactly one cell).
"""
function gather_cells_to_polys!(poly_values, eg::ExchangeGrid, cell_values)
    backend = KernelAbstractions.get_backend(poly_values)
    launch_kernel!(
        _gather_cells_kernel!,
        backend,
        eg.n_poly,
        poly_values,
        eg.oc_of_poly,
        cell_values,
    )
    return poly_values
end

"""
    scatter_polys_to_nodes!(nodal_values, eg::ExchangeGrid, poly_values)

Aggregate per-polygon values onto SE nodes with the control-volume weights:
`F_n = Σ_k sweight F_k`, the area integral of the values over node `n`'s
control volume divided by its area. The result is *coverage-weighted*: scattering a
constant yields `constant × node_cov`. Follow with `weighted_dss!` on the
receiving field. For a per-unit-wet-area result use
[`scatter_polys_to_nodes_normalized!`](@ref).
"""
scatter_polys_to_nodes!(nodal_values, eg::ExchangeGrid, poly_values) =
    _csr_matvec!(nodal_values, eg.snode_ptr, eg.spoly, eg.sweight, poly_values)

@kernel function _csr_matvec_normalized_kernel!(dst, ptr, col, w, src, cov)
    r = @index(Global)
    @inbounds begin
        c = cov[r]
        if c > 0
            acc = zero(eltype(dst))
            for p in ptr[r]:(ptr[r + 1] - 1)
                acc += w[p] * src[col[p]]
            end
            dst[r] = acc / c
        else
            dst[r] = 0
        end
    end
end

"""
    scatter_polys_to_nodes_normalized!(nodal_values, eg::ExchangeGrid, poly_values)

Like [`scatter_polys_to_nodes!`](@ref) but normalized by the nodal wet
coverage, yielding a per-unit-wet-area value: a convex combination of the
polygon values at every node with any wet coverage, and 0 at nodes with none.
"""
function scatter_polys_to_nodes_normalized!(nodal_values, eg::ExchangeGrid, poly_values)
    backend = KernelAbstractions.get_backend(nodal_values)
    launch_kernel!(
        _csr_matvec_normalized_kernel!,
        backend,
        eg.n_nodes,
        nodal_values,
        eg.snode_ptr,
        eg.spoly,
        eg.sweight,
        poly_values,
        eg.node_cov,
    )
    return nodal_values
end

@kernel function _scatter_cells_kernel!(dst, ptr, polys, area, cell_area, src)
    c = @index(Global)
    @inbounds begin
        aw = cell_area[c]
        if aw > 0
            acc = zero(eltype(dst))
            for p in ptr[c]:(ptr[c + 1] - 1)
                k = polys[p]
                acc += area[k] * src[k]
            end
            dst[c] = acc / aw
        else
            dst[c] = 0
        end
    end
end

"""
    scatter_polys_to_cells!(cell_values, eg::ExchangeGrid, poly_values,
                            poly_area = eg.b_area, cell_area = eg.oc_area)

Per-polygon values onto each FV cell: `F_c = Σ_{k ∈ c} poly_area_k F_k / cell_area_c`.
With the defaults, `F_c` is a flux per unit FV cell area whose integral `Σ_c Az_c F_c`
equals the SE-side integral `Σ_k b_area_k F_k`; pass `eg.area, eg.oc_wet_area` for the
area-weighted mean of a state. Cells with zero `cell_area` (dry, fold shadows) are set
to 0; run [`mirror_fold_partners!`](@ref) afterwards to fill the shadow copies.
"""
function scatter_polys_to_cells!(
    cell_values,
    eg::ExchangeGrid,
    poly_values,
    poly_area = eg.b_area,
    cell_area = eg.oc_area,
)
    backend = KernelAbstractions.get_backend(cell_values)
    launch_kernel!(
        _scatter_cells_kernel!,
        backend,
        eg.n_oc,
        cell_values,
        eg.soc_ptr,
        eg.soc_poly,
        poly_area,
        cell_area,
        poly_values,
    )
    return cell_values
end

#=
# Tripolar fold mirroring

On a `RightCenterFolded` `TripolarGrid` the fold-row cell `(i, Ny)` is the
same physical cell as `(Nx + 1 - i, Ny)`; the exchange grid keeps one copy as
a real cell and leaves the other a degenerate shadow with no polygons.
=#

"""
    mirror_fold_partners!(cell_values, grid)

Copy each fold-row primary cell's value into its shadow partner slot on a
`RightCenterFolded` tripolar grid; no-op for other grids. `cell_values` is a
flat vector over the `Nx × Ny` surface cells. Delegates to the internal
`mirror_fold_partners!` in `ConservativeRegriddingOceananigansExt`.
"""
mirror_fold_partners!(cell_values, grid::OC.ImmersedBoundaryGrid) =
    mirror_fold_partners!(cell_values, grid.underlying_grid)

function mirror_fold_partners!(cell_values, grid)
    CROCExt = get_ConservativeRegriddingOCExt()
    @assert !isnothing(CROCExt)
    CROCExt.mirror_fold_partners!(cell_values, grid)
    return cell_values
end

#=
# Wet-ocean surface fraction
=#

"""
    wet_ocean_fraction_field(boundary_space, eg::ExchangeGrid)

Build the wet-ocean surface fraction as a `CC.Fields.Field` on
`boundary_space` from a CPU-resident [`ExchangeGrid`](@ref): the nodal wet
coverage `node_cov` (the wet fraction of each node's control volume), made
continuous across shared element-boundary nodes with `weighted_dss!` and
clamped to [0, 1]. The complement is the land fraction; sea ice and open ocean
partition the wet fraction itself (see `FieldExchanger.align_surface_fractions!`).
"""
function wet_ocean_fraction_field(
    boundary_space,
    eg::ExchangeGrid{FT, Vector{Int32}},
) where {FT}
    CRExt = get_ConservativeRegriddingCCExt()
    field = CC.Fields.zeros(boundary_space)
    device_array_type = ClimaComms.array_type(ClimaComms.device(boundary_space))
    CRExt.vec_to_se_field!(field, Adapt.adapt(device_array_type, eg.node_cov))

    dss_buffer = Utilities.init_dss_buffer(field)
    Utilities.apply_dss!(field, dss_buffer)
    @. field = clamp(field, FT(0), FT(1))
    return field
end
