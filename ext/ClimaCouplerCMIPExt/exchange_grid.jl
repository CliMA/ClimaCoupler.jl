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

Each process holds the polygons of the SE elements it stores, so the SE-side
gathers and scatters touch only its own nodes; the FV side exchanges cells
with the other processes through a [`CellExchange`](@ref).
=#

"""
    ExchangeGrid{FT, VI, VF, CX}

Sparse coupling between the SE elements of this process, the FV (Oceananigans)
surface cells, and the intersection polygons that tile their overlap, stored in
compressed-sparse-row (CSR) form for the direction each coupling is consumed
in, so every per-step gather/scatter is a race-free segmented reduction.

SE-side weights come in three flavors over each polygon `Ω_k`:

  - gather (polygon-major): `f̄_k = Σ_n gweight[k,n] f_n` with the SEM basis
    integrals `B_{kn} = ∫_{Ω_k} ϕ_n \\, dA`
    (`ConservativeRegriddingClimaCoreExt.accumulate_principled_b`) and
    `gweight = B_{kn} / Σ_n B_{kn}` (rows sum to 1; constants preserved;
    high-order interpolation of the atmospheric state);
  - scatter (node-major): `F_n = Σ_k sweight[n,k] F_k` with
    `sweight = B̃_{kn} / (Jw)_n`, where `B̃_{kn} ≥ 0` is the area of `Ω_k`
    inside node `n`'s control volume, whose area is `(Jw)_n`; follow with
    `weighted_dss!`;
  - transpose of the scatter (polygon-major, on the gather's sparsity):
    `ḡ_k = Σ_n tweight[k,n] g_n` with `tweight = B̃_{kn} / b_area_k`, for the
    fluxes the atmosphere hands down (radiation, precipitation). For any
    polygon weights `w`, `Σ_k b_area_k w_k ḡ_k = Σ_n (Jw)_n (Σ_k sweight[n,k] w_k) g_n`;
  - `node_cov[n] = Σ_{k wet} sweight[n,k]`: the fraction of node `n`'s area
    over wet ocean.

Flat index conventions, over the elements and cells of this process: SE node
`n = (e-1) Nq² + (j-1) Nq + i` (i fastest within each element, then the
elements this process stores — the order of
`ConservativeRegriddingClimaCoreExt.flat_nodal_data` /
`Fields.field2array` on ClimaCore 0.15 `(v, i, j, h)` / `VIJFH` storage);
FV cells by the slot of the [`CellExchange`](@ref) `cells`, whose stored slots
`c = (j-1) Nx + i` follow `vec(OC.interior(field, :, :, Nz))` of the local grid.
On one process the slots are the cells.

# Fields
- `n_poly`, `n_nodes`, `n_elem`, `n_oc`: entity counts; `n_oc` counts the stored cells
- `gpoly_ptr`, `gnode`, `gweight`: polygon-major CSR of the node→polygon gather
- `tweight`: transpose weights `B̃_{kn} / b_area_k` on the gather's CSR
- `snode_ptr`, `spoly`, `sweight`: node-major CSR of the polygon→node scatter
- `node_cov`: nodal wet coverage fractions
- `elem_of_poly`, `oc_of_poly`: SE element / FV cell slot owning each polygon
- `soc_ptr`, `soc_poly`: slot-major CSR of the polygon→cell scatter
- `area`: geometric polygon areas [m²]
- `b_area`: SE-side scatter areas `Σ_n B̃_{kn} = Σ_n sweight[n,k] (Jw)_n` [m²]; equals
  `area` up to quadrature error and makes SE-side conservation statements
  exact
- `oc_wet_area`: retained (wet) polygon area, of the polygons of every process, per
  stored FV cell [m²]; zero for dry (immersed) cells and tripolar fold shadow cells
- `oc_area`: the FV grid's own cell area `Az` [m²] where `oc_wet_area > 0`;
  zero elsewhere
- `cells`: the [`CellExchange`](@ref) of the FV cells
"""
struct ExchangeGrid{FT, VI <: AbstractVector{Int32}, VF <: AbstractVector{FT}, CX}
    n_poly::Int
    n_nodes::Int
    n_elem::Int
    n_oc::Int
    gpoly_ptr::VI
    gnode::VI
    gweight::VF
    tweight::VF
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
    cells::CX
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
    single_process_space(boundary_space)

A space on the CPU of one process with the whole mesh, element order and
quadrature of `boundary_space`, whose element `e` is global element `e` of
`boundary_space`. Built rather than adapted: adapting a ClimaCore space to
`Array` gives it the default comms context, which is MPI in an MPI job.
"""
function single_process_space(boundary_space)
    grid = CC.Spaces.grid(boundary_space)
    topology = grid.topology
    cpu_topology = CC.Topologies.Topology2D(
        ClimaComms.SingletonCommsContext(ClimaComms.CPUSingleThreaded()),
        topology.mesh,
        topology.elemorder,
        nothing,
        topology.orderindex,
    )
    cpu_grid = CC.Grids.SpectralElementGrid2D(
        cpu_topology,
        grid.quadrature_style;
        enable_bubble = grid.enable_bubble,
        discretization = grid.discretization,
    )
    return CC.Spaces.SpectralElementSpace2D(cpu_grid)
end

"""
    build_exchange_grid(boundary_space, grid_oc)

Construct the [`ExchangeGrid`](@ref) of this process between a ClimaCore
cubed-sphere `boundary_space` and an Oceananigans `grid_oc`, on the CPU in
`Float64`: intersect the SE elements of the whole mesh with the FV cells of the
whole grid (fold-aware on a `TripolarGrid`, so fold shadow cells never produce
polygons), keep the polygons in the elements this process stores over wet
(not immersed) cells, then integrate the SEM basis (`accumulate_principled_b`,
gather weights) and the control volumes (`control_volume_weights`, scatter and
transpose weights and nodal coverages) over each polygon. Collective over the
processes of `boundary_space`. Move the result to the device with
[`on_device`](@ref).
"""
function build_exchange_grid(boundary_space, grid_oc)
    CRExt = get_ConservativeRegriddingCCExt()
    @assert !isnothing(CRExt) "ConservativeRegriddingClimaCoreExt must be loaded"
    GO = CRExt.GO

    context = ClimaComms.context(boundary_space)
    mesh_space_cpu = single_process_space(boundary_space)
    global_grid_oc_cpu = OC.on_architecture(
        OC.CPU(),
        OC.DistributedComputations.reconstruct_global_grid(grid_oc),
    )
    grid_oc_underlying_cpu = underlying_grid(global_grid_oc_cpu)

    FT = CC.Spaces.undertype(mesh_space_cpu)
    R = Float64(CC.Spaces.topology(mesh_space_cpu).mesh.domain.radius)
    manifold = CR.Spherical(; radius = R)

    # 1. SE-element × FV-cell intersection polygons.
    dst_tree = CR.Trees.treeify(manifold, mesh_space_cpu)
    src_tree = CR.Trees.treeify(manifold, grid_oc_underlying_cpu)
    intersections = CR.intersection_areas(
        manifold,
        CR.False(),
        dst_tree,
        src_tree;
        intersection_operator = CR.IntersectionGridOperator(manifold),
    )
    elem_of_poly, oc_of_poly, polys = SparseArrays.findnz(intersections)
    n_oc = size(intersections, 2)

    # 2. Keep the polygons in the elements of this process over wet cells.
    local_elems = CC.Spaces.topology(boundary_space).local_elem_gidx
    keep = [e in local_elems for e in elem_of_poly]
    if global_grid_oc_cpu isa OC.ImmersedBoundaryGrid
        Nx_oc, _, Nz_oc = size(global_grid_oc_cpu)
        for k in eachindex(keep)
            c = oc_of_poly[k]
            i, j = mod1(c, Nx_oc), (c - 1) ÷ Nx_oc + 1
            keep[k] &= !OC.ImmersedBoundaries.immersed_cell(i, j, Nz_oc, global_grid_oc_cpu)
        end
    end
    elem_of_poly, oc_of_poly, polys = elem_of_poly[keep], oc_of_poly[keep], polys[keep]
    area = [GO.area(manifold, poly) for poly in polys]

    # 3. Per-polygon weights: the SEM basis integrals B_kn for the gather and
    #    the control-volume weights B̃_kn for the scatter and the coverages.
    #    Nodes are numbered over the elements of this process.
    qs = CC.Spaces.quadrature_style(mesh_space_cpu)
    Nq = CC.Quadratures.degrees_of_freedom(qs)
    triangle_quad_degree = 2 * (Nq - 1)
    Jw = Array(CRExt.se_node_weights(boundary_space)) # flat nodal, this process
    n_nodes = length(Jw)
    volumes = ControlVolumes(mesh_space_cpu)

    g_poly, g_node, g_b, g_t = Int[], Int[], Float64[], Float64[] # gather COO: B_kn, B̃_kn
    s_poly, s_node, s_w = Int[], Int[], Float64[] # scatter COO: B̃_kn / A_n
    b_sum = zeros(Float64, length(area))
    b_area = zeros(Float64, length(area))
    for k in eachindex(area)
        elem = elem_of_poly[k]
        B = CRExt.accumulate_principled_b(
            manifold,
            mesh_space_cpu,
            elem,
            polys[k];
            triangle_quad_degree,
        )
        B̃ = control_volume_weights(volumes, mesh_space_cpu, elem, polys[k])
        node_offset = (elem - first(local_elems)) * Nq^2
        for j in 1:Nq, i in 1:Nq
            n = node_offset + (j - 1) * Nq + i
            w = B̃[i, j] / Jw[n]
            if w != 0
                push!(s_poly, k)
                push!(s_node, n)
                push!(s_w, w)
            end
            Bij = B[i, j]
            if Bij != 0 || B̃[i, j] != 0
                push!(g_poly, k)
                push!(g_node, n)
                push!(g_b, Bij)
                push!(g_t, B̃[i, j])
                b_sum[k] += Bij
                b_area[k] += B̃[i, j]
            end
        end
    end
    elem_of_poly = elem_of_poly .- (first(local_elems) - 1)

    # Degenerate polygons (b_sum ≤ 0 or b_area ≤ 0 can only come from
    # round-off on slivers) cannot be normalized for the gathers; drop their
    # COO entries in both directions and zero their area so they are inert in
    # every reduction.
    degenerate = findall((b_sum .<= 0) .| (b_area .<= 0))
    if !isempty(degenerate)
        is_degenerate = falses(length(area))
        is_degenerate[degenerate] .= true
        g_keep = .!is_degenerate[g_poly]
        g_poly, g_node, g_b, g_t = g_poly[g_keep], g_node[g_keep], g_b[g_keep], g_t[g_keep]
        s_keep = .!is_degenerate[s_poly]
        s_poly, s_node, s_w = s_poly[s_keep], s_node[s_keep], s_w[s_keep]
        area[degenerate] .= 0
        b_sum[degenerate] .= 0
        b_area[degenerate] .= 0
    end

    # 4. CSR assembly for both directions.
    n_poly = length(area)
    gpoly_ptr, gperm = _build_csr(g_poly, n_poly)
    gnode = Int32.(g_node[gperm])
    gweight = [g_b[p] / b_sum[g_poly[p]] for p in gperm]
    tweight = [g_t[p] / b_area[g_poly[p]] for p in gperm]

    snode_ptr, sperm = _build_csr(s_node, n_nodes)
    spoly = Int32.(s_poly[sperm])
    sweight = s_w[sperm]

    node_cov = zeros(Float64, n_nodes)
    for p in eachindex(s_node)
        node_cov[s_node[p]] += s_w[p]
    end

    # 5. The cells: a slot for each cell this process stores or its polygons use,
    #    and the wet area and the grid area of the stored cells, the wet area
    #    summed over the polygons of every process.
    slot, owned, plan = plan_cell_exchange(context, grid_oc, n_oc, oc_of_poly)
    n_local = length(owned)
    oc_of_poly = [slot[c] for c in oc_of_poly]
    soc_ptr, soc_poly = _build_csr(oc_of_poly, length(slot))
    oc_wet_area = zeros(Float64, n_local)
    _scatter_polys_to_cells!(
        oc_wet_area,
        CellExchange(Float64, context, plan, n_local),
        soc_ptr,
        soc_poly,
        ones(n_poly),
        area,
        ones(n_local),
    )
    Nx_oc, _, Nz_oc = size(grid_oc_underlying_cpu)
    oc_area = zeros(Float64, n_local)
    for s in findall(>(0), oc_wet_area)
        c = owned[s]
        i, j = mod1(c, Nx_oc), (c - 1) ÷ Nx_oc + 1
        oc_area[s] = OC.Operators.Azᶜᶜᶜ(i, j, Nz_oc, grid_oc_underlying_cpu)
    end
    cells = CellExchange(FT, context, plan, n_local)

    return ExchangeGrid{FT, Vector{Int32}, Vector{FT}, typeof(cells)}(
        n_poly,
        n_nodes,
        length(local_elems),
        n_local,
        gpoly_ptr,
        gnode,
        FT.(gweight),
        FT.(tweight),
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
        cells,
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
    gather_nodes_to_polys!(poly_values, eg::ExchangeGrid, nodal_values, weight = eg.gweight)

Average a flat SE nodal vector onto each polygon: `f̄_k = Σ_n weight[k,n] f_n`.
Rows sum to 1, so constants are preserved exactly. The default `gweight`
interpolates the atmospheric state to high order; `eg.tweight` is the transpose of
the scatter weights.
"""
gather_nodes_to_polys!(poly_values, eg::ExchangeGrid, nodal_values, weight = eg.gweight) =
    _csr_matvec!(poly_values, eg.gpoly_ptr, eg.gnode, weight, nodal_values)

@kernel function _gather_cells_kernel!(dst, oc_of_poly, src, remote, n_local)
    k = @index(Global)
    @inbounds begin
        s = oc_of_poly[k]
        dst[k] = s <= n_local ? src[s] : remote[s - n_local]
    end
end

"""
    gather_cells_to_polys!(poly_values, eg::ExchangeGrid, cell_values)

Copy the owning FV cell's value onto each polygon (direct indexing; each
polygon lies inside exactly one cell). `cell_values` holds the cells this
process stores; the others come from their processes.
"""
function gather_cells_to_polys!(poly_values, eg::ExchangeGrid, cell_values)
    remote = pull_cells!(eg.cells, cell_values)
    backend = KernelAbstractions.get_backend(poly_values)
    launch_kernel!(
        _gather_cells_kernel!,
        backend,
        eg.n_poly,
        poly_values,
        eg.oc_of_poly,
        cell_values,
        remote,
        eg.n_oc,
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

# dst[r] = Σ_{k ∈ slot n_local + r} area_k src_k: the partial sums of the cells of other processes
@kernel function _remote_partials_kernel!(dst, ptr, polys, area, src, n_local)
    r = @index(Global)
    acc = zero(eltype(dst))
    @inbounds begin
        s = n_local + r
        for p in ptr[s]:(ptr[s + 1] - 1)
            k = polys[p]
            acc += area[k] * src[k]
        end
        dst[r] = acc
    end
end

@kernel function _scatter_cells_kernel!(
    dst,
    ptr,
    polys,
    area,
    cell_area,
    src,
    returned_ptr,
    returned_pos,
    returned,
)
    c = @index(Global)
    @inbounds begin
        A = cell_area[c]
        if A > 0
            acc = zero(eltype(dst))
            for p in ptr[c]:(ptr[c + 1] - 1)
                k = polys[p]
                acc += area[k] * src[k]
            end
            for p in returned_ptr[c]:(returned_ptr[c + 1] - 1)
                acc += returned[returned_pos[p]]
            end
            dst[c] = acc / A
        else
            dst[c] = 0
        end
    end
end

# Push the partial sums of the cells of other processes, then add the partial
# sums returned for the stored cells, in process order, to their own.
function _scatter_polys_to_cells!(
    cell_values,
    cx::CellExchange,
    soc_ptr,
    soc_poly,
    poly_values,
    poly_area,
    cell_area,
)
    backend = KernelAbstractions.get_backend(cell_values)
    isempty(cx.remote_partial) || launch_kernel!(
        _remote_partials_kernel!,
        backend,
        length(cx.remote_partial),
        cx.remote_partial,
        soc_ptr,
        soc_poly,
        poly_area,
        poly_values,
        cx.n_local,
    )
    returned = push_partials!(cx)
    launch_kernel!(
        _scatter_cells_kernel!,
        backend,
        cx.n_local,
        cell_values,
        soc_ptr,
        soc_poly,
        poly_area,
        cell_area,
        poly_values,
        cx.returned_ptr,
        cx.returned_pos,
        returned,
    )
    return cell_values
end

"""
    scatter_polys_to_cells!(cell_values, eg::ExchangeGrid, poly_values,
                            poly_area = eg.b_area, cell_area = eg.oc_area)

Per-polygon values onto each FV cell this process stores, summed over the
polygons of every process: `F_c = Σ_{k ∈ c} poly_area_k F_k / cell_area_c`.
With the defaults, `F_c` is a flux per unit FV cell area whose integral `Σ_c Az_c F_c`
equals the SE-side integral `Σ_k b_area_k F_k`; pass `eg.area, eg.oc_wet_area` for the
area-weighted mean of a state. Cells with zero `cell_area` (dry, fold shadows) are set
to 0; run [`mirror_fold_partners!`](@ref) afterwards to fill the shadow copies.
"""
scatter_polys_to_cells!(
    cell_values,
    eg::ExchangeGrid,
    poly_values,
    poly_area = eg.b_area,
    cell_area = eg.oc_area,
) = _scatter_polys_to_cells!(
    cell_values,
    eg.cells,
    eg.soc_ptr,
    eg.soc_poly,
    poly_values,
    poly_area,
    cell_area,
)

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
