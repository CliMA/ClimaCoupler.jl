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
    PolygonIntersectionOperator{M <: CR.Manifold}

`ConservativeRegridding` intersection operator that assembles a sparse matrix
of intersection *polygons* rather than scalar areas.
"""
struct PolygonIntersectionOperator{M <: CR.Manifold}
    manifold::M
end

function _intersection_polygon_type()
    CRExt = get_ConservativeRegriddingCCExt()
    @assert !isnothing(CRExt)
    GI = CRExt.GI
    GO = CRExt.GO
    return GI.Polygon{
        true,
        false,
        Vector{
            GI.LinearRing{
                true,
                false,
                Vector{GO.UnitSpherical.UnitSphericalPoint{Float64}},
                Nothing,
                Nothing,
            },
        },
        Nothing,
        Nothing,
    }
end

CR.IntersectionReturnStyle(::PolygonIntersectionOperator) = CR.OutOfPlaceSingleResult()
CR.output_eltype(::PolygonIntersectionOperator) = _intersection_polygon_type()
CR.output_eltype(op::PolygonIntersectionOperator, src_tree, dst_tree) =
    CR.output_eltype(op)

function CR.should_store_result(::PolygonIntersectionOperator, result)
    result === nothing && return false
    CRExt = get_ConservativeRegriddingCCExt()
    return result isa CRExt.GI.Polygon
end

function (op::PolygonIntersectionOperator)(src_cell, dst_cell)
    CRExt = get_ConservativeRegriddingCCExt()
    GO = CRExt.GO
    GI = CRExt.GI
    intersection_poly = GO.intersection(
        GO.ConvexConvexSutherlandHodgman(op.manifold),
        src_cell,
        dst_cell;
        target = GI.PolygonTrait(),
    )
    iszero(GO.area(op.manifold, intersection_poly)) && return nothing
    return intersection_poly
end

"""
    ExchangeGrid{FT, VI, VF}

Sparse coupling between the SE boundary space, the FV (Oceananigans) surface
cells, and the intersection polygons that tile their overlap, stored in
compressed-sparse-row (CSR) form for the direction each coupling is consumed
in, so every per-step gather/scatter is a race-free segmented reduction.

SE-side weights come from the SEM basis integrals
``B_{kn} = ∫_{Ω_k} ϕ_n \\, dA`` over each polygon `Ω_k`
(`ConservativeRegriddingClimaCoreExt.accumulate_principled_b`):

  - gather (polygon-major): `f̄_k = Σ_n gweight[k,n] f_n` with
    `gweight = B_{kn} / Σ_n B_{kn}` (rows sum to 1; constants preserved);
  - scatter (node-major): `F_n = Σ_k sweight[n,k] F_k` with
    `sweight = B_{kn} / (Jw)_n` — the per-element L2 projection, to be
    followed by `weighted_dss!`;
  - `node_cov[n] = Σ_{k wet} B_{kn} / (Jw)_n` — the wet-ocean coverage of
    node `n`: scattering a constant yields `constant × node_cov`. It carries
    Gibbs overshoot and systematic geometric deficits, so the quantity to use
    as a wet fraction is `node_cov / node_cov_total` (clamped to [0, 1]),
    where `node_cov_total` is the same sum over *all* polygons before the
    wet-mask filter — the ratio cancels the geometric error.

Flat index conventions: SE node `n = (e-1) Nq² + (j-1) Nq + i` (the `IJFH`
layout flattened by `data2array`); FV cell `c = (j-1) Nx + i`
(`vec(OC.interior(field, :, :, Nz))` and the CR column index).

# Fields
- `n_poly`, `n_nodes`, `n_elem`, `n_oc`: entity counts
- `gpoly_ptr`, `gnode`, `gweight`: polygon-major CSR of the node→polygon gather
- `snode_ptr`, `spoly`, `sweight`: node-major CSR of the polygon→node scatter
- `node_cov`, `node_cov_total`: nodal wet / geometric coverage
- `elem_of_poly`, `oc_of_poly`: SE element / FV cell owning each polygon
- `soc_ptr`, `soc_poly`: FV-cell-major CSR of the polygon→cell scatter
- `area`: geometric polygon areas [m²]
- `b_area`: SE-side quadrature areas `Σ_n B_{kn}` [m²]; equals `area` up to
  quadrature error and makes SE-side conservation statements exact
- `oc_wet_area`: retained (wet) polygon area per FV cell [m²]; zero for dry
  (immersed) cells and tripolar fold shadow cells
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
    node_cov_total::VF
    elem_of_poly::VI
    oc_of_poly::VI
    soc_ptr::VI
    soc_poly::VI
    area::VF
    b_area::VF
    oc_wet_area::VF
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
    build_exchange_grid(boundary_space, grid_oc; sliver_rtol = 0.0)

Construct an [`ExchangeGrid`](@ref) between a ClimaCore cubed-sphere
`boundary_space` and an Oceananigans `grid_oc`, on the CPU in `Float64`:
intersect SE elements with FV cells (fold-aware on a `TripolarGrid`, so fold
shadow cells never produce polygons), drop polygons over dry (immersed) cells
and slivers with `area < sliver_rtol * mean(area)` (guards `Float32` weight
underflow; disabled by default), then integrate the SEM basis over each
polygon (`accumulate_principled_b`) to obtain the gather/scatter weights and
nodal coverages. Move the result to the device with [`on_device`](@ref).
"""
function build_exchange_grid(boundary_space, grid_oc; sliver_rtol = 0.0)
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
        intersection_operator = PolygonIntersectionOperator(manifold),
    )
    elem_of_poly, oc_of_poly, polys = SparseArrays.findnz(intersections)
    n_elem, n_oc = size(intersections)
    area = [GO.area(manifold, poly) for poly in polys]

    # 2. Mark polygons over dry (immersed) surface cells and slivers. They are
    #    dropped from the exchange grid but still contribute to the geometric
    #    coverage `node_cov_total`.
    keep = trues(length(area))
    if grid_oc isa OC.ImmersedBoundaryGrid
        grid_with_mask_cpu = OC.on_architecture(OC.CPU(), grid_oc)
        Nx_oc, Ny_oc, Nz_oc = size(grid_with_mask_cpu)
        for k in eachindex(keep)
            c = oc_of_poly[k]
            i, j = mod1(c, Nx_oc), (c - 1) ÷ Nx_oc + 1
            keep[k] =
                !OC.ImmersedBoundaries.immersed_cell(i, j, Nz_oc, grid_with_mask_cpu)
        end
    end
    if sliver_rtol > 0
        cutoff = sliver_rtol * sum(area) / length(area)
        keep .&= area .> cutoff
    end

    # 3. SEM basis integrals B_kn over each polygon. All polygons contribute
    #    to node_cov_total; only kept (wet) polygons produce COO weights,
    #    indexed by their position in the kept numbering.
    qs = CC.Spaces.quadrature_style(boundary_space_cpu)
    Nq = CC.Quadratures.degrees_of_freedom(qs)
    triangle_quad_degree = 2 * (Nq - 1)
    Jw = CRExt.se_node_weights(boundary_space_cpu) # Float64, flat nodal
    n_nodes = length(Jw)

    coo_poly = Int[]
    coo_node = Int[]
    coo_b = Float64[]
    node_cov_total = zeros(Float64, n_nodes)
    b_sum = zeros(Float64, count(keep))
    kept_id = 0
    for k in eachindex(area)
        elem = elem_of_poly[k]
        B = CRExt.accumulate_principled_b(
            manifold,
            boundary_space_cpu,
            elem,
            polys[k];
            triangle_quad_degree,
        )
        kept = keep[k]
        kept && (kept_id += 1)
        node_offset = (elem - 1) * Nq^2
        for j in 1:Nq, i in 1:Nq
            Bij = B[i, j]
            Bij == 0 && continue
            n = node_offset + (j - 1) * Nq + i
            node_cov_total[n] += Bij / Jw[n]
            if kept
                push!(coo_poly, kept_id)
                push!(coo_node, n)
                push!(coo_b, Bij)
                b_sum[kept_id] += Bij
            end
        end
    end
    elem_of_poly, oc_of_poly, area = elem_of_poly[keep], oc_of_poly[keep], area[keep]

    # Degenerate polygons (b_sum ≤ 0 can only come from quadrature round-off
    # on slivers) cannot be normalized; drop their COO entries and zero their
    # area so they are inert in every reduction.
    degenerate = findall(<=(0), b_sum)
    if !isempty(degenerate)
        is_degenerate = falses(length(area))
        is_degenerate[degenerate] .= true
        coo_keep = .!is_degenerate[coo_poly]
        coo_poly, coo_node, coo_b = coo_poly[coo_keep], coo_node[coo_keep], coo_b[coo_keep]
        area[degenerate] .= 0
        b_sum[degenerate] .= 0
    end

    # 4. CSR assembly for both directions.
    n_poly = length(area)
    gpoly_ptr, gperm = _build_csr(coo_poly, n_poly)
    gnode = Int32.(coo_node[gperm])
    gweight = [coo_b[p] / b_sum[coo_poly[p]] for p in gperm]

    snode_ptr, sperm = _build_csr(coo_node, n_nodes)
    spoly = Int32.(coo_poly[sperm])
    sweight = [coo_b[p] / Jw[coo_node[p]] for p in sperm]

    node_cov = zeros(Float64, n_nodes)
    for p in eachindex(coo_node)
        node_cov[coo_node[p]] += coo_b[p] / Jw[coo_node[p]]
    end

    soc_ptr, soc_poly = _build_csr(Vector{Int}(oc_of_poly), n_oc)
    oc_wet_area = zeros(Float64, n_oc)
    for k in eachindex(area)
        oc_wet_area[oc_of_poly[k]] += area[k]
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
        FT.(node_cov_total),
        Int32.(elem_of_poly),
        Int32.(oc_of_poly),
        soc_ptr,
        soc_poly,
        FT.(area),
        FT.(b_sum),
        FT.(oc_wet_area),
    )
end

#=
# Gather/scatter operations

Segmented reductions over the CSR structures above: race-free, deterministic,
allocation-free. Each wrapper runs a serial loop on the CPU or a
KernelAbstractions kernel on the GPU, selected by the destination's backend.
=#

import KernelAbstractions

const _KA_WORKGROUP = 256

@inline _is_cpu(x) = KernelAbstractions.get_backend(x) isa KernelAbstractions.CPU

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
    if _is_cpu(dst)
        @inbounds for r in eachindex(dst)
            acc = zero(eltype(dst))
            for p in ptr[r]:(ptr[r + 1] - 1)
                acc += w[p] * src[col[p]]
            end
            dst[r] = acc
        end
    else
        backend = KernelAbstractions.get_backend(dst)
        _csr_matvec_kernel!(backend, _KA_WORKGROUP)(dst, ptr, col, w, src; ndrange = length(dst))
    end
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
    if _is_cpu(poly_values)
        @inbounds for k in eachindex(poly_values)
            poly_values[k] = cell_values[eg.oc_of_poly[k]]
        end
    else
        backend = KernelAbstractions.get_backend(poly_values)
        _gather_cells_kernel!(backend, _KA_WORKGROUP)(
            poly_values,
            eg.oc_of_poly,
            cell_values;
            ndrange = eg.n_poly,
        )
    end
    return poly_values
end

"""
    scatter_polys_to_nodes!(nodal_values, eg::ExchangeGrid, poly_values)

Project per-polygon values onto SE nodes via the per-element L2 projection:
`F_n = Σ_k sweight F_k`. The result is *coverage-weighted*: scattering a
constant yields `constant × node_cov`. Follow with `weighted_dss!` on the
receiving field. For a per-unit-wet-area result use
[`scatter_polys_to_nodes_normalized!`](@ref).
"""
scatter_polys_to_nodes!(nodal_values, eg::ExchangeGrid, poly_values) =
    _csr_matvec!(nodal_values, eg.snode_ptr, eg.spoly, eg.sweight, poly_values)

@kernel function _csr_matvec_normalized_kernel!(dst, ptr, col, w, src, cov, cutoff)
    r = @index(Global)
    @inbounds begin
        c = cov[r]
        if c > cutoff
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
    scatter_polys_to_nodes_normalized!(nodal_values, eg::ExchangeGrid, poly_values,
                                       cov_cutoff)

Like [`scatter_polys_to_nodes!`](@ref) but normalized by the nodal wet
coverage, yielding a per-unit-wet-area value: constants are reproduced exactly
wherever `node_cov > cov_cutoff`; nodes at or below the cutoff are set to 0.
"""
function scatter_polys_to_nodes_normalized!(
    nodal_values,
    eg::ExchangeGrid,
    poly_values,
    cov_cutoff,
)
    if _is_cpu(nodal_values)
        @inbounds for n in eachindex(nodal_values)
            c = eg.node_cov[n]
            if c > cov_cutoff
                acc = zero(eltype(nodal_values))
                for p in eg.snode_ptr[n]:(eg.snode_ptr[n + 1] - 1)
                    acc += eg.sweight[p] * poly_values[eg.spoly[p]]
                end
                nodal_values[n] = acc / c
            else
                nodal_values[n] = 0
            end
        end
    else
        backend = KernelAbstractions.get_backend(nodal_values)
        _csr_matvec_normalized_kernel!(backend, _KA_WORKGROUP)(
            nodal_values,
            eg.snode_ptr,
            eg.spoly,
            eg.sweight,
            poly_values,
            eg.node_cov,
            cov_cutoff;
            ndrange = eg.n_nodes,
        )
    end
    return nodal_values
end

@kernel function _scatter_cells_kernel!(dst, ptr, polys, area, wet_area, src)
    c = @index(Global)
    @inbounds begin
        aw = wet_area[c]
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
    scatter_polys_to_cells!(cell_values, eg::ExchangeGrid, poly_values)

Area-weighted average of per-polygon values onto each FV cell:
`F_c = Σ_{k ∈ c} area_k F_k / oc_wet_area_c`, conserving the area integral
exactly. Cells with zero wet area (dry, fold shadows) are set to 0; run
[`mirror_fold_partners!`](@ref) afterwards to fill the shadow copies.
"""
function scatter_polys_to_cells!(cell_values, eg::ExchangeGrid, poly_values)
    if _is_cpu(cell_values)
        @inbounds for c in eachindex(cell_values)
            aw = eg.oc_wet_area[c]
            if aw > 0
                acc = zero(eltype(cell_values))
                for p in eg.soc_ptr[c]:(eg.soc_ptr[c + 1] - 1)
                    k = eg.soc_poly[p]
                    acc += eg.area[k] * poly_values[k]
                end
                cell_values[c] = acc / aw
            else
                cell_values[c] = 0
            end
        end
    else
        backend = KernelAbstractions.get_backend(cell_values)
        _scatter_cells_kernel!(backend, _KA_WORKGROUP)(
            cell_values,
            eg.soc_ptr,
            eg.soc_poly,
            eg.area,
            eg.oc_wet_area,
            poly_values;
            ndrange = eg.n_oc,
        )
    end
    return cell_values
end

#=
# Tripolar fold mirroring

On a `RightCenterFolded` `TripolarGrid` the fold-row cell `(i, Ny)` is the
same physical cell as `(Nx + 1 - i, Ny)`; the exchange grid keeps one copy as
a real cell and leaves the other a degenerate shadow with no polygons, so
shadow slots hold 0 after a scatter. Mirror each primary (fold locals
`1..Nx÷4` and `Nx÷2 + 1 .. Nx÷2 + Nx÷4`; partner of local `r` is
`Nx + 1 - r`) into its partner slot. Replicates the internal
`mirror_fold_partners!` in `ConservativeRegriddingOceananigansExt`.
=#

@kernel function _mirror_fold_kernel!(dst, Nx, Nh, Nquarter, fold_offset)
    i = @index(Global)
    r = ifelse(i <= Nquarter, i, i + Nh - Nquarter)
    @inbounds dst[Nx + 1 - r + fold_offset] = dst[r + fold_offset]
end

"""
    mirror_fold_partners!(cell_values, grid)

Copy each fold-row primary cell's value into its shadow partner slot on a
`RightCenterFolded` tripolar grid; no-op for other grids. `cell_values` is a
flat vector over the `Nx × Ny` surface cells.
"""
mirror_fold_partners!(cell_values, grid) = cell_values

mirror_fold_partners!(cell_values, grid::OC.ImmersedBoundaryGrid) =
    mirror_fold_partners!(cell_values, grid.underlying_grid)

function mirror_fold_partners!(
    cell_values,
    grid::OC.Grids.ZRegOrthogonalSphericalShellGrid{<:Number, <:Any, OC.RightCenterFolded},
)
    Nx, Ny, _ = size(grid)
    Nh = Nx ÷ 2
    Nquarter = Nx ÷ 4
    fold_offset = (Ny - 1) * Nx
    if _is_cpu(cell_values)
        @inbounds for i in 1:(2 * Nquarter)
            r = ifelse(i <= Nquarter, i, i + Nh - Nquarter)
            cell_values[Nx + 1 - r + fold_offset] = cell_values[r + fold_offset]
        end
    else
        backend = KernelAbstractions.get_backend(cell_values)
        _mirror_fold_kernel!(backend, _KA_WORKGROUP)(
            cell_values,
            Nx,
            Nh,
            Nquarter,
            fold_offset;
            ndrange = 2 * Nquarter,
        )
    end
    return cell_values
end

#=
# Wet-ocean surface fraction
=#

"""
    wet_ocean_fraction_field(boundary_space, eg::ExchangeGrid;
                             topography_damping_factor = 5)

Build the wet-ocean surface fraction as a `CC.Fields.Field` on
`boundary_space` from a CPU-resident [`ExchangeGrid`](@ref): the nodal ratio
`node_cov / node_cov_total` clamped to [0, 1], made continuous with
`weighted_dss!`, then smoothed with the same diffusion recipe ClimaAtmos
applies to its orography (`κ = 0.05 Δh²`, `dt = 1`,
`maxiter = round(log(damping_factor)/0.05)`; see
`ClimaAtmos.make_hybrid_spaces`), so the coupler never sees land-sea
contrasts sharper than the atmosphere's smoothed topography.
`topography_damping_factor` must match the ClimaAtmos option of the same
name. The complement is the land fraction; sea ice and open ocean partition
the wet fraction itself (see `FieldExchanger.align_surface_fractions!`).
"""
function wet_ocean_fraction_field(
    boundary_space,
    eg::ExchangeGrid{FT, Vector{Int32}};
    topography_damping_factor = 5,
) where {FT}
    frac_nodal = similar(eg.node_cov)
    @. frac_nodal = ifelse(
        eg.node_cov_total > 0,
        clamp(eg.node_cov / eg.node_cov_total, FT(0), FT(1)),
        FT(0),
    )

    CRExt = get_ConservativeRegriddingCCExt()
    field = CC.Fields.zeros(boundary_space)
    device_array_type = ClimaComms.array_type(ClimaComms.device(boundary_space))
    CRExt.vec_to_se_field!(field, Adapt.adapt(device_array_type, frac_nodal))

    # Reconcile shared element-boundary nodes, then filter exactly like the
    # atmosphere's Earth orography.
    dss_buffer = Utilities.init_dss_buffer(field)
    Utilities.apply_dss!(field, dss_buffer)
    @. field = clamp(field, FT(0), FT(1))
    Δh = CC.Spaces.node_horizontal_length_scale(boundary_space)
    maxiter = round(Int, log(topography_damping_factor) / 0.05)
    CC.Hypsography.diffuse_surface_elevation!(
        field;
        κ = FT(0.05) * FT(Δh)^2,
        dt = FT(1),
        maxiter,
    )
    @. field = clamp(field, FT(0), FT(1))
    return field
end
