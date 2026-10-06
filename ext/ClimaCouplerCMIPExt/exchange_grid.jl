#=
# Exchange (intersection) grid

The spectral-element (SE) boundary space is treated as a finite-volume mesh of
*node cells*: GLL node `(i, j)` of an element owns the image on the sphere of
the reference box `[ζ_{i-1}, ζ_i] × [ζ_{j-1}, ζ_j]`, with `ζ_0 = -1` and
`ζ_i = ζ_{i-1} + w_i` (1D widths equal to the GLL weights). Lines of constant
reference coordinate map to great circles under the cubed-sphere element maps,
so node cells are spherical quadrilaterals that tile the sphere.

Intersecting the node cells with the wet cells of an Oceananigans grid gives
the exchange-grid polygons. Each polygon lies in exactly one node cell and one
ocean cell. Turbulent fluxes are computed per polygon and returned to both
sides by area-weighted sums; the ocean and sea-ice area fractions are the
shares of each node cell covered by those polygons. Geometry is built once on
the CPU in `Float64`, cast to the simulation float type, and moved to the
device, where every per-step operation is a sparse gather/scatter.
=#

"""
    ExchangeGrid{FT, VI, VF}

Sparse coupling between SE nodes, ocean cells, and the polygons where their
cells overlap, stored as compressed-sparse-row (CSR) structures in the
direction each is used.

Flat indices: SE node `n = (e - 1) Nq² + (j - 1) Nq + i` (element-local storage
order of `Fields.field2array`); ocean cell `c = (j - 1) Nx + i`.

# Fields
- `n_poly`, `n_nodes`, `n_elem`, `n_oc`: entity counts
- `node_of_poly`, `elem_of_poly`, `oc_of_poly`: SE node, SE element and ocean
  cell containing each polygon
- `snode_ptr`, `spoly`, `sweight`: node-major CSR of the polygon → node
  scatter; `sweight` is the share of the node cell a polygon covers
- `node_area`: area of each node cell [m²]
- `wet_share`: share of each node cell covered by wet ocean (before DSS)
- `soc_ptr`, `soc_poly`: ocean-cell-major CSR of the polygon → cell scatter
- `area`: polygon areas [m²]
- `oc_wet_area`: polygon area per ocean cell [m²]; zero for dry cells and
  tripolar fold shadow cells
"""
struct ExchangeGrid{FT, VI <: AbstractVector{Int32}, VF <: AbstractVector{FT}}
    n_poly::Int
    n_nodes::Int
    n_elem::Int
    n_oc::Int
    node_of_poly::VI
    elem_of_poly::VI
    oc_of_poly::VI
    snode_ptr::VI
    spoly::VI
    sweight::VF
    node_area::VF
    wet_share::VF
    soc_ptr::VI
    soc_poly::VI
    area::VF
    oc_wet_area::VF
end

Adapt.@adapt_structure ExchangeGrid

function Base.show(io::IO, eg::ExchangeGrid{FT}) where {FT}
    print(
        io,
        "ExchangeGrid{$FT}: $(eg.n_poly) polygons from $(eg.n_nodes) SE node cells × $(eg.n_oc) FV cells",
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

#=
## Node cells
=#

# Unit-sphere corner of the node-cell mesh of one cube face, at corner index
# `(I₀, J₀) ∈ 0:ne·Nq` along the face; `ζ` are the GLL subcell breakpoints.
function _node_cell_corner(CRExt, mesh, ζ, Nq, face, I₀, J₀)
    ie = min(I₀ ÷ Nq + 1, mesh.ne)
    je = min(J₀ ÷ Nq + 1, mesh.ne)
    ξ = ζ[I₀ - (ie - 1) * Nq + 1]
    η = ζ[J₀ - (je - 1) * Nq + 1]
    x = CC.Geometry.components(
        CC.Meshes.coordinates(mesh, CartesianIndex(ie, je, face), (ξ, η)),
    )
    r = sqrt(x[1]^2 + x[2]^2 + x[3]^2)
    return CRExt.UnitSphericalPoint(x[1] / r, x[2] / r, x[3] / r)
end

"""
    node_cell_grids(CRExt, manifold, space)

One `CellBasedGrid` per cube face whose cells are the node cells of that face,
an `(ne Nq) × (ne Nq)` structured quadrilateral mesh.
"""
function node_cell_grids(CRExt, manifold, space)
    mesh = CC.Spaces.topology(space).mesh
    _, ws = CC.Quadratures.quadrature_points(Float64, CC.Spaces.quadrature_style(space))
    Nq = length(ws)
    ζ = collect([-1.0; cumsum(ws) .- 1])
    ζ[end] = 1.0
    M = mesh.ne * Nq
    return [
        CR.Trees.CellBasedGrid(
            manifold,
            [_node_cell_corner(CRExt, mesh, ζ, Nq, face, I₀, J₀) for I₀ in 0:M, J₀ in 0:M],
        ) for face in 1:6
    ]
end

"""
    node_cell_nodes_and_areas(CRExt, manifold, topology, Nq, grids)

For the node-cell meshes `grids` (from [`node_cell_grids`](@ref)), return the
flat SE node index of every node cell, in the global cell numbering used by
the intersection tree, and the area of every node cell indexed by SE node.
"""
function node_cell_nodes_and_areas(CRExt, manifold, topology, Nq, grids)
    M = topology.mesh.ne * Nq
    n_elem = CC.Topologies.nlocalelems(topology)
    elem_of = Dict(CRExt.element_face_local_indices(topology, e) => e for e in 1:n_elem)
    node_of_cell = zeros(Int, 6 * M^2)
    node_area = zeros(Float64, n_elem * Nq^2)
    for face in 1:6, c in 1:(M ^ 2)
        I, J = Tuple(CR.Trees.linear_to_cartesian_idx(grids[face], c))
        ie, i = (I - 1) ÷ Nq + 1, (I - 1) % Nq + 1
        je, j = (J - 1) ÷ Nq + 1, (J - 1) % Nq + 1
        n = (elem_of[(face, ie, je)] - 1) * Nq^2 + (j - 1) * Nq + i
        node_of_cell[(face - 1) * M ^ 2 + c] = n
        node_area[n] = CRExt.GO.area(manifold, CR.Trees.getcell(grids[face], c))
    end
    return node_of_cell, node_area
end

# Mask of polygons over wet (non-immersed) surface cells of `grid_oc`.
function _wet_polygon_mask(grid_oc, oc_of_poly)
    grid_oc isa OC.ImmersedBoundaryGrid || return trues(length(oc_of_poly))
    grid_cpu = OC.on_architecture(OC.CPU(), grid_oc)
    Nx, _, Nz = size(grid_cpu)
    return [
        !OC.ImmersedBoundaries.immersed_cell(mod1(c, Nx), (c - 1) ÷ Nx + 1, Nz, grid_cpu)
        for c in oc_of_poly
    ]
end

"""
    build_exchange_grid(boundary_space, grid_oc)

Construct an [`ExchangeGrid`](@ref) between a ClimaCore cubed-sphere
`boundary_space` and an Oceananigans `grid_oc`, on the CPU in `Float64`:
intersect the node cells with the FV cells (fold-aware on a `TripolarGrid`,
so fold shadow cells produce no polygons) and keep polygons over wet cells.
Move the result to the device with [`on_device`](@ref).
"""
function build_exchange_grid(boundary_space, grid_oc)
    CRExt = get_ConservativeRegriddingCCExt()
    @assert !isnothing(CRExt) "ConservativeRegriddingClimaCoreExt must be loaded"

    boundary_space_cpu = CC.Adapt.adapt(Array, boundary_space)
    grid_oc_underlying_cpu = OC.on_architecture(OC.CPU(), underlying_grid(grid_oc))

    FT = CC.Spaces.undertype(boundary_space_cpu)
    topology = CC.Spaces.topology(boundary_space_cpu)
    manifold = CR.Spherical(; radius = Float64(topology.mesh.domain.radius))
    Nq = CC.Quadratures.degrees_of_freedom(CC.Spaces.quadrature_style(boundary_space_cpu))
    n_elem = CC.Topologies.nlocalelems(topology)
    n_nodes = n_elem * Nq^2

    # 1. Node-cell × FV-cell intersection polygons.
    grids = node_cell_grids(CRExt, manifold, boundary_space_cpu)
    M = topology.mesh.ne * Nq
    dst_tree = CR.Trees.CubedSphereToplevelTree([
        CR.Trees.IndexOffsetQuadtreeCursor(grids[face], (face - 1) * M^2) for face in 1:6
    ])
    src_tree = CR.Trees.treeify(manifold, grid_oc_underlying_cpu)
    intersections = CR.intersection_areas(
        manifold,
        CR.False(),
        dst_tree,
        src_tree;
        intersection_operator = CR.IntersectionGridOperator(manifold),
    )
    cell_of_poly, oc_of_poly, polys = SparseArrays.findnz(intersections)
    n_oc = size(intersections, 2)

    # 2. Keep polygons over wet cells.
    keep = _wet_polygon_mask(grid_oc, oc_of_poly)
    cell_of_poly, oc_of_poly, polys = cell_of_poly[keep], oc_of_poly[keep], polys[keep]
    area = [CRExt.GO.area(manifold, poly) for poly in polys]
    n_poly = length(polys)

    # 3. Owning node and element of each polygon; node-cell areas.
    node_of_cell, node_area =
        node_cell_nodes_and_areas(CRExt, manifold, topology, Nq, grids)
    node_of_poly = node_of_cell[cell_of_poly]
    elem_of_poly = (node_of_poly .- 1) .÷ Nq^2 .+ 1

    # 4. Scatter: share of its node cell each polygon covers.
    share = area ./ node_area[node_of_poly]
    snode_ptr, spoly = _build_csr(node_of_poly, n_nodes)
    sweight = share[spoly]
    wet_share = zeros(Float64, n_nodes)
    for k in 1:n_poly
        wet_share[node_of_poly[k]] += share[k]
    end

    # 5. FV side.
    soc_ptr, soc_poly = _build_csr(Vector{Int}(oc_of_poly), n_oc)
    oc_wet_area = zeros(Float64, n_oc)
    for k in 1:n_poly
        oc_wet_area[oc_of_poly[k]] += area[k]
    end

    return ExchangeGrid{FT, Vector{Int32}, Vector{FT}}(
        n_poly,
        n_nodes,
        n_elem,
        n_oc,
        Int32.(node_of_poly),
        Int32.(elem_of_poly),
        Int32.(oc_of_poly),
        snode_ptr,
        spoly,
        FT.(sweight),
        FT.(node_area),
        FT.(wet_share),
        soc_ptr,
        soc_poly,
        FT.(area),
        FT.(oc_wet_area),
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

# dst[k] = src[owner[k]]
@kernel function _gather_by_owner_kernel!(dst, owner, src)
    k = @index(Global)
    @inbounds dst[k] = src[owner[k]]
end

function _gather_by_owner!(dst, owner, src)
    backend = KernelAbstractions.get_backend(dst)
    launch_kernel!(_gather_by_owner_kernel!, backend, length(dst), dst, owner, src)
    return dst
end

"""
    gather_nodes_to_polys!(poly_values, eg::ExchangeGrid, nodal_values)

Copy the value of the SE node whose cell contains each polygon from a flat
nodal vector.
"""
gather_nodes_to_polys!(poly_values, eg::ExchangeGrid, nodal_values) =
    _gather_by_owner!(poly_values, eg.node_of_poly, nodal_values)

"""
    gather_cells_to_polys!(poly_values, eg::ExchangeGrid, cell_values)

Copy the owning FV cell's value onto each polygon.
"""
gather_cells_to_polys!(poly_values, eg::ExchangeGrid, cell_values) =
    _gather_by_owner!(poly_values, eg.oc_of_poly, cell_values)

"""
    scatter_polys_to_nodes!(nodal_values, eg::ExchangeGrid, poly_values)

Share-weighted sum of per-polygon values into their node cells,
`F_n = Σ_{k ∈ n} sweight_k F_k`. Scattering ones gives the covered share of
each node cell. Follow with weighted DSS on the receiving field.
"""
scatter_polys_to_nodes!(nodal_values, eg::ExchangeGrid, poly_values) =
    _csr_matvec!(nodal_values, eg.snode_ptr, eg.spoly, eg.sweight, poly_values)

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
`F_c = Σ_{k ∈ c} area_k F_k / oc_wet_area_c`. Cells with zero wet area (dry,
fold shadows) are set to 0; run [`mirror_fold_partners!`](@ref) afterwards to
fill the shadow copies.
"""
function scatter_polys_to_cells!(cell_values, eg::ExchangeGrid, poly_values)
    backend = KernelAbstractions.get_backend(cell_values)
    launch_kernel!(
        _scatter_cells_kernel!,
        backend,
        eg.n_oc,
        cell_values,
        eg.soc_ptr,
        eg.soc_poly,
        eg.area,
        eg.oc_wet_area,
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
