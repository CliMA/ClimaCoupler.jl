"""
    IntersectionGridOperator

Custom [`ConservativeRegridding`](https://github.com/JuliaGeo/ConservativeRegridding.jl)
intersection operator that assembles a sparse matrix of intersection *polygons*
rather than scalar areas. Implements the return-type API from
[ConservativeRegridding.jl#121](https://github.com/JuliaGeo/ConservativeRegridding.jl/pull/121).
"""
struct IntersectionGridOperator{M <: CR.Manifold}
    manifold::M
end

IntersectionGridOperator(manifold::CR.Spherical) = IntersectionGridOperator{typeof(manifold)}(manifold)

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

CR.IntersectionReturnStyle(::IntersectionGridOperator) = CR.OutOfPlaceSingleResult()
CR.output_eltype(::IntersectionGridOperator) = _intersection_polygon_type()
CR.output_eltype(op::IntersectionGridOperator, src_tree, dst_tree) = CR.output_eltype(op)

function CR.should_store_result(::IntersectionGridOperator, result)
    result === nothing && return false
    CRExt = get_ConservativeRegriddingCCExt()
    return result isa CRExt.GI.Polygon
end

function (op::IntersectionGridOperator)(src_cell, dst_cell)
    CRExt = get_ConservativeRegriddingCCExt()
    GO = CRExt.GO
    GI = CRExt.GI
    intersection_poly = GO.intersection(
        GO.ConvexConvexSutherlandHodgman(op.manifold),
        src_cell,
        dst_cell;
        target = GI.PolygonTrait(),
    )
    if iszero(GO.area(op.manifold, intersection_poly))
        return nothing
    else
        return intersection_poly
    end
end

"""
    IntersectionGrid

Data structure representing the intersection grid between a ClimaCore
cubed-sphere spectral element mesh and an Oceananigans grid (typically the
production `TripolarGrid` wrapped in an `ImmersedBoundaryGrid`, but any OC
grid supported by ConservativeRegridding.jl works — including the legacy
`LatitudeLongitudeGrid` setup). The intersection grid consists of all
polygons formed by the overlap of CC elements and OC cells.

With a `TripolarGrid` covering the full sphere the intersection polygons
tile the entire CC topology, so no polar mask or polar truncation band is
needed; the polar-flux workaround used by the legacy lat-long setup has
been removed. This enables flux calculations directly on the intersection
polygons for better coastline representation, rather than computing fluxes
on one grid and remapping.

# Fields
- `cc_indices`: CC element index for each intersection polygon
- `oc_indices`: OC cell index for each intersection polygon  
- `areas`: Area of each intersection polygon [m²]
- `cc_areas`: Wet-ocean intersection area of each CC element [m²] (after the
  immersed mask filter when `respect_immersed_mask = true`)
- `cc_total_areas`: Total CC↔OC intersection area per element [m²] before the
  immersed mask filter (used for bathymetry-aligned surface fractions)
- `oc_areas`: Total area of each OC cell [m²]
- `n_cc`: Number of CC elements
- `n_oc`: Number of OC cells
- `n_intersections`: Number of intersection polygons
- `n_nodes`: Number of SE GLL nodes (`Nq² · Nh`) on the boundary space
- `node_gather_polygon`: COO polygon indices for [`gather_cc_nodal_to_intersection!`](@ref)
- `node_gather_node`: COO flat nodal indices (`1:n_nodes`)
- `node_gather_weight`: COO weights `Bᵢⱼ / areaₖ` integrating the SEM basis over polygon `k`

# Usage
The intersection grid is constructed from a boundary space + Oceananigans
grid pair via [`extract_intersection_grid`](@ref), which assembles the
element-level intersection matrix through `ConservativeRegridding.jl` using
[`IntersectionGridOperator`](@ref).
"""
struct IntersectionGrid{FT, IT, VFT <: AbstractVector{FT}, VIT <: AbstractVector{IT}}
    cc_indices::VIT
    oc_indices::VIT
    areas::VFT
    cc_areas::VFT
    cc_total_areas::VFT
    oc_areas::VFT
    n_cc::Int
    n_oc::Int
    n_intersections::Int
    n_nodes::Int
    node_gather_polygon::VIT
    node_gather_node::VIT
    node_gather_weight::VFT
end

function Base.show(io::IO, ig::IntersectionGrid{FT}) where {FT}
    print(io, "IntersectionGrid{$FT}: $(ig.n_intersections) polygons from $(ig.n_cc) CC elements × $(ig.n_oc) OC cells")
end

_oc_underlying(grid::OC.ImmersedBoundaryGrid) = grid.underlying_grid
_oc_underlying(grid) = grid

"""
    _se_n_nodes(boundary_space) -> Int

Number of flat GLL nodal degrees of freedom on a spectral-element boundary space.
"""
function _se_n_nodes(boundary_space)
    qs = CC.Spaces.quadrature_style(boundary_space)
    Nq = CC.Quadratures.degrees_of_freedom(qs)
    Nh = CC.Meshes.nelements(CC.Topologies.mesh(CC.Spaces.topology(boundary_space)))
    return Nq^2 * Nh
end

"""
    _normalize_intersection_polys(intersection_result)

Normalize `GeometryOps.intersection` output to an iterable of polygons.
"""
function _normalize_intersection_polys(intersection_result)
    intersection_result === nothing && return ()
    intersection_result isa AbstractVector && return intersection_result
    return (intersection_result,)
end

"""
    build_polygon_nodal_gather_weights(boundary_space_cpu, cc_indices,
                                       intersection_polys, areas)

Precompute COO weights that area-average a nodal SEM field onto each
intersection polygon:

    f̄ₖ = Σₙ wₖₙ fₙ,   wₖₙ = Bᵢⱼ(Ωₖ) / areaₖ

where `Bᵢⱼ` is the principled basis integral from
`ConservativeRegriddingClimaCoreExt.accumulate_principled_b`.

`intersection_polys` are the polygon values from a sparse intersection
matrix assembled with [`IntersectionGridOperator`](@ref).
"""
function build_polygon_nodal_gather_weights(
    boundary_space_cpu,
    cc_indices,
    intersection_polys,
    areas,
)
    CRExt = get_ConservativeRegriddingCCExt()
    @assert !isnothing(CRExt)

    FT = CC.Spaces.undertype(boundary_space_cpu)
    R = CC.Spaces.topology(boundary_space_cpu).mesh.domain.radius
    manifold = CR.Spherical(; radius = FT(R))

    qs = CC.Spaces.quadrature_style(boundary_space_cpu)
    Nq = CC.Quadratures.degrees_of_freedom(qs)
    triangle_quad_degree = 2 * (Nq - 1)

    node_gather_polygon = Int[]
    node_gather_node = Int[]
    node_gather_weight = FT[]

    @inbounds for k in eachindex(cc_indices)
        elem_idx = cc_indices[k]
        area_k = areas[k]
        area_k == 0 && continue

        B_total = zeros(Float64, Nq, Nq)
        for ipoly in _normalize_intersection_polys(intersection_polys[k])
            B_total .+= CRExt.accumulate_principled_b(
                manifold,
                boundary_space_cpu,
                elem_idx,
                ipoly;
                triangle_quad_degree,
            )
        end

        node_offset = (elem_idx - 1) * Nq^2
        inv_area = 1 / area_k
        for j in 1:Nq, i in 1:Nq
            Bᵢⱼ = B_total[i, j]
            Bᵢⱼ == 0 && continue
            push!(node_gather_polygon, k)
            push!(node_gather_node, node_offset + CC.Utilities.linear_ind((Nq, Nq), (i, j)))
            push!(node_gather_weight, FT(Bᵢⱼ * inv_area))
        end
    end

    return node_gather_polygon, node_gather_node, node_gather_weight
end

"""
    extract_intersection_grid(boundary_space, grid_oc;
                              respect_immersed_mask = true)

Build an `IntersectionGrid` directly from a ClimaCore boundary space and an
Oceananigans grid by computing the CC-element × OC-cell intersection matrix
via `ConservativeRegridding.intersection_areas` with an
[`IntersectionGridOperator`](@ref) that stores intersection polygons (see
[ConservativeRegridding.jl#121](https://github.com/JuliaGeo/ConservativeRegridding.jl/pull/121)).

This is the correct constructor for CC ↔ OC regridding because:

* The matrix it produces is indexed by SE *elements* on the CC side and
  by FV *cells* on the OC side — exactly what the gather/scatter and
  flux-on-intersection routines expect. The matrix exposed by the
  `CR.Regridder` built in `construct_remapper`, by contrast, is
  SE-*node*-indexed (`Nq² · Nh`), because `ConservativeRegriddingClimaCoreExt`
  bakes the L2 / principled SE operators into a node-level sparse matrix
  for direct use in `CR.regrid!`. Using that node-level matrix as if it
  were element-level would walk off the cubed-sphere tree.
* On a `TripolarGrid`, `CR.Trees.treeify(manifold, grid_oc_cpu)` routes
  through CR's Oceananigans extension and returns a fold-aware
  `PaddedTreeWrapper`, so the fold row is handled correctly without any
  caller-side bookkeeping.

The intersection grid is constructed against the underlying geometric grid -
the immersed boundary information is not used at this stage. Thus, the intersection
grid needs information about the dry/wet condition of the cell to correctly
compute ocean fluxes. Land fluxes are computed afterwards with the appropriate area-weighting.

Set `respect_immersed_mask = false` to keep every polygon, including
those that geometrically overlap immersed cells — useful for plotting
the raw CC↔OC geometry without coupling semantics layered on top.

# Arguments
- `boundary_space::CC.Spaces.AbstractSpace`: ClimaCore cubed-sphere boundary space.
- `grid_oc`: Oceananigans grid (`TripolarGrid`, `LatitudeLongitudeGrid`,
  `OrthogonalSphericalShellGrid`, or any of these wrapped in an
  `ImmersedBoundaryGrid`).

# Keyword arguments
- `respect_immersed_mask::Bool`: drop polygons whose parent OC cell is
  immersed. Defaults to `true`; no-op for non-`ImmersedBoundaryGrid`
  inputs.

# Returns
- `IntersectionGrid` containing the CC-element × OC-cell sparse intersection.
"""
function extract_intersection_grid(
    boundary_space,
    grid_oc;
    respect_immersed_mask::Bool = true,
)
    boundary_space_cpu = CC.Adapt.adapt(Array, boundary_space)
    grid_oc_cpu = OC.on_architecture(OC.CPU(), underlying_grid(grid_oc))

    FT = CC.Spaces.undertype(boundary_space_cpu)
    R = CC.Spaces.topology(boundary_space_cpu).mesh.domain.radius
    manifold = CR.Spherical(; radius = FT(R))

    dst_tree = CR.Trees.treeify(manifold, boundary_space_cpu)
    src_tree = CR.Trees.treeify(manifold, grid_oc_cpu)

    intersection_op = IntersectionGridOperator(manifold)
    intersections = CR.intersection_areas(
        manifold,
        CR.False(),
        dst_tree,
        src_tree;
        intersection_operator = intersection_op,
    )

    cc_indices, oc_indices, intersection_polys = SparseArrays.findnz(intersections)
    CRExt = get_ConservativeRegriddingCCExt()
    GO = CRExt.GO
    areas = [GO.area(manifold, poly) for poly in intersection_polys]

    n_cc, n_oc = size(intersections)
    cc_total_areas_vec = zeros(FT, n_cc)
    @inbounds for k in eachindex(areas)
        cc_total_areas_vec[cc_indices[k]] += areas[k]
    end

    # Optionally drop polygons whose parent OC cell is immersed (land /
    # dry under the bathymetry mask). We honor the live model layout — a
    # flat `(j - 1) * Nx + i` ordering matching `vec(OC.interior(...))`,
    # which is also what CR uses as the column index of the FV cell.
    if respect_immersed_mask && grid_oc isa OC.ImmersedBoundaryGrid
        grid_with_mask_cpu = OC.on_architecture(OC.CPU(), grid_oc)
        Nx_oc, Ny_oc, Nz_oc = size(grid_with_mask_cpu)
        is_dry = Vector{Bool}(undef, Nx_oc * Ny_oc)

        @inbounds for j in 1:Ny_oc, i in 1:Nx_oc
            is_dry[(j - 1) * Nx_oc + i] =
                OC.ImmersedBoundaries.immersed_cell(i, j, Nz_oc, grid_with_mask_cpu)
        end
        keep = .!is_dry[oc_indices]
        cc_indices = cc_indices[keep]
        oc_indices = oc_indices[keep]
        intersection_polys = intersection_polys[keep]
        areas = areas[keep]
    end

    # Recompute per-row / per-column totals from the (possibly filtered)
    # entries so that `ig.cc_areas[i]` and `ig.oc_areas[j]` agree with
    # the live `areas` vector. This keeps `scatter_to_cc!` / `scatter_to_oc!`
    # area-conserving on the polygons we actually retain.
    cc_areas_vec = zeros(FT, n_cc)
    oc_areas_vec = zeros(FT, n_oc)
    @inbounds for k in eachindex(areas)
        cc_areas_vec[cc_indices[k]] += areas[k]
        oc_areas_vec[oc_indices[k]] += areas[k]
    end

    n_intersections = length(areas)

    node_gather_polygon, node_gather_node, node_gather_weight =
        build_polygon_nodal_gather_weights(
            boundary_space_cpu,
            cc_indices,
            intersection_polys,
            areas,
        )
    n_nodes = _se_n_nodes(boundary_space_cpu)

    return IntersectionGrid(
        cc_indices,
        oc_indices,
        FT.(areas),
        cc_areas_vec,
        cc_total_areas_vec,
        oc_areas_vec,
        n_cc,
        n_oc,
        n_intersections,
        n_nodes,
        node_gather_polygon,
        node_gather_node,
        node_gather_weight,
    )
end

"""
    extract_intersection_grid(boundary_space, grid_oc, ArrayType;
                              respect_immersed_mask = true)

Build an `IntersectionGrid` (see the two-argument method) and convert its
backing arrays to `ArrayType` (e.g. `CuArray` for GPU computation).
"""
function extract_intersection_grid(
    boundary_space,
    grid_oc,
    arch;
    respect_immersed_mask::Bool = true,
)
    ig_cpu = extract_intersection_grid(
        boundary_space,
        grid_oc;
        respect_immersed_mask,
    )
    to_arch = x -> OC.on_architecture(arch, x)

    return IntersectionGrid(
        to_arch(ig_cpu.cc_indices),
        to_arch(ig_cpu.oc_indices),
        to_arch(ig_cpu.areas),
        to_arch(ig_cpu.cc_areas),
        to_arch(ig_cpu.cc_total_areas),
        to_arch(ig_cpu.oc_areas),
        ig_cpu.n_cc,
        ig_cpu.n_oc,
        ig_cpu.n_intersections,
        ig_cpu.n_nodes,
        to_arch(ig_cpu.node_gather_polygon),
        to_arch(ig_cpu.node_gather_node),
        to_arch(ig_cpu.node_gather_weight),
    )
end

"""
    intersection_array_type(arch) -> Type

Device array type for intersection-grid storage (`Array` on CPU, e.g. `CuArray` on GPU).
Matches the pattern used by [`ConservativeRegriddingClimaCoreExt`](https://juliageo.org/ConservativeRegridding.jl/stable/extensions/climacore/)
(`se_field_to_vec` / `regrid!` flat vectors on the compute device).
"""
intersection_array_type(::OC.CPU) = Array
function intersection_array_type(arch)
    return typeof(OC.on_architecture(arch, zeros(Float64, 1)))
end

"""
    extract_intersection_grid_on_arch(boundary_space, grid_oc, arch; kwargs...)

Build an [`IntersectionGrid`](@ref) and place backing arrays on the device for `arch`.
"""
function extract_intersection_grid_on_arch(boundary_space, grid_oc, arch; kwargs...)
    if arch isa OC.CPU
        return extract_intersection_grid(boundary_space, grid_oc; kwargs...)
    end
    return extract_intersection_grid(
        boundary_space,
        grid_oc,
        arch;
        kwargs...,
    )
end

zeros_on_arch(FT, ::OC.CPU, n) = zeros(FT, n)
zeros_on_arch(FT, arch, n) = OC.on_architecture(arch, zeros(FT, n))

"""Return a CPU vector for host-side intersection loops (no-op for `Array`)."""
host_intersection_vector(v::Array) = v
host_intersection_vector(v) = Array(v)

"""Copy a host vector back to its device counterpart when needed."""
function sync_intersection_vector!(device_v, host_v)
    device_v === host_v && return device_v
    copyto!(device_v, host_v)
    return device_v
end

on_gpu_intersection_array(v) = !(v isa Array)

function _host_intersection_namedtuple(nt::NamedTuple)
    return NamedTuple{keys(nt)}(map(host_intersection_vector, values(nt)))
end

function _sync_intersection_namedtuple!(device_nt, host_nt::NamedTuple)
    for name in keys(device_nt)
        sync_intersection_vector!(device_nt[name], host_nt[name])
    end
    return nothing
end

function _host_intersection_struct(s)
    names = fieldnames(typeof(s))
    return (; zip(names, map(n -> host_intersection_vector(getfield(s, n)), names))...)
end

function _sync_intersection_struct!(device_s, host_nt::NamedTuple)
    for name in keys(host_nt)
        sync_intersection_vector!(getfield(device_s, name), host_nt[name])
    end
    return nothing
end

_oc_surface_scratch(FT, arch, n) = (;
    T = zeros_on_arch(FT, arch, n),
    z0m = zeros_on_arch(FT, arch, n),
    z0b = zeros_on_arch(FT, arch, n),
    h = zeros_on_arch(FT, arch, n),
)

# `(cc_atmos_src, flux_state_dst)` pairs for atmosphere gather kernels.
const _ATMOS_FLUX_GATHER_PAIRS = (
    (:T, :atmos_T),
    (:q_tot, :atmos_q_tot),
    (:q_liq, :atmos_q_liq),
    (:q_ice, :atmos_q_ice),
    (:ρ, :atmos_ρ),
    (:u, :atmos_u),
    (:v, :atmos_v),
    (:h, :atmos_h),
)

const _SURFACE_FLUX_GATHER_PAIRS = (
    (:T, :surface_T),
    (:z0m, :surface_z0m),
    (:z0b, :surface_z0b),
    (:h, :surface_h),
)

const _CC_FLUX_AGGREGATE_PAIRS = (
    (:F_sh, :flux_sh),
    (:F_lh, :flux_lh),
    (:F_τx, :flux_τx),
    (:F_τy, :flux_τy),
    (:F_evap, :flux_evap),
)

"""
    allocate_intersection_flux_scratch(FT, arch, ig, n_oc_layout)

Allocate per-polygon and nodal scratch for [`compute_intersection_grid_fluxes!`](@ref).
"""
function allocate_intersection_flux_scratch(FT, arch, ig::IntersectionGrid, n_oc_layout)
    n_int = ig.n_intersections
    n_nodes = ig.n_nodes
    ArrayType = intersection_array_type(arch)

    return (;
        intersection_flux_state = IntersectionFluxState(FT, n_int, ArrayType),
        ice_intersection_flux_state = IntersectionFluxState(FT, n_int, ArrayType),
        cc_atmos_temp = (;
            T = zeros_on_arch(FT, arch, n_nodes),
            q_tot = zeros_on_arch(FT, arch, n_nodes),
            q_liq = zeros_on_arch(FT, arch, n_nodes),
            q_ice = zeros_on_arch(FT, arch, n_nodes),
            ρ = zeros_on_arch(FT, arch, n_nodes),
            u = zeros_on_arch(FT, arch, n_nodes),
            v = zeros_on_arch(FT, arch, n_nodes),
            h = zeros_on_arch(FT, arch, n_nodes),
        ),
        oc_surface_temp = _oc_surface_scratch(FT, arch, n_oc_layout),
        ice_cc_balance_nodal = (;
            SW_d = zeros_on_arch(FT, arch, n_nodes),
            LW_d = zeros_on_arch(FT, arch, n_nodes),
            δ = zeros_on_arch(FT, arch, n_nodes),
            T_i = zeros_on_arch(FT, arch, n_nodes),
            ϵ = zeros_on_arch(FT, arch, n_nodes),
            α_albedo = zeros_on_arch(FT, arch, n_nodes),
        ),
        ice_balance_at_int = (;
            κ = zeros_on_arch(FT, arch, n_int),
            SW_d = zeros_on_arch(FT, arch, n_int),
            LW_d = zeros_on_arch(FT, arch, n_int),
            δ = zeros_on_arch(FT, arch, n_int),
            T_i = zeros_on_arch(FT, arch, n_int),
            ϵ = zeros_on_arch(FT, arch, n_int),
            α_albedo = zeros_on_arch(FT, arch, n_int),
        ),
        ice_surface_temp = _oc_surface_scratch(FT, arch, n_oc_layout),
    )
end

"""
    wet_ocean_cc_fractions(ig::IntersectionGrid) -> AbstractVector

Per-CC-element wet-ocean fraction of the boundary column.

Computed as `ig.cc_areas[i] / ig.cc_total_areas[i]`:
- `cc_areas` — sum of polygon areas over **wet** OC cells (after the immersed mask filter).
- `cc_total_areas` — sum of all polygon areas for element `i` before the immersed mask
  filter (i.e. the full CC↔OC geometric overlap area of the element).

Both quantities share the same `GO.area` computation (on the same `Spherical` manifold),
so the ratio is independent of units or sphere-radius normalisation.
Values lie in `[0, 1]`.
"""
function wet_ocean_cc_fractions(ig::IntersectionGrid)
    FT = eltype(ig.areas)
    return clamp.(
        ifelse.(ig.cc_total_areas .> zero(FT), ig.cc_areas ./ ig.cc_total_areas, zero(FT)),
        zero(FT),
        one(FT),
    )
end

"""
    wet_ocean_fraction_field!(field, ig::IntersectionGrid, boundary_space)

Broadcast per-element wet-ocean fractions (see [`wet_ocean_cc_fractions`](@ref))
onto all GLL nodes of `field`.
"""
function wet_ocean_fraction_field!(field, ig::IntersectionGrid, boundary_space)
    return _element_values_to_se_field!(
        field,
        wet_ocean_cc_fractions(ig),
        boundary_space,
    )
end

"""
    gather_cc_to_intersection!(intersection_values, ig::IntersectionGrid, cc_values)

Gather values from CC elements to the intersection grid.

For each intersection polygon, look up the value from its parent CC element.

# Arguments
- `intersection_values`: Output vector of length `ig.n_intersections`
- `ig`: The intersection grid
- `cc_values`: Vector of values per CC element (length `ig.n_cc`)
"""
function gather_cc_to_intersection!(intersection_values, ig::IntersectionGrid, cc_values)
    if on_gpu_intersection_array(intersection_values)
        gather_cc_to_intersection_gpu!(
            intersection_values,
            ig,
            cc_values,
            get_backend(intersection_values),
        )
        return nothing
    end
    cc_host = host_intersection_vector(cc_values)
    cc_indices = host_intersection_vector(ig.cc_indices)
    @inbounds for k in 1:ig.n_intersections
        intersection_values[k] = cc_host[cc_indices[k]]
    end
    return nothing
end

"""
    gather_cc_nodal_to_intersection!(intersection_values, ig::IntersectionGrid, nodal_values)

Gather nodal SEM field values to the intersection grid using the precomputed
polygon-averaging weights (`node_gather_*`).

Each output entry is the area average of the SEM interpolant over the
corresponding intersection polygon:

    intersection_values[k] = Σₙ wₖₙ nodal_values[n]

# Arguments
- `intersection_values`: Output vector of length `ig.n_intersections`
- `ig`: `IntersectionGrid` with populated `node_gather_*` COO arrays
- `nodal_values`: Flat nodal vector of length `ig.n_nodes`
"""
function gather_cc_nodal_to_intersection!(
    intersection_values,
    ig::IntersectionGrid,
    nodal_values,
)
    iv_host = zeros(eltype(intersection_values), length(intersection_values))

    node_gather_polygon = host_intersection_vector(ig.node_gather_polygon)
    node_gather_node = host_intersection_vector(ig.node_gather_node)
    node_gather_weight = host_intersection_vector(ig.node_gather_weight)
    nodal_host = host_intersection_vector(nodal_values)
    @inbounds for idx in eachindex(node_gather_polygon)
        k = node_gather_polygon[idx]
        n = node_gather_node[idx]
        iv_host[k] += node_gather_weight[idx] * nodal_host[n]
    end
    sync_intersection_vector!(intersection_values, iv_host)
    return nothing
end

"""
    gather_oc_to_intersection!(intersection_values, ig::IntersectionGrid, oc_values)

Gather values from OC cells to the intersection grid.

For each intersection polygon, look up the value from its parent OC cell.

# Arguments
- `intersection_values`: Output vector of length `ig.n_intersections`
- `ig`: The intersection grid
- `oc_values`: Vector of values per OC cell (length `ig.n_oc`)
"""
function gather_oc_to_intersection!(intersection_values, ig::IntersectionGrid, oc_values)
    if on_gpu_intersection_array(intersection_values)
        gather_oc_to_intersection_gpu!(
            intersection_values,
            ig,
            oc_values,
            get_backend(intersection_values),
        )
        return nothing
    end
    oc_host = host_intersection_vector(oc_values)
    oc_indices = host_intersection_vector(ig.oc_indices)
    @inbounds for k in 1:ig.n_intersections
        intersection_values[k] = oc_host[oc_indices[k]]
    end
    return nothing
end

"""
    scatter_to_cc!(cc_values, ig::IntersectionGrid, intersection_values)

Scatter values from the intersection grid back to CC elements with area weighting.

For each CC element, compute the area-weighted average of values from all
its intersection polygons:

    cc_values[i] = Σ_k (intersection_values[k] × areas[k]) / cc_areas[i]

where the sum is over all polygons k belonging to CC element i.

# Arguments
- `cc_values`: Output vector of length `ig.n_cc`
- `ig`: The intersection grid
- `intersection_values`: Vector of values per intersection polygon
"""
function scatter_to_cc!(cc_values, ig::IntersectionGrid, intersection_values)
    cc_host = zeros(eltype(cc_values), ig.n_cc)
    iv = host_intersection_vector(intersection_values)
    cc_indices = host_intersection_vector(ig.cc_indices)
    areas = host_intersection_vector(ig.areas)
    cc_areas = host_intersection_vector(ig.cc_areas)

    fill!(cc_host, zero(eltype(cc_host)))
    @inbounds for k in 1:ig.n_intersections
        i = cc_indices[k]
        cc_host[i] += iv[k] * areas[k]
    end
    @inbounds for i in 1:ig.n_cc
        if cc_areas[i] > 0
            cc_host[i] /= cc_areas[i]
        end
    end
    sync_intersection_vector!(cc_values, cc_host)
    return nothing
end

"""
    scatter_to_oc!(oc_values, ig::IntersectionGrid, intersection_values)

Scatter values from the intersection grid back to OC cells with area weighting.

For each OC cell, compute the area-weighted average of values from all
its intersection polygons:

    oc_values[j] = Σ_k (intersection_values[k] × areas[k]) / oc_areas[j]

where the sum is over all polygons k belonging to OC cell j.

# Arguments
- `oc_values`: Output vector of length `ig.n_oc`
- `ig`: The intersection grid
- `intersection_values`: Vector of values per intersection polygon
"""
function scatter_to_oc!(oc_values, ig::IntersectionGrid, intersection_values)
    oc_host = zeros(eltype(oc_values), ig.n_oc)
    iv = host_intersection_vector(intersection_values)
    oc_indices = host_intersection_vector(ig.oc_indices)
    areas = host_intersection_vector(ig.areas)
    oc_areas = host_intersection_vector(ig.oc_areas)

    fill!(oc_host, zero(eltype(oc_host)))
    @inbounds for k in 1:ig.n_intersections
        j = oc_indices[k]
        oc_host[j] += iv[k] * areas[k]
    end
    @inbounds for j in 1:ig.n_oc
        if oc_areas[j] > 0
            oc_host[j] /= oc_areas[j]
        end
    end
    sync_intersection_vector!(oc_values, oc_host)
    return nothing
end

"""
    scatter_flux_to_cc!(cc_flux, ig::IntersectionGrid, intersection_flux)

Scatter flux densities from the intersection grid to CC elements.

For intensive quantities like flux densities (W/m²), we compute the area-weighted
mean, which preserves the total flux when integrated over the element.

This is equivalent to `scatter_to_cc!` but named explicitly for flux exchange.
"""
scatter_flux_to_cc!(cc_flux, ig::IntersectionGrid, intersection_flux) =
    scatter_to_cc!(cc_flux, ig, intersection_flux)

"""
    scatter_flux_to_oc!(oc_flux, ig::IntersectionGrid, intersection_flux)

Scatter flux densities from the intersection grid to OC cells.

For intensive quantities like flux densities (W/m²), we compute the area-weighted
mean, which preserves the total flux when integrated over the cell.

This is equivalent to `scatter_to_oc!` but named explicitly for flux exchange.
"""
scatter_flux_to_oc!(oc_flux, ig::IntersectionGrid, intersection_flux) =
    scatter_to_oc!(oc_flux, ig, intersection_flux)

# GPU kernel versions using KernelAbstractions
@kernel function gather_cc_to_intersection_kernel!(intersection_values, cc_indices, cc_values)
    k = @index(Global)
    @inbounds intersection_values[k] = cc_values[cc_indices[k]]
end

@kernel function gather_oc_to_intersection_kernel!(intersection_values, oc_indices, oc_values)
    k = @index(Global)
    @inbounds intersection_values[k] = oc_values[oc_indices[k]]
end

"""
    gather_cc_to_intersection_gpu!(intersection_values, ig::IntersectionGrid, cc_values, backend)

GPU-accelerated version of `gather_cc_to_intersection!`.
"""
function gather_cc_to_intersection_gpu!(intersection_values, ig::IntersectionGrid, cc_values, backend)
    kernel! = gather_cc_to_intersection_kernel!(backend)
    kernel!(intersection_values, ig.cc_indices, cc_values; ndrange=ig.n_intersections)
    return nothing
end

"""
    gather_oc_to_intersection_gpu!(intersection_values, ig::IntersectionGrid, oc_values, backend)

GPU-accelerated version of `gather_oc_to_intersection!`.
"""
function gather_oc_to_intersection_gpu!(intersection_values, ig::IntersectionGrid, oc_values, backend)
    kernel! = gather_oc_to_intersection_kernel!(backend)
    kernel!(intersection_values, ig.oc_indices, oc_values; ndrange=ig.n_intersections)
    return nothing
end

# Note: GPU scatter operations require atomic operations due to race conditions
# when multiple intersection polygons write to the same CC/OC element.
# For now, scatter is performed on CPU. A GPU version would need:
# - Atomic add operations for the accumulation
# - Separate normalization kernel

"""
    IntersectionFluxState{FT, VFT}

Temporary storage for flux calculations on the intersection grid.

# Fields
- `atmos_T`: Atmosphere temperature at each intersection [K]
- `atmos_q_tot`: Atmosphere total specific humidity at each intersection [kg/kg]
- `atmos_q_liq`: Atmosphere liquid specific humidity [kg/kg]
- `atmos_q_ice`: Atmosphere ice specific humidity [kg/kg]
- `atmos_ρ`: Atmosphere density [kg/m³]
- `atmos_u`: Atmosphere u-velocity [m/s]
- `atmos_v`: Atmosphere v-velocity [m/s]
- `atmos_h`: Atmosphere height [m]
- `surface_T`: Surface temperature [K]
- `surface_q`: Surface specific humidity [kg/kg]
- `surface_z0m`: Surface momentum roughness length [m]
- `surface_z0b`: Surface buoyancy roughness length [m]
- `surface_h`: Surface height [m]
- `flux_sh`: Sensible heat flux [W/m²]
- `flux_lh`: Latent heat flux [W/m²]
- `flux_τx`: Momentum flux x-component [N/m²]
- `flux_τy`: Momentum flux y-component [N/m²]
- `flux_evap`: Evaporation [kg/m²/s]
"""
struct IntersectionFluxState{FT, VFT <: AbstractVector{FT}}
    atmos_T::VFT
    atmos_q_tot::VFT
    atmos_q_liq::VFT
    atmos_q_ice::VFT
    atmos_ρ::VFT
    atmos_u::VFT
    atmos_v::VFT
    atmos_h::VFT
    surface_T::VFT
    surface_q::VFT
    surface_z0m::VFT
    surface_z0b::VFT
    surface_h::VFT
    flux_sh::VFT
    flux_lh::VFT
    flux_τx::VFT
    flux_τy::VFT
    flux_evap::VFT
end

"""
    IntersectionFluxState(FT, n_intersections; ArrayType = Array)

Allocate zeroed per-polygon scratch vectors for every field of
[`IntersectionFluxState`](@ref). `ArrayType` is `Array` on CPU or e.g. `CuArray`
on GPU (see [`intersection_array_type`](@ref)).
"""
function IntersectionFluxState(
    FT::Type,
    n_intersections::Int,
    ArrayType = Array,
)
    nfields = length(fieldnames(IntersectionFluxState))
    return IntersectionFluxState(
        ntuple(_ -> ArrayType(zeros(FT, n_intersections)), nfields)...,
    )
end

import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import StaticArrays

"""
    compute_surface_fluxes_on_intersection!(
        flux_state::IntersectionFluxState,
        ig::IntersectionGrid,
        cc_atmos_state,
        oc_surface_state,
        surface_fluxes_params,
        thermo_params,
    )

Compute turbulent surface fluxes on the intersection grid.

Atmosphere state is gathered from nodal SEM values using the
precomputed polygon-averaging weights in `ig.node_gather_*`. When those
weights are empty (e.g. a manually constructed degenerate `IntersectionGrid`),
the legacy per-element lookup via `cc_indices` is used instead.

# Arguments
- `flux_state`: `IntersectionFluxState` to store computed fluxes
- `ig`: `IntersectionGrid` defining the CC/OC intersection polygons
- `cc_atmos_state`: NamedTuple of atmosphere state vectors. With polygon
  weights populated, each field must be a flat nodal vector of length
  `ig.n_nodes`; otherwise per-CC-element vectors of length `ig.n_cc`:
  - `T`: Temperature [K]
  - `q_tot`: Total specific humidity [kg/kg]
  - `q_liq`: Liquid specific humidity [kg/kg]
  - `q_ice`: Ice specific humidity [kg/kg]
  - `ρ`: Density [kg/m³]
  - `u`, `v`: Wind components [m/s]
  - `h`: Height [m]
- `oc_surface_state`: NamedTuple of OC surface state vectors (per OC cell):
  - `T`: Surface temperature [K]
  - `z0m`: Momentum roughness length [m]
  - `z0b`: Buoyancy roughness length [m]
  - `h`: Surface height [m] (typically 0)
- `surface_fluxes_params`: SurfaceFluxes.Parameters.SurfaceFluxesParameters
- `thermo_params`: Thermodynamics parameters for humidity calculations

# Keyword arguments (sea ice)
- `ice_balance_at_int`: optional `NamedTuple` of per-intersection vectors
  (`κ`, `δ`, `T_i`, `ϵ`, `α_albedo`, `SW_d`, `LW_d`) for the conductive
  skin-temperature balance. When provided, [`update_T_sfc`](@ref) is invoked
  at each polygon instead of using a prescribed surface temperature.
- `σ`, `T_melt`: scalars required when `ice_balance_at_int` is provided.
- `ice_active`: optional `BitVector` / `Vector{Bool}` of length
  `ig.n_intersections`; inactive polygons receive zero flux.
"""
function compute_surface_fluxes_on_intersection!(
    flux_state::IntersectionFluxState,
    ig::IntersectionGrid,
    cc_atmos_state::NamedTuple,
    oc_surface_state::NamedTuple,
    surface_fluxes_params,
    thermo_params;
    ice_balance_at_int = nothing,
    σ = nothing,
    T_melt = nothing,
    ice_active = nothing,
)
    _gather_cc_atmos_to_intersection!(flux_state, ig, cc_atmos_state)

    for (src, dst) in _SURFACE_FLUX_GATHER_PAIRS
        gather_oc_to_intersection!(getfield(flux_state, dst), ig, getfield(oc_surface_state, src))
    end

    flux = _host_intersection_struct(flux_state)
    ice_balance = isnothing(ice_balance_at_int) ? nothing : _host_intersection_namedtuple(ice_balance_at_int)

    # Compute surface humidity from surface temperature and density
    # Use atmosphere density extrapolated to surface as approximation
    FT = eltype(flux.atmos_T)
    @inbounds for k in 1:ig.n_intersections
        T_sfc = flux.surface_T[k]
        h_int = flux.atmos_h[k]
        h_sfc = flux.surface_h[k]
        Δz = h_int - h_sfc

        # Estimate surface density using barometric formula approximation
        ρ_sfc = SF.surface_density(
            surface_fluxes_params,
            flux.atmos_T[k],
            flux.atmos_ρ[k],
            T_sfc,
            Δz,
            flux.atmos_q_tot[k],
            FT(0),  # q_liq at surface
            FT(0),  # q_ice at surface
        )

        # Compute saturation specific humidity at surface
        flux.surface_q[k] = TD.q_vap_saturation(thermo_params, T_sfc, ρ_sfc, FT(0), FT(0))
    end

    use_ice_balance = !isnothing(ice_balance_at_int)
    use_ice_balance && (@assert !isnothing(σ) && !isnothing(T_melt))

    # Compute fluxes at each intersection polygon
    @inbounds for k in 1:ig.n_intersections
        if !isnothing(ice_active) && !ice_active[k]
            flux.flux_sh[k] = zero(FT)
            flux.flux_lh[k] = zero(FT)
            flux.flux_τx[k] = zero(FT)
            flux.flux_τy[k] = zero(FT)
            flux.flux_evap[k] = zero(FT)
            continue
        end

        # Build wind vector
        uv_int = StaticArrays.SVector(flux.atmos_u[k], flux.atmos_v[k])
        uv_sfc = StaticArrays.SVector(FT(0), FT(0))  # Surface velocity = 0

        # Build roughness parameters
        roughness_params = SF.ConstantRoughnessParams(
            flux.surface_z0m[k],
            flux.surface_z0b[k],
        )

        # Build flux configuration
        gustiness = FT(1)  # Default gustiness
        config = SF.SurfaceFluxConfig(roughness_params, SF.ConstantGustinessSpec(gustiness))

        # Compute heights
        h_sfc = flux.surface_h[k]
        h_int = flux.atmos_h[k]
        Φ_sfc = SFP.grav(surface_fluxes_params) * h_sfc
        Δz = h_int - h_sfc

        update_T_sfc_cb = if use_ice_balance
            ib = ice_balance
            update_T_sfc(
                ib.κ[k],
                ib.δ[k],
                ib.T_i[k],
                σ,
                ib.ϵ[k],
                ib.SW_d[k],
                ib.LW_d[k],
                ib.α_albedo[k],
                T_melt,
            )
        else
            nothing
        end

        # Call SurfaceFluxes
        outputs = SF.surface_fluxes(
            surface_fluxes_params,
            flux.atmos_T[k],
            flux.atmos_q_tot[k],
            flux.atmos_q_liq[k],
            flux.atmos_q_ice[k],
            flux.atmos_ρ[k],
            flux.surface_T[k],
            flux.surface_q[k],
            Φ_sfc,
            Δz,
            FT(0),  # displacement height d
            uv_int,
            uv_sfc,
            roughness_params,
            config,
            SF.PointValueScheme(),
            nothing,  # solver_opts
            nothing,  # flux_specs
            update_T_sfc_cb,
            nothing,  # update_q_vap_sfc
        )

        # Store flux outputs
        flux.flux_sh[k] = outputs.shf
        flux.flux_lh[k] = outputs.lhf
        flux.flux_τx[k] = outputs.ρτxz
        flux.flux_τy[k] = outputs.ρτyz
        flux.flux_evap[k] = outputs.evaporation
        if use_ice_balance
            flux.surface_T[k] = outputs.T_sfc
        end
    end

    _sync_intersection_struct!(flux_state, flux)
    if use_ice_balance
        _sync_intersection_namedtuple!(ice_balance_at_int, ice_balance)
    end

    return nothing
end

"""
    gather_ice_balance_to_intersection!(balance_at_int, ig, balance_nodal, κ)

Polygon-average ice skin-balance inputs from flat nodal SEM vectors onto the
intersection grid. `balance_nodal` holds `SW_d`, `LW_d`, `δ`, `T_i`, `ϵ`,
`α_albedo` nodal vectors; `κ` may be a scalar or nodal vector.
"""
function gather_ice_balance_to_intersection!(
    balance_at_int::NamedTuple,
    ig::IntersectionGrid,
    balance_nodal::NamedTuple,
    κ,
)
    for f in (:SW_d, :LW_d, :δ, :T_i, :ϵ, :α_albedo)
        gather_cc_nodal_to_intersection!(
            getfield(balance_at_int, f),
            ig,
            getfield(balance_nodal, f),
        )
    end
    if κ isa Number
        fill!(balance_at_int.κ, κ)
    else
        gather_cc_nodal_to_intersection!(balance_at_int.κ, ig, κ)
    end
    return nothing
end

"""
    ice_intersection_active_mask(ig, sic_oc) -> Vector{Bool}

Return a per-polygon mask that is `true` when the parent OC cell carries ice.
"""
function ice_intersection_active_mask(ig::IntersectionGrid, sic_oc::AbstractVector)
    sic_host = host_intersection_vector(sic_oc)
    oc_indices = host_intersection_vector(ig.oc_indices)
    return [sic_host[oc_indices[k]] > 0 for k in 1:ig.n_intersections]
end

"""
    _gather_cc_atmos_to_intersection!(flux_state, ig, cc_atmos_state)

Route atmosphere fields to the intersection grid via nodal polygon averages
when `ig.node_gather_*` is populated, otherwise fall back to per-element lookup.
"""
function _gather_cc_atmos_to_intersection!(flux_state, ig::IntersectionGrid, cc_atmos_state)
    gather! =
        isempty(ig.node_gather_polygon) ?
        gather_cc_to_intersection! :
        gather_cc_nodal_to_intersection!
    for (src, dst) in _ATMOS_FLUX_GATHER_PAIRS
        gather!(getfield(flux_state, dst), ig, getfield(cc_atmos_state, src))
    end
    return nothing
end

"""
    intersection_fluxes_to_boundary_fields(boundary_space, intersection_grid, intersection_flux_state)

Area-average polygon flux densities onto CC elements and broadcast to boundary-space
Fields for `update_flux_fields!`.
"""
function intersection_fluxes_to_boundary_fields(
    boundary_space,
    intersection_grid,
    intersection_flux_state,
)
    FT = CC.Spaces.undertype(boundary_space)
    n_cc = intersection_grid.n_cc

    cc_F_sh = zeros(FT, n_cc)
    cc_F_lh = zeros(FT, n_cc)
    cc_F_τx = zeros(FT, n_cc)
    cc_F_τy = zeros(FT, n_cc)
    cc_F_evap = zeros(FT, n_cc)
    aggregate_fluxes_to_cc!(
        (; F_sh = cc_F_sh, F_lh = cc_F_lh, F_τx = cc_F_τx, F_τy = cc_F_τy, F_evap = cc_F_evap),
        intersection_flux_state,
        intersection_grid,
    )

    F_sh = CC.Fields.zeros(boundary_space)
    F_lh = CC.Fields.zeros(boundary_space)
    F_turb_ρτxz = CC.Fields.zeros(boundary_space)
    F_turb_ρτyz = CC.Fields.zeros(boundary_space)
    F_turb_moisture = CC.Fields.zeros(boundary_space)

    _element_values_to_se_field!(F_sh, cc_F_sh, boundary_space)
    _element_values_to_se_field!(F_lh, cc_F_lh, boundary_space)
    _element_values_to_se_field!(F_turb_ρτxz, cc_F_τx, boundary_space)
    _element_values_to_se_field!(F_turb_ρτyz, cc_F_τy, boundary_space)
    _element_values_to_se_field!(F_turb_moisture, cc_F_evap, boundary_space)

    return (; F_turb_ρτxz, F_turb_ρτyz, F_lh, F_sh, F_turb_moisture)
end

"""
    _element_values_to_se_field!(field, element_values, boundary_space)

Broadcast per-SE-element scalars to all GLL nodes of each element.
"""
function _element_values_to_se_field!(field, element_values::AbstractVector, boundary_space)
    # Per-element vectors are small; host copy avoids GPU scalar indexing here.
    element_values = Array(element_values)
    p = parent(CC.Fields.field_values(field))
    if ndims(p) == 4
        @inbounds for e in eachindex(element_values)
            p[:, :, 1, e] .= element_values[e]
        end
    elseif ndims(p) == 5
        @inbounds for e in eachindex(element_values)
            p[:, :, 1, 1, e] .= element_values[e]
        end
    else
        error("Unsupported field parent dimensions $(ndims(p)) for element broadcast")
    end
    return field
end

"""
    aggregate_fluxes_to_cc!(cc_fluxes, flux_state, ig)

Aggregate fluxes from the intersection grid to CC elements.

# Arguments
- `cc_fluxes`: NamedTuple with vectors for each flux type (length `ig.n_cc`)
- `flux_state`: `IntersectionFluxState` containing computed fluxes
- `ig`: `IntersectionGrid`
"""
function aggregate_fluxes_to_cc!(cc_fluxes::NamedTuple, flux_state::IntersectionFluxState, ig::IntersectionGrid)
    for (dst, src) in _CC_FLUX_AGGREGATE_PAIRS
        scatter_flux_to_cc!(getfield(cc_fluxes, dst), ig, getfield(flux_state, src))
    end
    return nothing
end

"""
    aggregate_fluxes_to_oc!(oc_fluxes, flux_state, ig)

Aggregate fluxes from the intersection grid to OC cells.

# Arguments
- `oc_fluxes`: NamedTuple with vectors for each flux type (length `ig.n_oc`)
- `flux_state`: `IntersectionFluxState` containing computed fluxes
- `ig`: `IntersectionGrid`
"""
function aggregate_fluxes_to_oc!(oc_fluxes::NamedTuple, flux_state::IntersectionFluxState, ig::IntersectionGrid)
    for (dst, src) in _CC_FLUX_AGGREGATE_PAIRS
        scatter_flux_to_oc!(getfield(oc_fluxes, dst), ig, getfield(flux_state, src))
    end
    return nothing
end
