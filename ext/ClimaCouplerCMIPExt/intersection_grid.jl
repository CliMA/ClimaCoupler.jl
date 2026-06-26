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
        Nx_oc, Ny_oc, _ = size(grid_with_mask_cpu)
        is_dry = Vector{Bool}(undef, Nx_oc * Ny_oc)
        @inbounds for j in 1:Ny_oc, i in 1:Nx_oc
            is_dry[(j - 1) * Nx_oc + i] =
                OC.ImmersedBoundaries.immersed_cell(i, j, 1, grid_with_mask_cpu)
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
    ArrayType;
    respect_immersed_mask::Bool = true,
)
    ig_cpu = extract_intersection_grid(
        boundary_space,
        grid_oc;
        respect_immersed_mask,
    )

    return IntersectionGrid(
        ArrayType(ig_cpu.cc_indices),
        ArrayType(ig_cpu.oc_indices),
        ArrayType(ig_cpu.areas),
        ArrayType(ig_cpu.cc_areas),
        ArrayType(ig_cpu.cc_total_areas),
        ArrayType(ig_cpu.oc_areas),
        ig_cpu.n_cc,
        ig_cpu.n_oc,
        ig_cpu.n_intersections,
        ig_cpu.n_nodes,
        ArrayType(ig_cpu.node_gather_polygon),
        ArrayType(ig_cpu.node_gather_node),
        ArrayType(ig_cpu.node_gather_weight),
    )
end

"""
    cc_element_areas(boundary_space) -> Vector

Spherical area [m²] of each spectral element on `boundary_space`, integrated
with the SEM quadrature weights (`local_geometry.WJ`).
"""
function cc_element_areas(boundary_space)
    space_cpu = CC.Adapt.adapt(Array, boundary_space)
    FT = CC.Spaces.undertype(space_cpu)
    lg = CC.Fields.local_geometry_field(CC.Fields.ones(space_cpu))
    wj = parent(CC.Fields.field_values(lg.WJ))
    n_elem = size(wj, ndims(wj))
    areas = zeros(FT, n_elem)
    @inbounds for e in 1:n_elem
        if ndims(wj) == 3
            areas[e] = sum(@view wj[:, :, e])
        elseif ndims(wj) == 4
            areas[e] = sum(@view wj[:, :, :, e])
        elseif ndims(wj) == 5
            areas[e] = sum(@view wj[:, :, :, :, e])
        else
            error("Unsupported WJ parent dimensions $(ndims(wj)) for element areas")
        end
    end
    return areas
end

"""
    wet_ocean_cc_fractions(ig::IntersectionGrid, cc_element_areas) -> Vector

Per-CC-element wet-ocean fraction of the boundary column: wet intersection
area divided by the full spectral-element area. Values lie in `[0, 1]`.
"""
function wet_ocean_cc_fractions(ig::IntersectionGrid, element_areas::AbstractVector)
    FT = eltype(element_areas)
    fractions = zeros(FT, ig.n_cc)
    @inbounds for i in 1:ig.n_cc
        if element_areas[i] > zero(FT)
            fractions[i] = ig.cc_areas[i] / element_areas[i]
        end
    end
    return clamp.(fractions, zero(FT), one(FT))
end

"""
    wet_ocean_fraction_field!(field, ig::IntersectionGrid, element_areas, boundary_space)

Broadcast per-element wet-ocean fractions onto all GLL nodes of `field`.
"""
function wet_ocean_fraction_field!(
    field,
    ig::IntersectionGrid,
    element_areas::AbstractVector,
    boundary_space,
)
    return _element_values_to_se_field!(
        field,
        wet_ocean_cc_fractions(ig, element_areas),
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
    @inbounds for k in 1:ig.n_intersections
        intersection_values[k] = cc_values[ig.cc_indices[k]]
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
    fill!(intersection_values, zero(eltype(intersection_values)))
    @inbounds for idx in eachindex(ig.node_gather_polygon)
        k = ig.node_gather_polygon[idx]
        n = ig.node_gather_node[idx]
        intersection_values[k] += ig.node_gather_weight[idx] * nodal_values[n]
    end
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
    @inbounds for k in 1:ig.n_intersections
        intersection_values[k] = oc_values[ig.oc_indices[k]]
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
    fill!(cc_values, zero(eltype(cc_values)))

    @inbounds for k in 1:ig.n_intersections
        i = ig.cc_indices[k]
        cc_values[i] += intersection_values[k] * ig.areas[k]
    end

    @inbounds for i in 1:ig.n_cc
        if ig.cc_areas[i] > 0
            cc_values[i] /= ig.cc_areas[i]
        end
    end
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
    fill!(oc_values, zero(eltype(oc_values)))

    @inbounds for k in 1:ig.n_intersections
        j = ig.oc_indices[k]
        oc_values[j] += intersection_values[k] * ig.areas[k]
    end

    @inbounds for j in 1:ig.n_oc
        if ig.oc_areas[j] > 0
            oc_values[j] /= ig.oc_areas[j]
        end
    end
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
    IntersectionFluxState(FT, n_intersections)

Allocate an `IntersectionFluxState` for `n_intersections` polygons.
"""
function IntersectionFluxState(FT::Type, n_intersections::Int)
    return IntersectionFluxState(
        zeros(FT, n_intersections),  # atmos_T
        zeros(FT, n_intersections),  # atmos_q_tot
        zeros(FT, n_intersections),  # atmos_q_liq
        zeros(FT, n_intersections),  # atmos_q_ice
        zeros(FT, n_intersections),  # atmos_ρ
        zeros(FT, n_intersections),  # atmos_u
        zeros(FT, n_intersections),  # atmos_v
        zeros(FT, n_intersections),  # atmos_h
        zeros(FT, n_intersections),  # surface_T
        zeros(FT, n_intersections),  # surface_q
        zeros(FT, n_intersections),  # surface_z0m
        zeros(FT, n_intersections),  # surface_z0b
        zeros(FT, n_intersections),  # surface_h
        zeros(FT, n_intersections),  # flux_sh
        zeros(FT, n_intersections),  # flux_lh
        zeros(FT, n_intersections),  # flux_τx
        zeros(FT, n_intersections),  # flux_τy
        zeros(FT, n_intersections),  # flux_evap
    )
end

"""
    IntersectionFluxState(FT, n_intersections, ArrayType)

Allocate an `IntersectionFluxState` for GPU computation using the specified `ArrayType`.
"""
function IntersectionFluxState(FT::Type, n_intersections::Int, ArrayType)
    return IntersectionFluxState(
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
        ArrayType(zeros(FT, n_intersections)),
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

    # Gather surface state to intersection grid
    gather_oc_to_intersection!(flux_state.surface_T, ig, oc_surface_state.T)
    gather_oc_to_intersection!(flux_state.surface_z0m, ig, oc_surface_state.z0m)
    gather_oc_to_intersection!(flux_state.surface_z0b, ig, oc_surface_state.z0b)
    gather_oc_to_intersection!(flux_state.surface_h, ig, oc_surface_state.h)

    # Compute surface humidity from surface temperature and density
    # Use atmosphere density extrapolated to surface as approximation
    FT = eltype(flux_state.atmos_T)
    @inbounds for k in 1:ig.n_intersections
        T_sfc = flux_state.surface_T[k]
        h_int = flux_state.atmos_h[k]
        h_sfc = flux_state.surface_h[k]
        Δz = h_int - h_sfc

        # Estimate surface density using barometric formula approximation
        ρ_sfc = SF.surface_density(
            surface_fluxes_params,
            flux_state.atmos_T[k],
            flux_state.atmos_ρ[k],
            T_sfc,
            Δz,
            flux_state.atmos_q_tot[k],
            FT(0),  # q_liq at surface
            FT(0),  # q_ice at surface
        )

        # Compute saturation specific humidity at surface
        flux_state.surface_q[k] = TD.q_vap_saturation(thermo_params, T_sfc, ρ_sfc, FT(0), FT(0))
    end

    use_ice_balance = !isnothing(ice_balance_at_int)
    use_ice_balance && (@assert !isnothing(σ) && !isnothing(T_melt))

    # Compute fluxes at each intersection polygon
    @inbounds for k in 1:ig.n_intersections
        if !isnothing(ice_active) && !ice_active[k]
            flux_state.flux_sh[k] = zero(FT)
            flux_state.flux_lh[k] = zero(FT)
            flux_state.flux_τx[k] = zero(FT)
            flux_state.flux_τy[k] = zero(FT)
            flux_state.flux_evap[k] = zero(FT)
            continue
        end

        # Build wind vector
        uv_int = StaticArrays.SVector(flux_state.atmos_u[k], flux_state.atmos_v[k])
        uv_sfc = StaticArrays.SVector(FT(0), FT(0))  # Surface velocity = 0

        # Build roughness parameters
        roughness_params = SF.ConstantRoughnessParams(
            flux_state.surface_z0m[k],
            flux_state.surface_z0b[k],
        )

        # Build flux configuration
        gustiness = FT(1)  # Default gustiness
        config = SF.SurfaceFluxConfig(roughness_params, SF.ConstantGustinessSpec(gustiness))

        # Compute heights
        h_sfc = flux_state.surface_h[k]
        h_int = flux_state.atmos_h[k]
        Φ_sfc = SFP.grav(surface_fluxes_params) * h_sfc
        Δz = h_int - h_sfc

        update_T_sfc_cb = if use_ice_balance
            ib = ice_balance_at_int
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
            flux_state.atmos_T[k],
            flux_state.atmos_q_tot[k],
            flux_state.atmos_q_liq[k],
            flux_state.atmos_q_ice[k],
            flux_state.atmos_ρ[k],
            flux_state.surface_T[k],
            flux_state.surface_q[k],
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
        flux_state.flux_sh[k] = outputs.shf
        flux_state.flux_lh[k] = outputs.lhf
        flux_state.flux_τx[k] = outputs.ρτxz
        flux_state.flux_τy[k] = outputs.ρτyz
        flux_state.flux_evap[k] = outputs.evaporation
        if use_ice_balance
            flux_state.surface_T[k] = outputs.T_sfc
        end
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
    gather_cc_nodal_to_intersection!(balance_at_int.SW_d, ig, balance_nodal.SW_d)
    gather_cc_nodal_to_intersection!(balance_at_int.LW_d, ig, balance_nodal.LW_d)
    gather_cc_nodal_to_intersection!(balance_at_int.δ, ig, balance_nodal.δ)
    gather_cc_nodal_to_intersection!(balance_at_int.T_i, ig, balance_nodal.T_i)
    gather_cc_nodal_to_intersection!(balance_at_int.ϵ, ig, balance_nodal.ϵ)
    gather_cc_nodal_to_intersection!(balance_at_int.α_albedo, ig, balance_nodal.α_albedo)
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
    return [sic_oc[ig.oc_indices[k]] > 0 for k in 1:ig.n_intersections]
end

"""
    _gather_cc_atmos_to_intersection!(flux_state, ig, cc_atmos_state)

Route atmosphere fields to the intersection grid via nodal polygon averages
when `ig.node_gather_*` is populated, otherwise fall back to per-element lookup.
"""
function _gather_cc_atmos_to_intersection!(flux_state, ig::IntersectionGrid, cc_atmos_state)
    if isempty(ig.node_gather_polygon)
        gather_cc_to_intersection!(flux_state.atmos_T, ig, cc_atmos_state.T)
        gather_cc_to_intersection!(flux_state.atmos_q_tot, ig, cc_atmos_state.q_tot)
        gather_cc_to_intersection!(flux_state.atmos_q_liq, ig, cc_atmos_state.q_liq)
        gather_cc_to_intersection!(flux_state.atmos_q_ice, ig, cc_atmos_state.q_ice)
        gather_cc_to_intersection!(flux_state.atmos_ρ, ig, cc_atmos_state.ρ)
        gather_cc_to_intersection!(flux_state.atmos_u, ig, cc_atmos_state.u)
        gather_cc_to_intersection!(flux_state.atmos_v, ig, cc_atmos_state.v)
        gather_cc_to_intersection!(flux_state.atmos_h, ig, cc_atmos_state.h)
    else
        gather_cc_nodal_to_intersection!(flux_state.atmos_T, ig, cc_atmos_state.T)
        gather_cc_nodal_to_intersection!(flux_state.atmos_q_tot, ig, cc_atmos_state.q_tot)
        gather_cc_nodal_to_intersection!(flux_state.atmos_q_liq, ig, cc_atmos_state.q_liq)
        gather_cc_nodal_to_intersection!(flux_state.atmos_q_ice, ig, cc_atmos_state.q_ice)
        gather_cc_nodal_to_intersection!(flux_state.atmos_ρ, ig, cc_atmos_state.ρ)
        gather_cc_nodal_to_intersection!(flux_state.atmos_u, ig, cc_atmos_state.u)
        gather_cc_nodal_to_intersection!(flux_state.atmos_v, ig, cc_atmos_state.v)
        gather_cc_nodal_to_intersection!(flux_state.atmos_h, ig, cc_atmos_state.h)
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
    scatter_flux_to_cc!(cc_fluxes.F_sh, ig, flux_state.flux_sh)
    scatter_flux_to_cc!(cc_fluxes.F_lh, ig, flux_state.flux_lh)
    scatter_flux_to_cc!(cc_fluxes.F_τx, ig, flux_state.flux_τx)
    scatter_flux_to_cc!(cc_fluxes.F_τy, ig, flux_state.flux_τy)
    scatter_flux_to_cc!(cc_fluxes.F_evap, ig, flux_state.flux_evap)
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
    scatter_flux_to_oc!(oc_fluxes.F_sh, ig, flux_state.flux_sh)
    scatter_flux_to_oc!(oc_fluxes.F_lh, ig, flux_state.flux_lh)
    scatter_flux_to_oc!(oc_fluxes.F_τx, ig, flux_state.flux_τx)
    scatter_flux_to_oc!(oc_fluxes.F_τy, ig, flux_state.flux_τy)
    scatter_flux_to_oc!(oc_fluxes.F_evap, ig, flux_state.flux_evap)
    return nothing
end
