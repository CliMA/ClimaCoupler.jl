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
- `cc_areas`: Total area of each CC element [m²]
- `oc_areas`: Total area of each OC cell [m²]
- `n_cc`: Number of CC elements
- `n_oc`: Number of OC cells
- `n_intersections`: Number of intersection polygons

# Usage
The intersection grid is constructed from a boundary space + Oceananigans
grid pair via [`extract_intersection_grid`](@ref), which builds the
element-level intersection matrix through `ConservativeRegridding`'s
tree-level dual-DFS (fold-aware on `TripolarGrid` via CR's Oceananigans
extension). Values can then be gathered from CC/OC grids to the
intersection grid, flux calculations performed, and results scattered
back to component grids with proper area weighting.
"""
struct IntersectionGrid{FT, IT, VFT <: AbstractVector{FT}, VIT <: AbstractVector{IT}}
    cc_indices::VIT
    oc_indices::VIT
    areas::VFT
    cc_areas::VFT
    oc_areas::VFT
    n_cc::Int
    n_oc::Int
    n_intersections::Int
end

function Base.show(io::IO, ig::IntersectionGrid{FT}) where {FT}
    print(io, "IntersectionGrid{$FT}: $(ig.n_intersections) polygons from $(ig.n_cc) CC elements × $(ig.n_oc) OC cells")
end

"""
    extract_intersection_grid(boundary_space, grid_oc;
                              respect_immersed_mask = true)

Build an `IntersectionGrid` directly from a ClimaCore boundary space and an
Oceananigans grid by computing the CC-element × OC-cell intersection matrix
via the tree-level `ConservativeRegridding.intersection_areas` call.

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

    intersections = CR.intersection_areas(manifold, CR.False(), dst_tree, src_tree)

    cc_indices, oc_indices, areas = SparseArrays.findnz(intersections)

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
        areas = areas[keep]
    end

    # Recompute per-row / per-column totals from the (possibly filtered)
    # entries so that `ig.cc_areas[i]` and `ig.oc_areas[j]` agree with
    # the live `areas` vector. This keeps `scatter_to_cc!` / `scatter_to_oc!`
    # area-conserving on the polygons we actually retain.
    n_cc, n_oc = size(intersections)
    cc_areas_vec = zeros(FT, n_cc)
    oc_areas_vec = zeros(FT, n_oc)
    @inbounds for k in eachindex(areas)
        cc_areas_vec[cc_indices[k]] += areas[k]
        oc_areas_vec[oc_indices[k]] += areas[k]
    end

    n_intersections = length(areas)

    return IntersectionGrid(
        cc_indices,
        oc_indices,
        FT.(areas),
        cc_areas_vec,
        oc_areas_vec,
        n_cc,
        n_oc,
        n_intersections,
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
        ArrayType(ig_cpu.oc_areas),
        ig_cpu.n_cc,
        ig_cpu.n_oc,
        ig_cpu.n_intersections,
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

This function gathers atmosphere state from CC elements and surface state from OC cells
to each intersection polygon, computes fluxes using SurfaceFluxes.jl, and stores the
results in `flux_state`.

# Arguments
- `flux_state`: `IntersectionFluxState` to store computed fluxes
- `ig`: `IntersectionGrid` defining the CC/OC intersection polygons
- `cc_atmos_state`: NamedTuple of CC atmosphere state vectors (per CC element):
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
"""
function compute_surface_fluxes_on_intersection!(
    flux_state::IntersectionFluxState,
    ig::IntersectionGrid,
    cc_atmos_state::NamedTuple,
    oc_surface_state::NamedTuple,
    surface_fluxes_params,
    thermo_params,
)
    # Gather atmosphere state to intersection grid
    gather_cc_to_intersection!(flux_state.atmos_T, ig, cc_atmos_state.T)
    gather_cc_to_intersection!(flux_state.atmos_q_tot, ig, cc_atmos_state.q_tot)
    gather_cc_to_intersection!(flux_state.atmos_q_liq, ig, cc_atmos_state.q_liq)
    gather_cc_to_intersection!(flux_state.atmos_q_ice, ig, cc_atmos_state.q_ice)
    gather_cc_to_intersection!(flux_state.atmos_ρ, ig, cc_atmos_state.ρ)
    gather_cc_to_intersection!(flux_state.atmos_u, ig, cc_atmos_state.u)
    gather_cc_to_intersection!(flux_state.atmos_v, ig, cc_atmos_state.v)
    gather_cc_to_intersection!(flux_state.atmos_h, ig, cc_atmos_state.h)

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

    # Compute fluxes at each intersection polygon
    @inbounds for k in 1:ig.n_intersections
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
            nothing,  # update_T_sfc_cb
            nothing,  # update_q_vap_sfc
        )

        # Store flux outputs
        flux_state.flux_sh[k] = outputs.shf
        flux_state.flux_lh[k] = outputs.lhf
        flux_state.flux_τx[k] = outputs.ρτxz
        flux_state.flux_τy[k] = outputs.ρτyz
        flux_state.flux_evap[k] = outputs.evaporation
    end

    return nothing
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
