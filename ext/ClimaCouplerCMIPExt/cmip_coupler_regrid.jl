# Conservative (finite-volume) regridding between the ClimaCore cubed-sphere exchange grid
# and Oceananigans `LatitudeLongitudeGrid`, via ConservativeRegridding.jl (see ClimaCoupler PR #1724).
#
# GPU: `LinearAlgebra.mul!` on `CUDA.CUSPARSE.CuSparseMatrixCSC` falls back to scalar indexing; we call
# `CUDA.CUSPARSE.mv!` instead when the intersection matrix lives on CUDA.CUSPARSE.
#
# Weights for remapping MUST come from `_conservative_weights_alias` / `_cmip_conservative_weights` below.
# Do not call the optional `ConservativeRegridMath` submodule on `ClimaCoupler` from this file — it is not
# guaranteed to exist in every precompiled or `] dev` checkout (UndefVarError at `construct_...`).

const _CUDA_PKGID = Base.PkgId(Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA")

function _cuda_cusparse_module()
    CUDA = get(Base.loaded_modules, _CUDA_PKGID, nothing)
    isnothing(CUDA) && return nothing
    return Base.getproperty(CUDA, :CUSPARSE)
end

"""
True if `A` is a CUDA CUSPARSE sparse matrix and we should use `CUSPARSE.mv!` instead of
`LinearAlgebra.mul!` (which scalar-indexes `CuSparseMatrixCSC` on the host).

Uses `A isa CUSPARSE.AbstractCuSparseMatrix` so this still works when `parentmodule(typeof(A))`
does not compare equal to the `CUSPARSE` module object (can happen across Julia/CUDA load orders).
"""
function _use_cusparse_spmv(A)
    CUS = _cuda_cusparse_module()
    isnothing(CUS) && return false
    if isdefined(CUS, :AbstractCuSparseMatrix)
        return A isa CUS.AbstractCuSparseMatrix
    end
    return parentmodule(typeof(A)) === CUS
end

"""`y = op(A) * x` with `op` either identity (`'N'`) or transpose (`'T'`); GPU uses CUSPARSE SpMV."""
function _intersection_spmv!(y::AbstractVector, A, x::AbstractVector, trans::Char)
    if _use_cusparse_spmv(A)
        CUS = _cuda_cusparse_module()::Module
        T = eltype(y)
        CUS.mv!(trans, one(T), A, x, zero(T), y, 'O')
    elseif trans == 'N'
        LinearAlgebra.mul!(y, A, x)
    else
        LinearAlgebra.mul!(y, LinearAlgebra.transpose(A), x)
    end
    return y
end

"""Alias regridder fields (equivalent to `precompute_from_regridder(...; copy_arrays = false)` in tests)."""
function _conservative_weights_alias(remapper::CR.Regridder)
    return (;
        A = remapper.intersections,
        dst_areas = remapper.dst_areas,
        src_areas = remapper.src_areas,
    )
end

"""Public name used from `construct_conservative_ocean_coupler_remapping` (single place to change if needed)."""
_cmip_conservative_weights(remapper::CR.Regridder) = _conservative_weights_alias(remapper)

function cmip_conservative_regridding_cc_ext()
    ext = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)
    isnothing(ext) && error(
        "ConservativeRegriddingClimaCoreExt is not available. Load `ConservativeRegridding` and ensure the ClimaCore extension is built.",
    )
    return ext
end

"""
    construct_conservative_ocean_coupler_remapping(grid, boundary_space)

Build area-conservative regridding operators between `boundary_space` (spectral-element
coupler grid) and the ocean model's horizontal `LatitudeLongitudeGrid` (underlying grid).

The returned named tuple includes:

  - `remapper_oc_to_cc`: `ConservativeRegridding.Regridder` on the ocean architecture (GPU/CPU),
    holding dense temporaries sized for sparse matvec;
  - `conservative_weights`: `(; A, dst_areas, src_areas)` **aliasing** the regridder's intersection
    matrix and area vectors on the same device (no duplicate sparse matrix on GPU).
"""
function construct_conservative_ocean_coupler_remapping(grid, boundary_space)
    grid_oc_underlying_cpu = OC.on_architecture(OC.CPU(), grid.underlying_grid)
    boundary_space_cpu = Adapt.adapt_structure(Array, boundary_space)
    # ConservativeRegridding requires identical manifold types for src/dst.
    # Force a shared spherical manifold to avoid Float32 vs Float64 radius mismatch.
    radius = try
        Float64(CC.Spaces.topology(boundary_space).mesh.domain.radius)
    catch
        6.371e6
    end
    remapper_oc_to_cc = CR.Regridder(
        CR.Spherical(; radius),
        boundary_space_cpu,
        grid_oc_underlying_cpu;
        normalize = false,
        threaded = false,
    )
    remapper_oc_to_cc = OC.on_architecture(OC.architecture(grid), remapper_oc_to_cc)
    conservative_weights = _cmip_conservative_weights(remapper_oc_to_cc)

    FT = CC.Spaces.undertype(boundary_space)
    field_ones_cc = CC.Fields.ones(boundary_space)
    ArrayType = ClimaComms.array_type(boundary_space)
    topo = CC.Spaces.topology(boundary_space)
    n_elem = CC.Topologies.nlocalelems(topo)
    value_per_element_cc = ArrayType(zeros(FT, n_elem))

    scratch_field_oc1 = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    scratch_field_oc2 = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    # Extra ocean surface scratch fields (e.g. `ocean_seaice_fluxes!` temporaries)
    scratch_cc1 = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    scratch_cc2 = OC.Field{OC.Center, OC.Center, Nothing}(grid)
    zsurf = size(scratch_field_oc1, 3)
    scratch_arr3 = similar(OC.interior(scratch_field_oc1, :, :, zsurf))
    temp_uv_vec = CC.Fields.Field(CC.Geometry.UVVector{FT}, boundary_space)

    polar_exclusion_flux_mask_centers =
        ocean_flux_highlat_mask(grid; location = (OC.Center(), OC.Center(), OC.Center()))
    polar_exclusion_flux_mask_u =
        ocean_flux_highlat_mask(grid; location = (OC.Face(), OC.Center(), OC.Center()))
    polar_exclusion_flux_mask_v =
        ocean_flux_highlat_mask(grid; location = (OC.Center(), OC.Face(), OC.Center()))

    return (;
        regridding = :conservative,
        remapper_oc_to_cc,
        conservative_weights,
        field_ones_cc,
        value_per_element_cc,
        scratch_field_oc1,
        scratch_field_oc2,
        scratch_cc1,
        scratch_cc2,
        scratch_arr3,
        temp_uv_vec,
        polar_exclusion_flux_mask_centers,
        polar_exclusion_flux_mask_u,
        polar_exclusion_flux_mask_v,
    )
end

function _dense_vec1(x)
    return x isa DenseArray{<:Any, 1}
end

"""
Ocean (source) → coupler per-element values: `dst_cc .= (A * src_oc) ./ a_cc`.

On GPU, CUSPARSE requires dense `CuVector` operands; we always route through `r.src_temp` /
`r.dst_temp` when using a CUDA sparse intersection matrix (also covers `ReshapedArray` / `SubArray` inputs).
"""
function _conservative_regrid_cc_from_oc!(
    dst_cc::AbstractVector,
    A,
    src_oc::AbstractVector,
    dst_areas,
    r::CR.Regridder,
)
    if _use_cusparse_spmv(A)
        r.src_temp .= src_oc
        _intersection_spmv!(r.dst_temp, A, r.src_temp, 'N')
        r.dst_temp ./= dst_areas
        dst_cc .= r.dst_temp
    elseif _dense_vec1(dst_cc) && _dense_vec1(src_oc)
        LinearAlgebra.mul!(dst_cc, A, src_oc)
        dst_cc ./= dst_areas
    else
        r.src_temp .= src_oc
        _intersection_spmv!(r.dst_temp, A, r.src_temp, 'N')
        r.dst_temp ./= dst_areas
        dst_cc .= r.dst_temp
    end
    return nothing
end

"""
Coupler per-element values → ocean cell means: `dst_oc .= (A' * src_cc) ./ a_oc`.
See [`_conservative_regrid_cc_from_oc!`](@ref) for GPU temporaries.
"""
function _conservative_regrid_oc_from_cc!(
    dst_oc::AbstractVector,
    A,
    src_cc::AbstractVector,
    src_areas,
    r::CR.Regridder,
)
    if _use_cusparse_spmv(A)
        r.dst_temp .= src_cc
        _intersection_spmv!(r.src_temp, A, r.dst_temp, 'T')
        r.src_temp ./= src_areas
        dst_oc .= r.src_temp
    elseif _dense_vec1(dst_oc) && _dense_vec1(src_cc)
        LinearAlgebra.mul!(dst_oc, LinearAlgebra.transpose(A), src_cc)
        dst_oc ./= src_areas
    else
        r.dst_temp .= src_cc
        _intersection_spmv!(r.src_temp, A, r.dst_temp, 'T')
        r.src_temp ./= src_areas
        dst_oc .= r.src_temp
    end
    return nothing
end

### `Interfacer.remap` extensions (Oceananigans Field ⟷ ClimaCore Field) with explicit regridding context

function Interfacer.remap!(target_field::OC.Field, source_field::CC.Fields.Field, remapping)
    regridding = remapping.regridding
    regridding === :conservative ||
        error(
            "remap!(::Oceananigans.Field, ::ClimaCore.Field, remapping) is only defined for `regridding === :conservative` (got $(repr(regridding))).",
        )
    CRX = cmip_conservative_regridding_cc_ext()
    CRX.get_value_per_element!(
        remapping.value_per_element_cc,
        source_field,
        remapping.field_ones_cc,
    )
    z = size(target_field, 3)
    dst = vec(OC.interior(target_field, :, :, z))
    src = remapping.value_per_element_cc
    w = remapping.conservative_weights
    _conservative_regrid_oc_from_cc!(dst, w.A, src, w.src_areas, remapping.remapper_oc_to_cc)
    return nothing
end

function Interfacer.remap!(target_field::CC.Fields.Field, source_field::OC.Field, remapping)
    regridding = remapping.regridding
    regridding === :conservative ||
        error("Unknown CMIP ocean regridding mode: $(repr(regridding))")
    z = size(source_field, 3)
    src = vec(OC.interior(source_field, :, :, z))
    dst = remapping.value_per_element_cc
    w = remapping.conservative_weights
    _conservative_regrid_cc_from_oc!(dst, w.A, src, w.dst_areas, remapping.remapper_oc_to_cc)
    CRX = cmip_conservative_regridding_cc_ext()
    CRX.set_value_per_element!(target_field, dst)
    return nothing
end

function Interfacer.remap!(
    target_field::CC.Fields.Field,
    operation::OC.AbstractOperations.AbstractOperation,
    remapping,
)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    Interfacer.remap!(target_field, evaluated_field, remapping)
    return nothing
end

function Interfacer.remap(target_space::CC.Spaces.AbstractSpace, source_field::OC.Field, remapping)
    if remapping.regridding === :spectral
        return Interfacer.remap(target_space, source_field)
    else
        target_field = CC.Fields.zeros(target_space)
        Interfacer.remap!(target_field, source_field, remapping)
        return target_field
    end
end

function Interfacer.remap(
    target_space::CC.Spaces.AbstractSpace,
    operation::OC.AbstractOperations.AbstractOperation,
    remapping,
)
    target_field = CC.Fields.zeros(target_space)
    Interfacer.remap!(target_field, operation, remapping)
    return target_field
end

Interfacer.remap!(target_field::CC.Fields.Field, source_field::CC.Fields.Field, _remapping) =
    Interfacer.remap!(target_field, source_field)

Interfacer.remap!(target_field::CC.Fields.Field, source_field::Number, _remapping) =
    Interfacer.remap!(target_field, source_field)

Interfacer.remap(target_space::CC.Spaces.AbstractSpace, source_num::Number, _remapping) =
    Interfacer.remap(target_space, source_num)
