# Conservative (finite-volume) regridding between the ClimaCore cubed-sphere exchange grid
# and Oceananigans `LatitudeLongitudeGrid`, via ConservativeRegridding.jl (see ClimaCoupler PR #1724).

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
    CR.regrid!(dst, transpose(remapping.remapper_oc_to_cc), src)
    return nothing
end

function Interfacer.remap!(target_field::CC.Fields.Field, source_field::OC.Field, remapping)
    regridding = remapping.regridding
    regridding === :conservative ||
        error("Unknown CMIP ocean regridding mode: $(repr(regridding))")
    z = size(source_field, 3)
    src = vec(OC.interior(source_field, :, :, z))
    dst = remapping.value_per_element_cc
    CR.regrid!(dst, remapping.remapper_oc_to_cc, src)
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
