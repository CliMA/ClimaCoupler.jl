### Helper functions to use ConservativeRemapping.jl with Oceananigans.jl
"""
    compute_cell_matrix(grid::Union{OC.OrthogonalSphericalShellGrid, OC.LatitudeLongitudeGrid})

Get a vector of vector of coordinate tuples, of the format expected by the
ConservativeRemapping.jl regridder.
"""
function compute_cell_matrix(
    grid::Union{OC.OrthogonalSphericalShellGrid, OC.LatitudeLongitudeGrid},
)
    Fx, Fy, _ = size(grid)
    # TODO is it ok to hardcode Center? Regridder is specifically for Center, Center fields so I think it's ok
    ℓx, ℓy = OC.Center(), OC.Center()

    if isnothing(ℓx) || isnothing(ℓy)
        error(
            "cell_matrix can only be computed for fields with non-nothing horizontal location.",
        )
    end

    arch = grid.architecture
    FT = eltype(grid)

    vertices_per_cell = 5 # convention: [sw, nw, ne, se, sw]
    ArrayType = OC.Architectures.array_type(arch)
    cell_matrix = ArrayType{Tuple{FT, FT}}(undef, vertices_per_cell, Fx * Fy)

    OC.Utils.launch!(
        arch,
        grid,
        (Fx, Fy),
        _compute_cell_matrix!,
        cell_matrix,
        Fx,
        ℓx,
        ℓy,
        grid,
    )

    return cell_matrix
end

flip(::OC.Face) = OC.Center()
flip(::OC.Center) = OC.Face()

left_index(i, ::OC.Center) = i
left_index(i, ::OC.Face) = i - 1
right_index(i, ::OC.Center) = i + 1
right_index(i, ::OC.Face) = i

@kernel function _compute_cell_matrix!(cell_matrix, Fx, ℓx, ℓy, grid)
    i, j = @index(Global, NTuple)

    vx = flip(ℓx)
    vy = flip(ℓy)

    isw = left_index(i, ℓx)
    jsw = left_index(j, ℓy)

    inw = left_index(i, ℓx)
    jnw = right_index(j, ℓy)

    ine = right_index(i, ℓx)
    jne = right_index(j, ℓy)

    ise = right_index(i, ℓx)
    jse = left_index(j, ℓy)

    xsw = OC.ξnode(isw, jsw, 1, grid, vx, vy, nothing)
    ysw = OC.ηnode(isw, jsw, 1, grid, vx, vy, nothing)

    xnw = OC.ξnode(inw, jnw, 1, grid, vx, vy, nothing)
    ynw = OC.ηnode(inw, jnw, 1, grid, vx, vy, nothing)

    xne = OC.ξnode(ine, jne, 1, grid, vx, vy, nothing)
    yne = OC.ηnode(ine, jne, 1, grid, vx, vy, nothing)

    xse = OC.ξnode(ise, jse, 1, grid, vx, vy, nothing)
    yse = OC.ηnode(ise, jse, 1, grid, vx, vy, nothing)

    linear_idx = i + (j - 1) * Fx
    @inbounds begin
        cell_matrix[1, linear_idx] = (xsw, ysw)
        cell_matrix[2, linear_idx] = (xnw, ynw)
        cell_matrix[3, linear_idx] = (xne, yne)
        cell_matrix[4, linear_idx] = (xse, yse)
        cell_matrix[5, linear_idx] = (xsw, ysw)
    end
end

### Extensions of Interfacer.jl functions for Oceananigans fields/grids
# Non-allocating ClimaCore -> Oceananigans remap
function Interfacer.remap!(dst_field::OC.Field, src_field::CC.Fields.Field, remapping)
    CC.Remapping.get_value_per_element!(
        remapping.value_per_element_cc,
        src_field,
        remapping.field_ones_cc,
    )

    # Get the index of the top level (surface); 1 for 2D fields, Nz for 3D fields
    z = size(dst_field, 3)
    dst = vec(OC.interior(dst_field, :, :, z))
    src = remapping.value_per_element_cc

    # Multiply by transpose of the matrix of intersection areas
    LA.mul!(dst, transpose(remapping.remapper_oc_to_cc.intersections), src)

    # Normalize by the destination (Oceananigans) element areas
    dst ./= remapping.remapper_oc_to_cc.src_areas # Oceananigans areas are source areas
    return nothing
end
# Allocating ClimaCore -> Oceananigans remap
function Interfacer.remap(
    src_field::CC.Fields.Field,
    remapping,
    dst_space::Union{OC.OrthogonalSphericalShellGrid, OC.LatitudeLongitudeGrid},
)
    dst_field = OC.Field{OC.Center, OC.Center, Nothing}(dst_space)
    remap!(dst_field, src_field, remapping)
    return dst_field
end

# Non-allocating Oceananigans -> ClimaCore remap
function Interfacer.remap!(dst_field::CC.Fields.Field, src_field::OC.Field, remapping)
    # Get the index of the top level (surface); 1 for 2D fields, Nz for 3D fields
    z = size(src_field, 3)
    # Store the remapped FV values in a vector of length equal to the number of elements in the target space
    dst = remapping.value_per_element_cc
    src = vec(OC.interior(src_field, :, :, z))
    LA.mul!(dst, remapping.remapper_oc_to_cc.intersections, src)

    # Normalize by the destination (ClimaCore) element areas
    dst ./= remapping.remapper_oc_to_cc.dst_areas # ClimaCore areas are destination areas

    # Convert the vector of remapped values to a ClimaCore Field with one value per element
    CC.Remapping.set_value_per_element!(dst_field, dst)
    return nothing
end
# Handle the case of remapping the area fraction field, which is a ClimaCore Field
Interfacer.remap!(dst_field::CC.Fields.Field, src_field::CC.Fields.Field, remapping) =
    Interfacer.remap!(dst_field, src_field)
Interfacer.remap!(dst_field::CC.Fields.Field, src_field::Number, remapping) =
    Interfacer.remap!(dst_field, src_field)
# Allocating Oceananigans -> ClimaCore remap
function Interfacer.remap(
    src_field::OC.Field,
    remapping,
    dst_space::CC.Spaces.AbstractSpace,
)
    dst_field = CC.Fields.zeros(dst_space)
    Interfacer.remap!(dst_field, src_field, remapping)
    return dst_field
end

# Handle the case of remapping a scalar number to a ClimaCore space
Interfacer.remap(num::Number, remapping, target_space::CC.Spaces.AbstractSpace) =
    Interfacer.remap(num, target_space)


function Interfacer.remap(
    operation::OC.AbstractOperations.AbstractOperation,
    remapping,
    target_space,
)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    return Interfacer.remap(evaluated_field, remapping, target_space)
end

function Interfacer.remap!(
    target_field,
    operation::OC.AbstractOperations.AbstractOperation,
    remapping,
)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    return Interfacer.remap!(target_field, evaluated_field, remapping)
end

"""
    construct_remappers(oc_grid, space_cc)

Given an Oceananigans LatitudeLongitudeGrid and a ClimaCore space, construct the
remappers needed to remap between the two grids in both directions.

Returns a remapper from the Oceananigans grid to the ClimaCore boundary space.
To regrid from Oceananigans to ClimaCore, use `LA.mul!(dest_vector, remapper_oc_to_cc, src_vector)`.
To regrid from ClimaCore to Oceananigans, use `LA.mul!(dest_vector, transpose(remapper_oc_to_cc), src_vector)`.
"""
function construct_remappers(grid_oc, space_cc)
    # Get the vector of polygons for Oceananigans and ClimaCore spaces
    vertices_oc = compute_cell_matrix(grid_oc.underlying_grid)
    vertices_cc = CC.Remapping.get_element_vertices(space_cc)

    remapper_oc_to_cc = CR.Regridder(vertices_cc, vertices_oc; normalize = false)

    # Create a field of ones on the boundary space so we can compute element areas
    field_ones_cc = CC.Fields.ones(space_cc)

    # Allocate a vector with length equal to the number of elements in the target space
    # To be used as a temp field for remapping
    value_per_element_cc = zeros(Float64, CC.Meshes.nelements(space_cc.grid.topology.mesh))

    # Construct two 2D Oceananigans Center/Center fields to use as scratch space while remapping
    scratch_field_oc1 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)
    scratch_field_oc2 = OC.Field{OC.Center, OC.Center, Nothing}(grid_oc)
    return (;
        remapper_oc_to_cc,
        field_ones_cc,
        value_per_element_cc,
        scratch_field_oc1,
        scratch_field_oc2,
    )
end
