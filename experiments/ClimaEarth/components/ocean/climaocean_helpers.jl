

"""
    surface_flux(f::OC.AbstractField)

Extract the top boundary conditions for the given field.
"""
function surface_flux(f::OC.AbstractField)
    top_bc = f.boundary_conditions.top
    if top_bc isa OC.BoundaryCondition{<:OC.BoundaryConditions.Flux}
        return top_bc.condition
    else
        return nothing
    end
end

function Interfacer.remap(operation::OC.AbstractOperations.AbstractOperation, target_space)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    return Interfacer.remap(evaluated_field, target_space)
end

"""
    set_from_extrinsic_vector!(vector, grid, u_cc, v_cc)

Given the extrinsic vector components `u_cc` and `v_cc` as `Center, Center`
fields, rotate them onto the target grid and remap to `Face, Center` and
`Center, Face` fields, respectively.
"""
function set_from_extrinsic_vector!(vector, grid, u_cc, v_cc)
    arch = OC.Architectures.architecture(grid)

    # Rotate vector components onto the grid
    OC.Utils.launch!(arch, grid, :xy, _rotate_vector!, u_cc, v_cc, grid)

    # Fill halo regions with the rotated vector components so we can use them to interpolate
    OC.fill_halo_regions!(u_cc)
    OC.fill_halo_regions!(v_cc)

    # Interpolate the vector components to face/center and center/face respectively
    OC.Utils.launch!(
        arch,
        grid,
        :xy,
        _interpolate_vector!,
        vector.u,
        vector.v,
        grid,
        u_cc,
        v_cc,
    )
    return nothing
end

"""
    _rotate_vector!(τx, τy, grid)

Rotate the velocities from the extrinsic coordinate system to the intrinsic
coordinate system.
"""
@kernel function _rotate_vector!(τx, τy, grid)
    # Use `k = 1` to index into the reduced Fields
    i, j = @index(Global, NTuple)
    # Rotate u, v from extrinsic to intrinsic coordinate system
    τxr, τyr = OC.Operators.intrinsic_vector(i, j, 1, grid, τx, τy)
    @inbounds begin
        τx[i, j, 1] = τxr
        τy[i, j, 1] = τyr
    end
end

"""
    _interpolate_vector!(τx, τy, grid, τx_cc, τy_cc)

Interpolate the input fluxes `τx_cc` and `τy_cc`, which are Center/Center
Fields to Face/Center and Center/Face coordinates, respectively.
"""
@kernel function _interpolate_vector!(τx, τy, grid, τx_cc, τy_cc)
    # Use `k = 1` to index into the reduced Fields
    i, j = @index(Global, NTuple)
    @inbounds begin
        τx[i, j, 1] = OC.Operators.ℑxᶠᵃᵃ(i, j, 1, grid, τx_cc)
        τy[i, j, 1] = OC.Operators.ℑyᵃᶠᵃ(i, j, 1, grid, τy_cc)
    end
end

### Helper functions to use ConservativeRemapping.jl with Oceananigans.jl
"""
    compute_cell_matrix(field::AbstractField)

Get a vector of vector of coordinate tuples, of the format expected by the
ConservativeRemapping.jl regridder.
"""
function compute_cell_matrix(field::AbstractField)
    Fx, Fy, _ = size(field)
    LX, LY, _ = Oceananigans.Fields.location(field)
    ℓx, ℓy = LX(), LY()

    if isnothing(ℓx) || isnothing(ℓy)
        error(
            "cell_matrix can only be computed for fields with non-nothing horizontal location.",
        )
    end

    grid = field.grid
    arch = grid.architecture
    FT = eltype(grid)

    vertices_per_cell = 5 # convention: [sw, nw, ne, se, sw]
    ArrayType = Oceananigans.Architectures.array_type(arch)
    cell_matrix = ArrayType{Tuple{FT, FT}}(undef, vertices_per_cell, Fx * Fy)

    arch = grid.architecture
    Oceananigans.Utils.launch!(
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

flip(::Face) = Center()
flip(::Center) = Face()

left_index(i, ::Center) = i
left_index(i, ::Face) = i - 1
right_index(i, ::Center) = i + 1
right_index(i, ::Face) = i

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

    xsw = ξnode(isw, jsw, 1, grid, vx, vy, nothing)
    ysw = ηnode(isw, jsw, 1, grid, vx, vy, nothing)

    xnw = ξnode(inw, jnw, 1, grid, vx, vy, nothing)
    ynw = ηnode(inw, jnw, 1, grid, vx, vy, nothing)

    xne = ξnode(ine, jne, 1, grid, vx, vy, nothing)
    yne = ηnode(ine, jne, 1, grid, vx, vy, nothing)

    xse = ξnode(ise, jse, 1, grid, vx, vy, nothing)
    yse = ηnode(ise, jse, 1, grid, vx, vy, nothing)

    linear_idx = i + (j - 1) * Fx
    @inbounds begin
        cell_matrix[1, linear_idx] = (xsw, ysw)
        cell_matrix[2, linear_idx] = (xnw, ynw)
        cell_matrix[3, linear_idx] = (xne, yne)
        cell_matrix[4, linear_idx] = (xse, yse)
        cell_matrix[5, linear_idx] = (xsw, ysw)
    end
end
