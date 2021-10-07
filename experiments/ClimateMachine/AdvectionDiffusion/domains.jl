using ClimateMachine, MPI
using ClimateMachine.Mesh.Grids
using ClimateMachine.Mesh.Topologies
using ClimateMachine.MPIStateArrays

import ClimateMachine.Mesh.Grids: DiscontinuousSpectralElementGrid
import ClimateMachine.Mesh.Topologies: StackedBrickTopology
export Interval, Periodic

import LinearAlgebra: ×
import Base: *, getindex

abstract type AbstractDomain end

struct IntervalDomain{S, B} <: AbstractDomain
    min::S
    max::S
    periodic::B

    function IntervalDomain(a, b, periodic)
        a, b = promote(a, b)

        return new{typeof(a), typeof(periodic)}(a, b, periodic)
    end
end

function Interval(a, b)
    return IntervalDomain(a, b, false)
end

function Periodic(a, b)
    return IntervalDomain(a, b, true)
end

function Base.show(io::IO, Ω::IntervalDomain)
    min = Ω.min
    max = Ω.max
    periodic = Ω.periodic

    printstyled(io, "[", color = 226)
    printstyled(io, "$min, $max", color = 7)
    if periodic
        printstyled(io, ")", color = 226)
    else
        printstyled(io, "]", color = 226)
    end

    return nothing
end

struct ProductDomain{T} <: AbstractDomain
    domains::T
end

function ×(a::IntervalDomain, b::IntervalDomain)
    return ProductDomain((a, b))
end

function ×(a::ProductDomain, b::IntervalDomain)
    return ProductDomain((a.domains..., b))
end

function ×(a::IntervalDomain, b::ProductDomain)
    return ProductDomain((a, b.domains...))
end

function *(a::AbstractDomain, b::AbstractDomain)
    return a × b
end

function Base.show(io::IO, Ω::ProductDomain)
    for (i, domain) in enumerate(Ω.domains)
        print(domain)

        if i != length(Ω.domains)
            printstyled(io, " × ", color = 118)
        end
    end

    return nothing
end


getindex(Ω::ProductDomain, i::Int) = Ω.domains[i]

ndims(Ω::IntervalDomain) = 1
ndims(Ω::ProductDomain) = +(ndims.(Ω.domains)...)

function periodicityof(Ω::ProductDomain)
    periodicity = ones(Bool, ndims(Ω))

    for i in 1:ndims(Ω)
        periodicity[i] = Ω[i].periodic
    end

    return Tuple(periodicity)
end


```
CubedSphere for GCM 
```

struct DeepSphericalShellDomain{S} <: AbstractDomain
    radius::S
    height::S

    function DeepSphericalShellDomain(; radius = nothing, height = nothing )
        radius, height = promote(radius, height)

        return new{typeof(radius)}(radius, height)
    end
end

struct OceanDomain{S} <: AbstractDomain
    radius::S
    depth::S

    function OceanDomain(; radius = nothing, depth = nothing )
        radius, height = promote(radius, depth)

        return new{typeof(radius)}(radius, depth)
    end
end


```
Sphere helper functions for GCM 
```
rad(x,y,z) = sqrt(x^2 + y^2 + z^2)
lat(x,y,z) = asin(z/rad(x,y,z)) # ϕ ∈ [-π/2, π/2] 
lon(x,y,z) = atan(y,x) # λ ∈ [-π, π) 

r̂ⁿᵒʳᵐ(x,y,z) = norm([x,y,z]) ≈ 0 ? 1 : norm([x, y, z])^(-1)
ϕ̂ⁿᵒʳᵐ(x,y,z) = norm([x,y,0]) ≈ 0 ? 1 : (norm([x, y, z]) * norm([x, y, 0]))^(-1)
λ̂ⁿᵒʳᵐ(x,y,z) = norm([x,y,0]) ≈ 0 ? 1 : norm([x, y, 0])^(-1)

r̂(x,y,z) = r̂ⁿᵒʳᵐ(x,y,z) * @SVector([x, y, z])
ϕ̂(x,y,z) = ϕ̂ⁿᵒʳᵐ(x,y,z) * @SVector [x*z, y*z, -(x^2 + y^2)]
λ̂(x,y,z) = λ̂ⁿᵒʳᵐ(x,y,z) * @SVector [-y, x, 0] 


"""
function DiscontinuousSpectralElementGrid(Ω::ProductDomain; elements = nothing, polynomialorder = nothing)
# Description 
Computes a DiscontinuousSpectralElementGrid as specified by a product domain
# Arguments
-`Ω`: A product domain object
# Keyword Arguments TODO: Add brickrange and topology as keyword arguments
-`elements`: A tuple of integers ordered by (Nx, Ny, Nz) for number of elements
-`polynomialorder`: A tupe of integers ordered by (npx, npy, npz) for polynomial order
-`FT`: floattype, assumed Float64 unless otherwise specified
-`mpicomm`: default = MPI.COMM_WORLD
-`array`: default = Array, but should generally be ArrayType
# Return 
A DiscontinuousSpectralElementGrid object
"""
function DiscontinuousSpectralElementGrid(
    Ω::ProductDomain{FT},
    elements,
    polynomialorder,
    mpicomm = MPI.COMM_WORLD,
    boundary = (1, 2),
    array = Array,
) where {FT}
    if elements == nothing
        error_message = "Please specify the number of elements as a tuple whose size is commensurate with the domain,"
        error_message = "e.g., a 3 dimensional domain would need a specification like elements = (10,10,10)."
        @error(error_message)
        return nothing
    end

    if polynomialorder == nothing
        error_message = "Please specify the polynomial order as a tuple whose size is commensurate with the domain,"
        error_message = "e.g., a 3 dimensional domain would need a specification like polynomialorder = (3,3,3)."
        @error(error_message)
        return nothing
    end

    dimension = ndims(Ω)

    if (dimension < 2) || (dimension > 3)
        error_message = "SpectralElementGrid only works with dimensions 2 or 3. "
        error_message *= "The current dimension is " * string(ndims(Ω))
        println("The domain is ", Ω)
        @error(error_message)
        return nothing
    end

    if ndims(Ω) != length(elements)
        @error("Specified too many elements for the dimension of the domain")
        return nothing
    end

    if ndims(Ω) != length(polynomialorder)
        @error("Specified too many polynomialorders for the dimension of the domain")
        return nothing
    end

    periodicity = periodicityof(Ω)
    tuple_ranges = []

    for i in 1:dimension
        push!(
            tuple_ranges,
            range(FT(Ω[i].min); length = elements[i] + 1, stop = FT(Ω[i].max)),
        )
    end

    brickrange = Tuple(tuple_ranges)
    if boundary == nothing
        boundary = (ntuple(j -> (1, 2), dimension - 1)..., (3, 4))
    end

    topology = StackedBrickTopology(
        mpicomm,
        brickrange;
        periodicity = periodicity,
        boundary = boundary,
    )

    grid = DiscontinuousSpectralElementGrid(
        topology,
        FloatType = FT,
        DeviceArray = array,
        polynomialorder = polynomialorder,
    )
    return grid
end

```
for CubedSphere
```
# using CLIMAParameters
# using CLIMAParameters.Planet: MSLP, R_d, day, grav, Omega, planet_radius
# struct EarthParameterSet <: AbstractEarthParameterSet end
# const param_set = EarthParameterSet()
# _a::Float64 = planet_radius(param_set)
# atmos_height::Float64 = 30e3

# ```
# grid = DiscontinuousSpectralElementGrid(Ω=, elements=, polynomialorder=)
# ```

function DiscontinuousSpectralElementGrid(
    Ω::DeepSphericalShellDomain{FT},
    elements,
    polynomialorder,
    mpicomm = MPI.COMM_WORLD,
    boundary = (1, 2),
    array = Array,
) where {FT}
    Rrange = grid1d(Ω.radius, Ω.radius + Ω.height, nelem = elements.vertical)

    topl = StackedCubedSphereTopology(
        mpicomm,
        elements.horizontal,
        Rrange,
        boundary = boundary, 
    )

    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = array,
        polynomialorder = (polynomialorder.horizontal, polynomialorder.vertical),
        meshwarp = ClimateMachine.Mesh.Topologies.equiangular_cubed_sphere_warp,
    )
    return grid
end
function DiscontinuousSpectralElementGrid(
    Ω::OceanDomain{FT},
    elements,
    polynomialorder,
    mpicomm = MPI.COMM_WORLD,
    boundary = (1, 2),
    array = Array,
) where {FT}
    Rrange = grid1d(Ω.radius + Ω.depth, Ω.radius, nelem = elements.vertical)

    topl = StackedCubedSphereTopology(
        mpicomm,
        elements.horizontal,
        Rrange,
        boundary = boundary, 
    )

    grid = DiscontinuousSpectralElementGrid(
        topl,
        FloatType = FT,
        DeviceArray = array,
        polynomialorder = (polynomialorder.horizontal, polynomialorder.vertical),
        meshwarp = ClimateMachine.Mesh.Topologies.equiangular_cubed_sphere_warp,
    )
    return grid
end


function is_surface(xc, yc, zc)
    # Sphere case - could dispatch on domain type, maybe?
    height_from_surface=(xc^2 + yc^2 + zc^2)^0.5 - planet_radius(param_set)
    return height_from_surface ≈ 0
end
# end of module
