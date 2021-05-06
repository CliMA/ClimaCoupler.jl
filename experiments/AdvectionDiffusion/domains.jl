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
Box domain and grid
```

name_it(Ne::NamedTuple{(:x, :y, :z)}) = Ne
name_it(Ne) = (x = Ne[1], y = Ne[2], z = Ne[3])
struct RectangularDomain{FT} <: AbstractDomain
    Np::Int
    Ne::NamedTuple{(:x, :y, :z), NTuple{3, Int}}
    L::NamedTuple{(:x, :y, :z), NTuple{3, FT}}
    x::NTuple{2, FT}
    y::NTuple{2, FT}
    z::NTuple{2, FT}
    periodicity::NamedTuple{(:x, :y, :z), NTuple{3, Bool}}
end
function RectangularDomain(
    FT = Float64;
    Ne,
    Np,
    x::Tuple{<:Number, <:Number},
    y::Tuple{<:Number, <:Number},
    z::Tuple{<:Number, <:Number},
    periodicity = (true, true, false),
)

    Ne = name_it(Ne)
    periodicity = name_it(periodicity)

    west, east = FT.(x)
    south, north = FT.(y)
    bottom, top = FT.(z)

    east > west || error("Domain x-limits must be increasing!")
    north > south || error("Domain y-limits must be increasing!")
    top > bottom || error("Domain z-limits must be increasing!")

    L = (x = east - west, y = north - south, z = top - bottom)

    return RectangularDomain(
        Np,
        Ne,
        L,
        (west, east),
        (south, north),
        (bottom, top),
        periodicity,
    )
end


Base.eltype(::RectangularDomain{FT}) where {FT} = FT

function DiscontinuousSpectralElementGrid(
    domain::RectangularDomain{FT};
    boundary_tags = ((0, 0), (0, 0), (1, 2)),
    array= Array,
    mpicomm = MPI.COMM_WORLD,
) where {FT}

    west, east = domain.x
    south, north = domain.y
    bottom, top = domain.z

    element_coordinates = (
        range(west, east, length = domain.Ne.x + 1),
        range(south, north, length = domain.Ne.y + 1),
        range(bottom, top, length = domain.Ne.z + 1),
    )

    topology = StackedBrickTopology(
        mpicomm,
        element_coordinates;
        periodicity = tuple(domain.periodicity...),
        boundary = boundary_tags,
    )

    grid = DiscontinuousSpectralElementGrid(
        topology,
        FloatType = FT,
        DeviceArray = array,
        polynomialorder = domain.Np,
    )

    return grid
end


```
CubedSphere domain and grid 
```

struct DeepSphericalShellDomain{S} <: AbstractDomain
    radius::S
    height::S

    function DeepSphericalShellDomain(; radius = nothing, height = nothing )
        radius, height = promote(radius, height)

        return new{typeof(radius)}(radius, height)
    end
end

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

```
Sphere helper functions
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


function is_surface(xc, yc, zc)
    # Sphere case - could dispatch on domain type, maybe?
    height_from_surface=(xc^2 + yc^2 + zc^2)^0.5 - planet_radius(param_set)
    return height_from_surface ≈ 0
end
# end of module
