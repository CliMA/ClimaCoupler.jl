# wireframe test
using GLMakie

x, y = collect(-8:0.5:8), collect(-8:0.5:8)
z = [sinc(√(X^2 + Y^2) / π) for X ∈ x, Y ∈ y]

wireframe(x, y, z, axis=(type=Axis3,), color=:black)

# Clima warp test - need the Project.toml in CC
using GeometryBasics
include("cubed_sphere_warp_helper.jl")

n = 8

r = range(-1,1,length=n+1)

xa = zeros(n+1,n+1)
xb = zeros(n+1,n+1)
xc = zeros(n+1,n+1)

xa .= r
xb .= r'
xc .= 1

a = [equiangular_cubed_shell_warp(a,b,c)[1] for (a,b,c) in zip(xa,xb,xc)]
b = [equiangular_cubed_shell_warp(a,b,c)[2] for (a,b,c) in zip(xa,xb,xc)]
c = [equiangular_cubed_shell_warp(a,b,c)[3] for (a,b,c) in zip(xa,xb,xc)]


sc = Scene(limits = GeometryBasics.HyperRectangle(GeometryBasics.Vec3f0(-2), GeometryBasics.Vec3f0(2)))

wireframe!(sc,
           a,b,c,
           show_axis=false,
           linewidth=1.2)
wireframe!(sc,
           a,b,-c,
           show_axis=false,
           linewidth=1.2)
wireframe!(sc,
           c,a,b,
           show_axis=false,
           linewidth=1.2)
wireframe!(sc,
           -c,a,b,
           show_axis=false,
           linewidth=1.2)
wireframe!(sc,
           b,c,a,
           show_axis=false,
           linewidth=1.2)
wireframe!(sc,
           b,-c,a,
           show_axis=false,
           linewidth=1.2)

# field plotting in ClimaCore
import ClimaCore
using IntervalSets
import ClimaCoreMakie: ClimaCoreMakie, Makie
R = 6.37122e6

domain = ClimaCore.Domains.SphereDomain(R)
mesh = ClimaCore.Meshes.EquiangularCubedSphere(domain, 6)
grid_topology = ClimaCore.Topologies.Topology2D(mesh)
quad = ClimaCore.Spaces.Quadratures.GLL{5}()
space = ClimaCore.Spaces.SpectralElementSpace2D(grid_topology, quad)
coords = ClimaCore.Fields.coordinate_field(space)

u = map(coords) do coord
    u0 = 20.0
    α0 = 45.0
    ϕ = coord.lat
    λ = coord.long

    uu = u0 * (cosd(α0) * cosd(ϕ) + sind(α0) * cosd(λ) * sind(ϕ))
    uv = -u0 * sind(α0) * sind(λ)
    ClimaCore.Geometry.UVVector(uu, uv)
end

fig = ClimaCoreMakie.viz(u.components.data.:1) 

lmesh = ClimaCoreMakie.level_mesh(u.components.data.:1,2).position

lmesh1 = parent([pt[1] for pt in lmesh])
lmesh2 = parent([pt[2] for pt in lmesh])
lmesh3 = parent([pt[3] for pt in lmesh])


# overlay mesh:
# wireframe!(sc,
#             lmesh1,lmesh2,lmesh3,
#            show_axis=false,
#            linewidth=1.2)
#Makie.mesh(lmesh1,lmesh2,lmesh3, colormap = Reverse(:Spectral)) # https://makie.juliaplots.org/stable/examples/plotting_functions/mesh/
i = 3456 # first face; NB: size(lmesh) = 20736 = 2 *  triangle_count = 2 *  Nh * (Nj - 1) * (Ni - 1) * 2 = 2 * (6*6*6) * (Nj - 1) * (Ni - 1) * 2 
Makie.linesegments(lmesh1[1:i],lmesh2[1:i],lmesh3[1:i], overdraw = true)
Makie.scatter(lmesh2[1:i],lmesh3[1:i], overdraw = true) # 3d scatter doesn't work - known bug - https://discourse.julialang.org/t/ploints-not-appearing-in-3d-scatter/64622


