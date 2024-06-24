# Approaches:
#  - Find index that's not matching
#  - output data to NC file, reload nad

import ClimaComms
ClimaComms.@import_required_backends
import ClimaCore: Domains, Meshes, Quadratures, Topologies, Spaces, Fields, DataLayouts

function make_space(comms_ctx, FT)
    # Set up parameters
    h_elem = 30
    radius = FT(6731e3)
    nh_poly = 3

    # Create space components
    domain = Domains.SphereDomain(radius)
    mesh = Meshes.EquiangularCubedSphere(domain, h_elem)
    quad = Quadratures.GLL{nh_poly + 1}()
    topology = Topologies.DistributedTopology2D(
        comms_ctx,
        mesh,
        Topologies.spacefillingcurve(mesh),
    )
    space = Spaces.SpectralElementSpace2D(
        topology,
        quad;
    )
    return space
end


# MRE:
FT = Float32
device = ClimaComms.device()
comms_ctx = ClimaComms.context(device)
ClimaComms.init(comms_ctx)
space = make_space(comms_ctx, FT)

# Create a field of 0s in top half, 1s in bottom half
f = Fields.ones(space)
ClimaComms.allowscalar(device) do
    parent(f) .= rand()
    # dims = size(parent(f))
    # m = dims[1]
    # parent(f)[1:(m ÷ 2), :, :, :] .= FT(0)
end

value₀ = FT(0)
value₁ = FT(1)

dl = Fields.field_values(f)
pf = parent(f)
fv = Fields.field_values(f)

local i
ClimaComms.allowscalar(device) do
	us = DataLayouts.universal_size(dl)
	ctu = CartesianIndices(map(x -> Base.OneTo(x), us))         # CartesianIndex((Ni,Nj,Nf,Nv,Nh))
	# ctp = CartesianIndices(map(x -> Base.OneTo(x), size(pf))) # CartesianIndex((Nv,Ni,Nj,Nf,Nh)) --> giving CartesianIndex((Ni, Nj, Nf/Nv, Nh)) instead
	i = findfirst(ctu) do I
		# IP = CartesianIndex((I.I[4],I.I[1],I.I[2],I.I[3],I.I[5]))
        IP = CartesianIndex((I.I[1],I.I[2],I.I[3],I.I[5]))
		fv[I]==value₀ && pf[IP]==value₁
	end
	@show i
	i
end


# Notes:
# land
# ocean
# ice

# f = land .+ ocean .+ ice

# pf = parent(land) .+ parent(ocean) .+ parent(ice) # any difference?
# pf = parent(land .+ ocean .+ ice)                 # any difference?

# cond = minimum(f) == value₀
# cond = minimum(parent(f)) == value₁
# bad_found = minimum(f) == value₀ && minimum(parent(f)) == value₁ && value₀ ≠ value₁

# Example:
# julia> size(Fields.field_values(cfield))
# (4, 4, 1, 64, 5400)

# julia> DataLayouts.universal_size(Fields.field_values(cfield))
# (4, 4, 16, 64, 5400)

# julia> size(parent(Fields.field_values(cfield)))
# (64, 4, 4, 16, 5400)


# Next step:
# g = similar(f)
# parent(g) .= CUDA.rand(...)
# ClimaComms.allowscalar(device) do
# 	g[i] = value₀
# end

# f = land .+ ocean .+ ice
# dl = Fields.field_values(f)
# pf = parent(land) .+ parent(ocean) .+ parent(ice)

# cond = minimum(g) == value₀
# cond = minimum(parent(g)) == value₁
