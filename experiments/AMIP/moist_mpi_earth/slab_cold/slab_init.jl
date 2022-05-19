# slab_rhs!

# TODO: rename slab cold
using ClimaCore

struct ThermalSlabParameters# <: CLIMAParameters.AbstractEarthParameterSet{F} 
    h::FT
    ρ::FT
    c::FT
    T_init::FT
end

# domain
function ShellDomain(; radius = 6371e3, Nel = 8, Nq = 2)
    domain = ClimaCore.Domains.SphereDomain(radius)
    mesh = ClimaCore.Meshes.EquiangularCubedSphere(domain, Nel)
    topology = ClimaCore.Topologies.Topology2D(mesh)
    quad = ClimaCore.Spaces.Quadratures.GLL{Nq}()
    space = ClimaCore.Spaces.SpectralElementSpace2D(topology, quad)
end

# init simulation
function slab_cold_space_init(::Type{FT}; radius = 6371e3, Nel = 8, Nq = 2, space = nothing) where {FT}

    # construct domain spaces - get only surface layer (NB: z should be zero, not z = first central height)
    # space = ShellDomain(radius = radius, Nel = Nel, Nq = Nq)
    coords = ClimaCore.Fields.coordinate_field(space)

    # initial condition
    T_sfc = map(coords) do coord
        T_sfc_0 = FT(275.0) #- FT(275) # close to the average of T_1 in atmos
        anom_ampl = FT(0)
        radlat = coord.lat / FT(180) * pi
        lat_0 = FT(60) / FT(180) * pi
        lon_0 = FT(-90) / FT(180) * pi
        radlon = coord.long / FT(180) * pi
        stdev = FT(5) / FT(180) * pi
        anom = anom_ampl * exp(-((radlat - lat_0)^2 / 2stdev^2 + (radlon - lon_0)^2 / 2stdev^2))
        T_sfc = T_sfc_0 + anom
    end

    # prognostic variable
    Y = ClimaCore.Fields.FieldVector(T_sfc = T_sfc)

    return Y, space
end

# ode
function slab_cold_rhs!(dY, Y, Ya, t)
    """
    Slab ocean:
    ∂_t T_sfc = F_sfc + G
    """
    p, F_sfc, F_rad, mask = Ya

    rhs = @. (F_sfc + F_rad) / (p.h * p.ρ * p.c)
    parent(dY.T_sfc) .= apply_mask_.(parent(mask), parent(rhs), FT(0)) 
end

apply_mask_(mask, yes, no) = mask < 0.5 ? yes : no


struct SlabSimulation{P, Y, D, I}
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function slab_cold_init(
    ::Type{FT},
    tspan;
    stepper = Euler(),
    nelements = 6,
    npolynomial = 4,
    dt = 0.02,
    saveat = 1.0e10,
    space = nothing,
    mask = nothing,
) where {FT}

    params = ThermalSlabParameters(40, 1500.0, 800.0, 260.0)

    Y, space = slab_space_init(FT, space = space)
    Ya = (params = params, F_sfc = ClimaCore.Fields.zeros(space), F_rad = ClimaCore.Fields.zeros(space), mask = mask) #auxiliary
    problem = OrdinaryDiffEq.ODEProblem(slab_cold_rhs!, Y, tspan, Ya)
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)

    SlabSimulation(params, Y, space, integrator)
end

get_slab_energy(slab_sim, T_sfc) = slab_sim.params.ρ .* slab_sim.params.c .* T_sfc .* slab_sim.params.h
