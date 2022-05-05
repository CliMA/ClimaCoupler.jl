# slab_rhs!
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
function slab_space_init(; radius = 6371e3, Nel = 8, Nq = 2)

    # construct domain spaces - get only surface layer (NB: z should be zero, not z = first central height)
    space = ShellDomain(radius = radius, Nel = Nel, Nq = Nq)
    coords = ClimaCore.Fields.coordinate_field(space)

    # initial condition
    T_sfc = map(coords) do coord
        T_sfc = 281.0 # close to the average of T_1 in atmos
        # x_rnge = xmax - xmin
        # T_sfc = 270.0 + 5.0 * coord.x / x_rnge
    end

    # prognostic variable
    Y = ClimaCore.Fields.FieldVector(T_sfc = T_sfc)

    return Y, space
end

# ode
function slab_rhs!(Y, dY, Ya, t)
    """
    Slab ocean:
    ∂_t T_sfc = F_sfc + G
    """
    p, F_sfc = Ya
    dT_sfc = dY.T_sfc

    dY.T_sfc .+= F_sfc ./ (p.h .* p.ρ .* p.c)
end

struct SlabSimulation{P, Y, D, I}
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function slab_init(
    ::Type{FT},
    tspan; # (0.0,10.0)
    stepper = Euler(),
    nelements = 6,
    npolynomial = 4,
    dt = 0.02,
    saveat = 1.0e10,
) where {FT}

    params = ThermalSlabParameters(0.5, 1500.0, 800.0, 260.0)

    Y, space = slab_space_init(Nel = nelements, Nq = npolynomial + 1)

    Ya = (params = params, F_sfc = ClimaCore.Fields.zeros(space), cp_d = [FT(0.0)]) #auxiliary
    problem = OrdinaryDiffEq.ODEProblem(slab_rhs!, Y, tspan, Ya)
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)

    simulation = SlabSimulation(params, Y, space, integrator)
end
