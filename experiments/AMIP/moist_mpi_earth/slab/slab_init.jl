# slab_rhs!
using ClimaCore

struct ThermalSlabParameters# <: CLIMAParameters.AbstractEarthParameterSet{F} 
    h::FT
    ρ::FT
    c::FT
    T_init::FT
    z0m::FT
    z0b::FT
end

# domain
function ShellDomain(; radius = 6371e3, Nel = 8, Nq = 2)
    domain = ClimaCore.Domains.SphereDomain(radius)
    mesh = ClimaCore.Meshes.EquiangularCubedSphere(domain, Nel)
    topology = ClimaCore.Topologies.Topology2D(mesh)
    quad = ClimaCore.Spaces.Quadratures.GLL{Nq}()
    space = ClimaCore.Spaces.SpectralElementSpace2D(topology, quad)
end

struct SlabSimulation{P, Y, D, I}
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

# init simulation
function slab_space_init(::Type{FT}, space, p; anomaly = false, hs_sfc = false) where {FT}

    coords = ClimaCore.Fields.coordinate_field(space)

    # initial condition
    T_sfc = map(coords) do coord
        T_sfc_0 = FT(p.T_init)
        radlat = coord.lat / FT(180) * pi
        ΔT = FT(0)
        if anomaly == true
            anom_ampl = FT(0)
            lat_0 = FT(60) / FT(180) * pi
            lon_0 = FT(-90) / FT(180) * pi
            radlon = coord.long / FT(180) * pi
            stdev = FT(5) / FT(180) * pi
            ΔT = anom_ampl * exp(-((radlat - lat_0)^2 / 2stdev^2 + (radlon - lon_0)^2 / 2stdev^2))
        elseif hs_sfc == true
            ΔT = -FT(60) * sin(radlat)^2
        end
        T_sfc_0 + ΔT
    end

    # prognostic variable
    Y = ClimaCore.Fields.FieldVector(T_sfc = T_sfc)

    return Y, space
end

get_slab_energy(slab_sim, T_sfc) = slab_sim.params.ρ .* slab_sim.params.c .* T_sfc .* slab_sim.params.h

function slab_rhs!(dY, Y, Ya, t)
    """
    Slab:
    ∂_t T_sfc = F_aero + G
    """
    p, F_aero, F_rad, mask = Ya

    rhs = @. (F_aero + F_rad) / (p.h * p.ρ * p.c)
    parent(dY.T_sfc) .= apply_mask.(parent(mask), >, parent(rhs), FT(0))
end

function slab_init(
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

    params = ThermalSlabParameters(FT(1), FT(1500.0), FT(800.0), FT(315.0), FT(1e-3), FT(1e-5)) # T_init close to the average of T_1 in atmos

    Y, space = slab_space_init(FT, space, params, hs_sfc = true)
    Ya = (params = params, F_aero = ClimaCore.Fields.zeros(space), F_rad = ClimaCore.Fields.zeros(space), mask = mask) #auxiliary
    problem = OrdinaryDiffEq.ODEProblem(slab_rhs!, Y, tspan, Ya)
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)

    SlabSimulation(params, Y, space, integrator)
end
