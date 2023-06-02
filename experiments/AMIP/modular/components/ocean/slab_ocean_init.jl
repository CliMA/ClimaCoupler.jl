using ClimaCore

"""
    SlabOceanSimulation{P, Y, D, I}

Equation:

    (h * ρ * c) dTdt = -(F_aero + F_rad)

"""
struct SlabOceanSimulation{P, Y, D, I} <: SurfaceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

# ocean parameters
struct OceanSlabParameters{FT <: AbstractFloat}
    h::FT
    ρ::FT
    c::FT
    T_init::FT
    z0m::FT
    z0b::FT
    α::FT
end

# init simulation
function slab_ocean_space_init(::Type{FT}, space, p) where {FT}

    coords = ClimaCore.Fields.coordinate_field(space)

    # initial condition
    T_sfc = map(coords) do coord
        T_sfc_0 = FT(p.T_init) #- FT(275) # close to the average of T_1 in atmos
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
function slab_ocean_rhs!(dY, Y, cache, t)
    p, F_aero, F_rad, area_fraction = cache
    FT = eltype(Y.T_sfc)
    rhs = @. -(F_aero + F_rad) / (p.h * p.ρ * p.c)
    parent(dY.T_sfc) .= parent(rhs .* Regridder.binary_mask.(area_fraction, threshold = eps()))
end

"""
    ocean_init(::Type{FT}; tspan, dt, saveat, space, area_fraction, stepper = Euler()) where {FT}

Initializes the `DiffEq` problem, and creates a Simulation-type object containing the necessary information for `step!` in the coupling loop.
"""
function ocean_init(::Type{FT}; tspan, dt, saveat, space, area_fraction, stepper = Euler()) where {FT}

    params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5), FT(0.06))

    Y, space = slab_ocean_space_init(FT, space, params)
    cache = (
        params = params,
        F_aero = ClimaCore.Fields.zeros(space),
        F_rad = ClimaCore.Fields.zeros(space),
        area_fraction = area_fraction,
    )
    problem = OrdinaryDiffEq.ODEProblem(slab_ocean_rhs!, Y, tspan, cache)
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)

    SlabOceanSimulation(params, Y, space, integrator)
end

# file specific
"""
    clean_sst(SST::FT, _info)
Ensures that the space of the SST struct matches that of the land_fraction, and converts the units to Kelvin (N.B.: this is dataset specific)
"""
clean_sst(SST, _info) = (swap_space!(zeros(axes(_info.land_fraction)), SST) .+ float_type_bcf(_info)(273.15))

# required by Interfacer
get_field(sim::SlabOceanSimulation, ::Val{:surface_temperature}) = sim.integrator.u.T_sfc
get_field(sim::SlabOceanSimulation, ::Val{:roughness_momentum}) = sim.integrator.p.params.z0m
get_field(sim::SlabOceanSimulation, ::Val{:roughness_buoyancy}) = sim.integrator.p.params.z0b
get_field(sim::SlabOceanSimulation, ::Val{:beta}) = convert(eltype(sim.integrator.u), 1.0)
get_field(sim::SlabOceanSimulation, ::Val{:albedo}) = sim.integrator.p.params.α
get_field(sim::SlabOceanSimulation, ::Val{:area_fraction}) = sim.integrator.p.area_fraction

function update_field!(sim::SlabOceanSimulation, ::Val{:area_fraction}, field::Fields.Field)
    sim.integrator.p.area_fraction .= field
end
