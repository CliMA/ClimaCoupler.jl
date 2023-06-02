# slab_rhs! (TODO: rm during clean up)
using ClimaCore

"""
    ThermalSlabParameters{FT<:AbstractFloat}

Parameters needed for a thermal slab model.

"""
struct ThermalSlabParameters{FT <: AbstractFloat}
    h::FT
    ρ::FT
    c::FT
    z0m::FT
    z0b::FT
    α::FT
end


"""
    slab_space_init

Creates the state vector and sets initial conditions.
"""
function slab_space_init(::Type{FT}, space, T_init; anomaly = false, hs_sfc = false) where {FT}

    coords = ClimaCore.Fields.coordinate_field(space)

    # initial condition
    T_sfc = map(coords) do coord
        T_sfc_0 = FT(T_init)
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

"""
    slab_rhs!(dY, Y, cache, t)

Computes the rhs of the slab model.
"""
function slab_rhs!(dY, Y, cache, t)
    p, F_aero, F_rad = cache
    FT = eltype(Y.T_sfc)
    rhs = @. -(F_aero + F_rad) / (p.h * p.ρ * p.c)
    parent(dY.T_sfc) .= parent(rhs)
end

"""
    slab_init

Initializes the slab simulation.
"""
function slab_init(::Type{FT}; tspan, dt, saveat, space, land_fraction, stepper = Euler()) where {FT}

    params = ThermalSlabParameters(FT(1), FT(1500.0), FT(800.0), FT(1e-3), FT(1e-5), FT(0.2))
    T_init = FT(315)
    Y, space = slab_space_init(FT, space, T_init, hs_sfc = true)
    cache = (
        params = params,
        F_aero = ClimaCore.Fields.zeros(space),
        F_rad = ClimaCore.Fields.zeros(space),
        land_fraction = land_fraction,
    )
    problem = OrdinaryDiffEq.ODEProblem(slab_rhs!, Y, tspan, cache)
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)

    SlabSimulation(params, Y, space, integrator)
end

"""
    get_land_energy(slab_sim::SlabSimulation, T_sfc)

Returns the energy of the slab; a method for the slab
when used as the land model.
"""
function get_land_energy(slab_sim::SlabSimulation, T_sfc)
    return get_slab_energy(slab_sim, T_sfc)
end

"""
    get_land_temp(slab_sim::SlabSimulation)

Returns the surface temperature of the slab;
a method for the slab
when used as the land model.
"""
function get_land_temp(slab_sim)
    return slab_sim.integrator.u.T_sfc
end

"""
    get_land_roughness(slab_sim::SlabSimulation)

Returns the roughness length parameters of the slab;
a method for the slab
when used as the land model.
"""
function get_land_roughness(slab_sim::SlabSimulation)
    return slab_sim.integrator.p.params.z0m, slab_sim.integrator.p.params.z0b
end

"""
   land_albedo(slab_sim::SlabSimulation)

Returns the surface albedo of the slab;
a method for the slab
when used as the land model.
"""
land_albedo(slab_sim::SlabSimulation) = slab_sim.params.α

"""
   land_beta(slab_sim::SlabSimulation)

Returns the beta factor of the slab;
a method for the slab
when used as the land model.
"""
land_beta(slab_sim::SlabSimulation) = slab_sim.params.β

"""
    get_land_q(slab_sim::SlabSimulation, atmos_sim, T_land, ρ_sfc)

Returns the surface specific humidity of the slab;
a method for the slab
when used as the land model.
"""
function get_land_q(slab_sim::SlabSimulation, atmos_sim, T_land, ρ_sfc)
    return TD.q_vap_saturation_generic.(atmos_sim.integrator.p.params, T_land, ρ_sfc, TD.Liquid())
end
