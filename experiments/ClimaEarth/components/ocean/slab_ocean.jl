import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaUtilities
import ClimaUtilities.TimeManager: date
import ClimaCoupler: Checkpointer, FluxCalculator, Interfacer, Utilities, FieldExchanger

###
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###
"""
    SlabOceanSimulation{P, I}

Equation:

    (h * ρ * c) dTdt = (-F_turb_energy + (1 - α) * SW_d + LW_d - LW_u)

"""
struct SlabOceanSimulation{P, I} <: Interfacer.OceanModelSimulation
    params::P
    integrator::I
end

# ocean parameters
Base.@kwdef struct OceanSlabParameters{FT <: AbstractFloat}
    h::FT                   # depth of the ocean [m]
    ρ::FT                   # density of the ocean [kg / m3]
    c::FT                   # specific heat of the ocean [J / kg / K]
    T_init::FT              # initial temperature of the ocean [K]
    z0m::FT                 # roughness length for momentum [m]
    z0b::FT                 # roughness length for heat [m]
    α::FT                   # albedo of the ocean [0, 1]
    ϵ::FT                   # emissivity of the ocean
    σ::FT                   # Stefan-Boltzmann constant [W / m2 / K4]
    evolving_switch::FT     # switch to turn off the evolution of the ocean temperature [0 or 1]
end

"""
    OceanSlabParameters{FT}(coupled_param_dict; h = FT(20), ρ = FT(1500),
                            c = FT(800), T_init = FT(271), z0m = FT(5e-4),
                            z0b = FT(5e-4), α = FT(0.38), ϵ = FT(1),
                            evolving_switch = FT(1))

Initialize the `OceanSlabParameters` object with the coupled parameters.

# Arguments
- `coupled_param_dict`: a dictionary of coupled parameters (required)
- `h`: depth of the ocean [m] (default: 20)
- `ρ`: density of the ocean [kg / m3] (default: 1500)
- `c`: specific heat of the ocean [J / kg / K] (default: 800)
- `T_init`: initial temperature of the ocean [K] (default: 271)
- `z0m`: roughness length for momentum [m] (default: 5e-4)
- `z0b`: roughness length for heat [m] (default: 5e-4)
- `α`: albedo of the ocean [0, 1] (default: 0.38)
- `ϵ`: emissivity of the ocean (default: 1)
- `evolving_switch`: switch to turn off the evolution of the ocean temperature [0 or 1] (default: 1)

# Returns
- `OceanSlabParameters{FT}`: an `OceanSlabParameters` object
"""
function OceanSlabParameters{FT}(
    coupled_param_dict;
    h = FT(20),
    ρ = FT(1500),
    c = FT(800),
    T_init = FT(271),
    z0m = FT(5e-4),
    z0b = FT(5e-4),
    α = FT(0.38),
    ϵ = FT(1),
    evolving_switch = FT(1),
) where {FT}
    return OceanSlabParameters{FT}(;
        h,
        ρ,
        c,
        T_init,
        z0m,
        z0b,
        α,
        ϵ,
        σ = coupled_param_dict["stefan_boltzmann_constant"],
        evolving_switch,
    )
end

"""
    slab_ocean_space_init(space, params)

Initialize the slab ocean prognostic variable (temperature), including an
anomaly in the tropics by default.
"""
function slab_ocean_space_init(space, params)
    FT = CC.Spaces.undertype(space)
    coords = CC.Fields.coordinate_field(space)

    # initial condition
    # FT(271) close to the average of T_1 in atmos
    T_sfc = params.T_init .+ temp_anomaly.(coords)

    # prognostic variable
    Y = CC.Fields.FieldVector(; T_sfc = T_sfc)

    return Y
end

"""
    Interfacer.OceanSimulation(::Type{FT}, ::Val{:slab}; kwargs...)

Extension of the generic OceanSimulation constructor for the slab ocean model.
"""
function Interfacer.OceanSimulation(::Type{FT}, ::Val{:slab}; kwargs...) where {FT}
    return SlabOceanSimulation(FT; kwargs...)
end

"""
    SlabOceanSimulation(::Type{FT};
        tspan,
        dt,
        saveat,
        space,
        coupled_param_dict,
        thermo_params,
        stepper = CTS.RK4(),
        evolving = true,
    ) where {FT}

Initializes the `DiffEq` problem, and creates a Simulation-type object containing the
necessary information for `step!` in the coupling loop.
"""
function SlabOceanSimulation(
    ::Type{FT};
    tspan,
    dt,
    saveat,
    boundary_space,
    coupled_param_dict,
    thermo_params,
    stepper = CTS.RK4(),
    evolving = true,
    extra_kwargs...,
) where {FT}
    # Create params with evolving_switch override
    evolving_switch = evolving ? FT(1) : FT(0)
    params = OceanSlabParameters{FT}(coupled_param_dict; evolving_switch)

    Y = slab_ocean_space_init(boundary_space, params)
    cache = (
        params = params,
        F_turb_energy = CC.Fields.zeros(boundary_space),
        SW_d = CC.Fields.zeros(boundary_space),
        LW_d = CC.Fields.zeros(boundary_space),
        area_fraction = ones(boundary_space),
        thermo_params = thermo_params,
        α_direct = CC.Fields.ones(boundary_space) .* params.α,
        α_diffuse = CC.Fields.ones(boundary_space) .* params.α,
        u_atmos = CC.Fields.zeros(boundary_space),
        v_atmos = CC.Fields.zeros(boundary_space),
        # add dss_buffer to cache to avoid runtime dss allocation
        dss_buffer = CC.Spaces.create_dss_buffer(Y),
    )

    ode_algo = CTS.ExplicitAlgorithm(stepper)
    ode_function = CTS.ClimaODEFunction(;
        T_exp! = slab_ocean_rhs!,
        dss! = (Y, p, t) -> CC.Spaces.weighted_dss!(Y, p.dss_buffer),
    )
    if typeof(dt) isa Number
        dt = Float64(dt)
        tspan = Float64.(tspan)
        saveat = Float64.(saveat)
    end
    problem = SciMLBase.ODEProblem(ode_function, Y, tspan, cache)
    integrator =
        SciMLBase.init(problem, ode_algo, dt = dt, saveat = saveat, adaptive = false)

    sim = SlabOceanSimulation(params, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

# extensions required by Interfacer
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:area_fraction}) =
    sim.integrator.p.area_fraction
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:roughness_buoyancy}) =
    sim.integrator.p.params.z0b
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:roughness_momentum}) =
    sim.integrator.p.params.z0m
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_direct_albedo}) =
    sim.integrator.p.α_direct
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_diffuse_albedo}) =
    sim.integrator.p.α_diffuse
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_temperature}) =
    sim.integrator.u.T_sfc

"""
    Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:energy})

Extension of Interfacer.get_field to get the energy of the ocean.
It multiplies the the slab temperature by the heat capacity, density, and depth.
"""
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:energy}) =
    sim.integrator.p.params.ρ .* sim.integrator.p.params.c .* sim.integrator.u.T_sfc .*
    sim.integrator.p.params.h

function Interfacer.update_field!(
    sim::SlabOceanSimulation,
    ::Val{:area_fraction},
    field::CC.Fields.Field,
)
    sim.integrator.p.area_fraction .= field
    return nothing
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:SW_d}, field)
    Interfacer.remap!(sim.integrator.p.SW_d, field)
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:LW_d}, field)
    Interfacer.remap!(sim.integrator.p.LW_d, field)
end
function Interfacer.update_field!(
    sim::SlabOceanSimulation,
    ::Val{:turbulent_energy_flux},
    field,
)
    Interfacer.remap!(sim.integrator.p.F_turb_energy, field)
end
function Interfacer.update_field!(
    sim::SlabOceanSimulation,
    ::Val{:surface_direct_albedo},
    field::CC.Fields.Field,
)
    Interfacer.remap!(sim.integrator.p.α_direct, field)
end
function Interfacer.update_field!(
    sim::SlabOceanSimulation,
    ::Val{:surface_diffuse_albedo},
    field::CC.Fields.Field,
)
    Interfacer.remap!(sim.integrator.p.α_diffuse, field)
end

# extensions required by FieldExchanger
Interfacer.step!(sim::SlabOceanSimulation, t) =
    Interfacer.step!(sim.integrator, t - sim.integrator.t, true)

function FluxCalculator.update_turbulent_fluxes!(
    sim::SlabOceanSimulation,
    fields::NamedTuple,
)
    Interfacer.update_field!(sim, Val(:turbulent_energy_flux), fields.F_lh .+ fields.F_sh)
    return nothing
end

"""
    Checkpointer.get_model_prog_state(sim::SlabOceanSimulation)

Extension of Checkpointer.get_model_prog_state to get the model state.
"""
function Checkpointer.get_model_prog_state(sim::SlabOceanSimulation)
    return sim.integrator.u
end

###
### Slab ocean model-specific functions (not explicitly required by ClimaCoupler.jl)
###
# ode
function slab_ocean_rhs!(dY, Y, cache, t)
    params, F_turb_energy, SW_d, LW_d = cache
    (; α, ϵ, σ, h, ρ, c) = params
    rhs = @. (-F_turb_energy + (1 - α) * SW_d + ϵ * (LW_d - σ * Y.T_sfc^4)) / (h * ρ * c)

    # Zero out tendencies where there is no ocean, so that temperature remains constant there
    @. rhs = ifelse(cache.area_fraction ≈ 0, zero(rhs), rhs)

    # Note that the area fraction has already been applied to the fluxes,
    #  so we don't need to multiply by it here.
    @. dY.T_sfc = rhs * params.evolving_switch
end

"""
    temp_anomaly(coord)

Calculates a an anomaly to be added to the initial condition for temperature.
This default case includes only an anomaly at the tropics.
"""
function temp_anomaly(coord)
    # include tropics anomaly
    FT = eltype(coord)
    anom = FT(29 * exp(-coord.lat^2 / (2 * 26^2)))
    return anom
end

"""
    dss_state!(sim::SlabOceanSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on slab ocean model sims.
"""
dss_state!(sim::SlabOceanSimulation) =
    CC.Spaces.weighted_dss!(sim.integrator.u, sim.integrator.p.dss_buffer)
