import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, FluxCalculator, Interfacer, Utilities

###
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###
"""
    SlabOceanSimulation{P, I}

Equation:

    (h * ρ * c) dTdt = -(F_turb_energy + F_radiative)

"""
struct SlabOceanSimulation{P, I} <: Interfacer.OceanModelSimulation
    params::P
    integrator::I
end
Interfacer.name(::SlabOceanSimulation) = "SlabOceanSimulation"

# ocean parameters
Base.@kwdef struct OceanSlabParameters{FT <: AbstractFloat}
    h::FT = 20              # depth of the ocean [m]
    ρ::FT = 1500            # density of the ocean [kg / m3]
    c::FT = 800             # specific heat of the ocean [J / kg / K]
    T_init::FT = 271        # initial temperature of the ocean [K]
    z0m::FT = 5e-4          # roughness length for momentum [m]
    z0b::FT = 5e-4          # roughness length for heat [m]
    α::FT = 0.38            # albedo of the ocean [0, 1]
    evolving_switch::FT = 1 # switch to turn off the evolution of the ocean temperature [0 or 1]
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

    return Y, space
end

"""
    SlabOceanSimulation(::Type{FT}; tspan, dt, saveat, space, area_fraction, stepper = CTS.RK4()) where {FT}

Initializes the `DiffEq` problem, and creates a Simulation-type object containing the necessary information for `step!` in the coupling loop.
"""
function SlabOceanSimulation(
    ::Type{FT};
    tspan,
    dt,
    saveat,
    space,
    area_fraction,
    thermo_params,
    stepper = CTS.RK4(),
    evolving = true,
) where {FT}

    evolving_switch = evolving ? FT(1) : FT(0)
    params = OceanSlabParameters{FT}(evolving_switch = evolving_switch)

    Y, space = slab_ocean_space_init(space, params)
    cache = (
        params = params,
        F_turb_energy = CC.Fields.zeros(space),
        F_radiative = CC.Fields.zeros(space),
        q_sfc = CC.Fields.zeros(space),
        ρ_sfc = CC.Fields.zeros(space),
        area_fraction = area_fraction,
        thermo_params = thermo_params,
        # add dss_buffer to cache to avoid runtime dss allocation
        dss_buffer = CC.Spaces.create_dss_buffer(Y),
        α_direct = CC.Fields.ones(space) .* params.α,
        α_diffuse = CC.Fields.ones(space) .* params.α,
    )

    ode_algo = CTS.ExplicitAlgorithm(stepper)
    ode_function =
        CTS.ClimaODEFunction(; T_exp! = slab_ocean_rhs!, dss! = (Y, p, t) -> CC.Spaces.weighted_dss!(Y, p.dss_buffer))
    if typeof(dt) isa Number
        dt = Float64(dt)
        tspan = Float64.(tspan)
        saveat = Float64.(saveat)
    end
    problem = SciMLBase.ODEProblem(ode_function, Y, tspan, cache)
    integrator = SciMLBase.init(problem, ode_algo, dt = dt, saveat = saveat, adaptive = false)

    sim = SlabOceanSimulation(params, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

# extensions required by Interfacer
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:roughness_buoyancy}) = sim.integrator.p.params.z0b
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:roughness_momentum}) = sim.integrator.p.params.z0m
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_direct_albedo}) = sim.integrator.p.α_direct
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_diffuse_albedo}) = sim.integrator.p.α_diffuse
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_temperature}) = sim.integrator.u.T_sfc

"""
    Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:energy})

Extension of Interfacer.get_field to get the energy of the ocean.
It multiplies the the slab temperature by the heat capacity, density, and depth.
"""
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:energy}) =
    sim.integrator.p.params.ρ .* sim.integrator.p.params.c .* sim.integrator.u.T_sfc .* sim.integrator.p.params.h

function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:area_fraction}, field::CC.Fields.Field)
    sim.integrator.p.area_fraction .= field
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:air_density}, field)
    parent(sim.integrator.p.ρ_sfc) .= parent(field)
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:radiative_energy_flux_sfc}, field)
    parent(sim.integrator.p.F_radiative) .= parent(field)
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:turbulent_energy_flux}, field)
    parent(sim.integrator.p.F_turb_energy) .= parent(field)
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:surface_direct_albedo}, field::CC.Fields.Field)
    sim.integrator.p.α_direct .= field
end
function Interfacer.update_field!(sim::SlabOceanSimulation, ::Val{:surface_diffuse_albedo}, field::CC.Fields.Field)
    sim.integrator.p.α_diffuse .= field
end

# extensions required by FieldExchanger
Interfacer.step!(sim::SlabOceanSimulation, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::SlabOceanSimulation) = Interfacer.reinit!(sim.integrator)

"""
Extend Interfacer.add_coupler_fields! to add the fields required for SlabOceanSimulation.

The fields added are:
- `:ρ_sfc` (for humidity calculation)
- `:F_radiative` (for radiation input)
"""
function Interfacer.add_coupler_fields!(coupler_field_names, ::SlabOceanSimulation)
    ocean_coupler_fields = [:ρ_sfc, :F_radiative]
    push!(coupler_field_names, ocean_coupler_fields...)
end

function FluxCalculator.update_turbulent_fluxes!(sim::SlabOceanSimulation, fields::NamedTuple)
    (; F_lh, F_sh) = fields
    @. sim.integrator.p.F_turb_energy = F_lh + F_sh
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
    p, F_turb_energy, F_radiative = cache
    rhs = @. -(F_turb_energy + F_radiative) / (p.h * p.ρ * p.c)

    # Note that the area fraction has already been applied to the fluxes,
    #  so we don't need to multiply by it here.
    @. dY.T_sfc = rhs * p.evolving_switch
    @. cache.q_sfc = TD.q_vap_saturation_generic.(cache.thermo_params, Y.T_sfc, cache.ρ_sfc, TD.Liquid())
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
dss_state!(sim::SlabOceanSimulation) = CC.Spaces.weighted_dss!(sim.integrator.u, sim.integrator.p.dss_buffer)
