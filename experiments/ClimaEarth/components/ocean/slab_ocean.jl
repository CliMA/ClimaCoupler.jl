import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: BCReader, Checkpointer, FluxCalculator, Interfacer, Utilities

include("../slab_utils.jl")

###
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###
"""
    SlabOceanSimulation{P, Y, D, I}

Equation:

    (h * ρ * c) dTdt = -(F_turb_energy + F_radiative)

"""
struct SlabOceanSimulation{P, Y, D, I} <: Interfacer.OceanModelSimulation
    params::P
    Y_init::Y
    domain::D
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

Interfacer.name(::SlabOceanSimulation) = "SlabOceanSimulation"

"""
    slab_ocean_space_init(space, params)

Initialize the slab ocean prognostic variable (temperature), including an
anomaly in the tropics by default.
"""
function slab_ocean_space_init(space, params)
    FT = CC.Spaces.undertype(space)
    coords = CC.Fields.coordinate_field(space)

    # initial condition
    T_sfc = CC.Fields.zeros(space) .+ params.T_init # FT(271) close to the average of T_1 in atmos
    @. T_sfc += temp_anomaly(coords)

    # prognostic variable
    Y = CC.Fields.FieldVector(T_sfc = T_sfc)

    return Y, space
end

"""
    ocean_init(::Type{FT}; tspan, dt, saveat, space, area_fraction, stepper = CTS.RK4()) where {FT}

Initializes the `DiffEq` problem, and creates a Simulation-type object containing the necessary information for `step!` in the coupling loop.
"""
function ocean_init(
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
        dss_buffer = CC.Spaces.create_dss_buffer(CC.Fields.zeros(space)),
    )

    ode_algo = CTS.ExplicitAlgorithm(stepper)
    ode_function = CTS.ClimaODEFunction(T_exp! = slab_ocean_rhs!, dss! = weighted_dss_slab!)

    problem = SciMLBase.ODEProblem(ode_function, Y, Float64.(tspan), cache)
    integrator = SciMLBase.init(problem, ode_algo, dt = Float64(dt), saveat = Float64(saveat), adaptive = false)

    sim = SlabOceanSimulation(params, Y, space, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

# extensions required by Interfacer
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:air_density}) = sim.integrator.p.ρ_sfc
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:beta}) = convert(eltype(sim.integrator.u), 1.0)
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:roughness_buoyancy}) = sim.integrator.p.params.z0b
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:roughness_momentum}) = sim.integrator.p.params.z0m
Interfacer.get_field(sim::SlabOceanSimulation, ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}}) =
    sim.integrator.p.params.α
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_humidity}) = sim.integrator.p.q_sfc
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:surface_temperature}) = sim.integrator.u.T_sfc
Interfacer.get_field(sim::SlabOceanSimulation, ::Val{:water}) = nothing

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

# extensions required by FieldExchanger
Interfacer.step!(sim::SlabOceanSimulation, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::SlabOceanSimulation) = Interfacer.reinit!(sim.integrator)

# extensions required by FluxCalculator (partitioned fluxes)
function FluxCalculator.update_turbulent_fluxes_point!(
    sim::SlabOceanSimulation,
    fields::NamedTuple,
    colidx::CC.Fields.ColumnIndex,
)
    (; F_turb_energy) = fields
    @. sim.integrator.p.F_turb_energy[colidx] = F_turb_energy
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
"""
    scale_sst(SST::FT, _info)
Ensures that the space of the SST struct matches that of the land_fraction, and converts the units to Kelvin (N.B.: this is dataset specific)
"""
scale_sst(SST, _info) =
    (Utilities.swap_space!(axes(_info.land_fraction), SST) .+ BCReader.float_type_bcf(_info)(273.15))

# ode
function slab_ocean_rhs!(dY, Y, cache, t)
    p, F_turb_energy, F_radiative, area_fraction = cache
    FT = eltype(Y.T_sfc)
    rhs = @. -(F_turb_energy + F_radiative) / (p.h * p.ρ * p.c)
    @. dY.T_sfc = rhs * Regridder.binary_mask(area_fraction) * p.evolving_switch
    @. cache.q_sfc = TD.q_vap_saturation_generic.(cache.thermo_params, Y.T_sfc, cache.ρ_sfc, TD.Liquid())
end

"""
    temp_anomaly(coord)

Calculates a an anomaly to be added to the initial condition for temperature.
This default case includes only an anomaly at the tropics.
"""
function temp_anomaly(coord)
    # include tropics anomaly
    anom = FT(29 * exp(-coord.lat^2 / (2 * 26^2)))
    return anom
end

"""
    dss_state!(sim::SlabOceanSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on slab ocean model sims.
"""
function dss_state!(sim::SlabOceanSimulation)
    Y = sim.integrator.u
    p = sim.integrator.p
    for key in propertynames(Y)
        field = getproperty(Y, key)
        buffer = get_dss_buffer(axes(field), p)
        CC.Spaces.weighted_dss!(field, buffer)
    end

end
