import ClimaTimeSteppers as CTS
import ClimaCoupler.Interfacer: SeaIceModelSimulation, get_field, update_field!, name
import ClimaCoupler.FieldExchanger: step!, reinit!
import ClimaCoupler.FluxCalculator: update_turbulent_fluxes_point!
using ClimaCoupler: Regridder
import Thermodynamics as TD

include("../slab_utils.jl")

"""
    PrescribedIceSimulation{P, Y, D, I}

Ice concentration is prescribed, and we solve the following energy equation:

    (h * ρ * c) d T_sfc dt = -(F_turb_energy + F_radiativead) + F_conductive

    with
    F_conductive = k_ice (T_base - T_sfc) / (h)

    The surface temperature (`T_sfc`) is the prognostic variable which is being
    modified by turbulent aerodynamic (`F_turb_energy`) and radiative (`F_turb_energy`) fluxes,
    as well as a conductive flux that depends on the temperature difference
    across the ice layer (with `T_base` being prescribed).

In the current version, the sea ice has a prescribed thickness, and we assume that it is not
sublimating. That contribution has been zeroed out in the atmos fluxes.

"""
struct PrescribedIceSimulation{P, Y, D, I} <: SeaIceModelSimulation
    params::P
    Y_init::Y
    domain::D
    integrator::I
end
name(::PrescribedIceSimulation) = "PrescribedIceSimulation"

# sea-ice parameters
struct IceSlabParameters{FT <: AbstractFloat}
    h::FT # ice thickness [m]
    ρ::FT # density of sea ice [kg / m3]
    c::FT # specific heat of sea ice [J / kg / K]
    T_base::FT # temperature of sea water at the ice base
    z0m::FT # roughness length for momentum [m]
    z0b::FT # roughness length for tracers [m]
    T_freeze::FT # freezing point of sea water [K]
    k_ice::FT # thermal condictivity of ice [W / m / K] (less in HM71)
    α::FT # sea ice albedo
end
name(::IceSlabParameters) = "IceSlabParameters"

# init simulation
function slab_ice_space_init(::Type{FT}, space, p) where {FT}
    Y = Fields.FieldVector(T_sfc = ones(space) .* p.T_freeze)
    return Y
end

"""
    ice_rhs!(du, u, p, _)

Rhs method in the form as required by `ClimeTimeSteppers`, with the tendency vector `du`,
the state vector `u` and the parameter vector, `p`, as input arguments.

This sea-ice energy formulation follows [Holloway and Manabe 1971](https://journals.ametsoc.org/view/journals/mwre/99/5/1520-0493_1971_099_0335_socbag_2_3_co_2.xml?tab_body=pdf),
where sea-ice concentrations and thicknes are prescribed, and the model solves for temperature (curbed at the freezing point).
"""
function ice_rhs!(du, u, p, _)
    dY = du
    Y = u
    FT = eltype(dY)

    params = p.params
    F_turb_energy = p.F_turb_energy
    F_radiative = p.F_radiative
    area_fraction = p.area_fraction
    T_freeze = params.T_freeze

    F_conductive = @. params.k_ice / (params.h) * (params.T_base - Y.T_sfc) # fluxes are defined to be positive when upward
    rhs = @. (-F_turb_energy - F_radiative + F_conductive) / (params.h * params.ρ * params.c)

    # do not count tendencies that lead to temperatures above freezing, and mask out no-ice areas
    area_mask = Regridder.binary_mask.(area_fraction, threshold = eps(FT))
    unphysical = @. Regridder.binary_mask.(T_freeze - (Y.T_sfc + FT(rhs) * FT(p.dt)), threshold = FT(0)) .* area_mask
    parent(dY.T_sfc) .= parent(rhs .* unphysical)

    @. p.q_sfc = TD.q_vap_saturation_generic.(p.thermo_params, Y.T_sfc, p.ρ_sfc, TD.Ice())
end

"""
    ice_init(::Type{FT}; tspan, dt, saveat, space, ice_fraction, stepper = CTS.RK4()) where {FT}

Initializes the `DiffEq` problem, and creates a Simulation-type object containing the necessary information for `step!` in the coupling loop.
"""
function ice_init(::Type{FT}; tspan, saveat, dt, space, area_fraction, thermo_params, stepper = CTS.RK4()) where {FT}

    params = IceSlabParameters(FT(2), FT(900.0), FT(2100.0), FT(271.2), FT(1e-3), FT(1e-5), FT(271.2), FT(2.0), FT(0.8))

    Y = slab_ice_space_init(FT, space, params)
    additional_cache = (;
        F_turb_energy = ClimaCore.Fields.zeros(space),
        F_radiative = ClimaCore.Fields.zeros(space),
        q_sfc = ClimaCore.Fields.zeros(space),
        ρ_sfc = ClimaCore.Fields.zeros(space),
        area_fraction = area_fraction,
        dt = dt,
        thermo_params = thermo_params,
        # add dss_buffer to cache to avoid runtime dss allocation
        dss_buffer = ClimaCore.Spaces.create_dss_buffer(ClimaCore.Fields.zeros(space)),
    )

    ode_algo = CTS.ExplicitAlgorithm(stepper)
    ode_function = CTS.ClimaODEFunction(T_exp! = ice_rhs!, dss! = weighted_dss_slab!)

    problem = ODEProblem(ode_function, Y, Float64.(tspan), (; additional_cache..., params = params))
    integrator = init(problem, ode_algo, dt = Float64(dt), saveat = Float64(saveat), adaptive = false)

    sim = PrescribedIceSimulation(params, Y, space, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

# file-specific
"""
    clean_sic(SIC, _info)
Ensures that the space of the SIC struct matches that of the mask, and converts the units from area % to area fraction.
"""
clean_sic(SIC, _info) = swap_space!(zeros(axes(_info.land_fraction)), SIC) ./ float_type_bcf(_info)(100.0)

# setting that SIC < 0.5 is counted as ocean if binary remapping.
get_ice_fraction(h_ice::FT, mono::Bool, threshold = 0.5) where {FT} =
    mono ? h_ice : Regridder.binary_mask(h_ice, threshold = FT(threshold))

# extensions required by Interfacer
get_field(sim::PrescribedIceSimulation, ::Val{:surface_temperature}) = sim.integrator.u.T_sfc
get_field(sim::PrescribedIceSimulation, ::Val{:surface_humidity}) = sim.integrator.p.q_sfc
get_field(sim::PrescribedIceSimulation, ::Val{:roughness_momentum}) = sim.integrator.p.params.z0m
get_field(sim::PrescribedIceSimulation, ::Val{:roughness_buoyancy}) = sim.integrator.p.params.z0b
get_field(sim::PrescribedIceSimulation, ::Val{:beta}) = convert(eltype(sim.integrator.u), 1.0)
get_field(sim::PrescribedIceSimulation, ::Val{:albedo}) = sim.integrator.p.params.α
get_field(sim::PrescribedIceSimulation, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
get_field(sim::PrescribedIceSimulation, ::Val{:air_density}) = sim.integrator.p.ρ_sfc

function update_field!(sim::PrescribedIceSimulation, ::Val{:area_fraction}, field::Fields.Field)
    sim.integrator.p.area_fraction .= field
end

function update_field!(sim::PrescribedIceSimulation, ::Val{:turbulent_energy_flux}, field)
    parent(sim.integrator.p.F_turb_energy) .= parent(field)
end
function update_field!(sim::PrescribedIceSimulation, ::Val{:radiative_energy_flux}, field)
    parent(sim.integrator.p.F_radiative) .= parent(field)
end
function update_field!(sim::PrescribedIceSimulation, ::Val{:air_density}, field)
    parent(sim.integrator.p.ρ_sfc) .= parent(field)
end

# extensions required by FieldExchanger
step!(sim::PrescribedIceSimulation, t) = step!(sim.integrator, t - sim.integrator.t, true)
reinit!(sim::PrescribedIceSimulation) = reinit!(sim.integrator)

# extensions required by FluxCalculator (partitioned fluxes)
function update_turbulent_fluxes_point!(sim::PrescribedIceSimulation, fields::NamedTuple, colidx::Fields.ColumnIndex)
    (; F_turb_energy) = fields
    @. sim.integrator.p.F_turb_energy[colidx] = F_turb_energy
end

"""
    get_model_state_vector(sim::PrescribedIceSimulation)

Extension of Checkpointer.get_model_state_vector to get the model state.
"""
function get_model_state_vector(sim::PrescribedIceSimulation)
    return sim.integrator.u
end

"""
    get_field(sim::PrescribedIceSimulation, ::Val{:energy})

Extension of Interfacer.get_field to get the energy of the ocean.
It multiplies the the slab temperature by the heat capacity, density, and depth.
"""
get_field(sim::PrescribedIceSimulation, ::Val{:energy}) =
    sim.integrator.p.params.ρ .* sim.integrator.p.params.c .* sim.integrator.u.T_sfc .* sim.integrator.p.params.h

get_field(sim::PrescribedIceSimulation, ::Val{:water}) = nothing

"""
    dss_state!(sim::PrescribedIceSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on prescribed ice simulations.
"""
function dss_state!(sim::PrescribedIceSimulation)
    Y = sim.integrator.u
    p = sim.integrator.p
    for key in propertynames(Y)
        field = getproperty(Y, key)
        buffer = get_dss_buffer(axes(field), p)
        Spaces.weighted_dss!(field, buffer)
    end
end
