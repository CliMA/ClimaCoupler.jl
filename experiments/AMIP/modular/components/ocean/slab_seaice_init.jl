"""
    PrescribedIceSimulation{P, Y, D, I}

Ice concentration is prescribed, and solving the following energy equation:

    F_conductive = k_ice (T_base - T_sfc) / (h)
    (h * ρ * c) dTdt = -(F_aero + F_rad) + F_conductive

"""
struct PrescribedIceSimulation{F, P, Y, D, I} <: SurfaceModelSimulation
    FT::F
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

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
    F_aero = p.F_aero
    F_rad = p.F_rad
    ice_fraction = p.area_fraction
    T_freeze = params.T_freeze

    F_conductive = @. params.k_ice / (params.h) * (params.T_base - Y.T_sfc) # fluxes are defined to be positive when upward
    rhs = @. (-F_aero - F_rad + F_conductive) / (params.h * params.ρ * params.c)

    # do not count tendencies that lead to temperatures above freezing, and mask out no-ice areas
    unphysical = @. Regridder.binary_mask.(T_freeze - (Y.T_sfc + FT(rhs) * p.dt), threshold = FT(0)) .*
       Regridder.binary_mask.(ice_fraction, threshold = eps())
    parent(dY.T_sfc) .= parent(rhs .* unphysical)
end

"""
    ice_init(::Type{FT}; tspan, dt, saveat, space, ice_fraction, stepper = Euler()) where {FT}

Initializes the `DiffEq` problem, and creates a Simulation-type object containing the necessary information for `step!` in the coupling loop.
"""
function ice_init(::Type{FT}; tspan, saveat, dt, space, ice_fraction, stepper = Euler()) where {FT}

    params = IceSlabParameters(FT(2), FT(900.0), FT(2100.0), FT(271.2), FT(1e-3), FT(1e-5), FT(271.2), FT(2.0), FT(0.8))

    Y = slab_ice_space_init(FT, space, params)
    additional_cache = (;
        F_aero = ClimaCore.Fields.zeros(space),
        F_rad = ClimaCore.Fields.zeros(space),
        dt = dt,
        area_fraction = ice_fraction,
    )

    problem = OrdinaryDiffEq.ODEProblem(ice_rhs!, Y, tspan, (; additional_cache..., params = params))
    integrator = OrdinaryDiffEq.init(problem, stepper, dt = dt, saveat = saveat)


    PrescribedIceSimulation(FT, params, Y, space, integrator)
end

# setting that SIC < 0.5 is counted as ocean if binary remapping.
get_ice_fraction(h_ice::FT, mono::Bool, threshold = 0.5) where {FT} =
    mono ? h_ice : Regridder.binary_mask(h_ice, threshold = FT(threshold))

function update_turbulent_fluxes_point!(sim::PrescribedIceSimulation, fields, colidx)
    (; F_shf, F_lhf) = fields
    @. sim.integrator.p.F_aero[colidx] = F_shf + F_lhf
end

function update!(sim::PrescribedIceSimulation, ::Val{:net_radiation}, field)
    @. sim.integrator.p.F_rad .= field
end



get_temperature(sim::PrescribedIceSimulation) = sim.integrator.u.T_sfc
get_z0m(sim::PrescribedIceSimulation) = sim.integrator.p.params.z0m
get_z0b(sim::PrescribedIceSimulation) = sim.integrator.p.params.z0b
get_beta(sim::PrescribedIceSimulation) = convert(eltype(sim.integrator.u), 1.0)
get_albedo(sim::PrescribedIceSimulation) = sim.integrator.p.params.α
get_area_fraction(sim::PrescribedIceSimulation) = sim.integrator.p.area_fraction
get_humidity_point(sim::PrescribedIceSimulation, thermo_params, T_sfc, ρ_sfc) = TD.q_vap_saturation_generic.(thermo_params, T_sfc, ρ_sfc, TD.Liquid()) # this assumes a saturated surface!!!

get_temperature_point(sim::PrescribedIceSimulation, colidx) = get_temperature(sim)[colidx]
get_z0m_point(sim::PrescribedIceSimulation, colidx) = get_z0m(sim)
get_z0b_point(sim::PrescribedIceSimulation, colidx) = get_z0b(sim)
get_beta_point(sim::PrescribedIceSimulation, colidx) = get_beta(sim)
get_albedo_point(sim::PrescribedIceSimulation, colidx) = get_albedo(sim)[colidx]
get_heat_transfer_coefficient_point(sim::PrescribedIceSimulation, colidx) = sim.integrator.p.params.Ch
get_drag_transfer_coefficient_point(sim::PrescribedIceSimulation, colidx) = sim.integrator.p.params.Cd

# file-specific (move!)
"""
    clean_sic(SIC, _info)
Ensures that the space of the SIC struct matches that of the mask, and converts the units from area % to area fraction.
"""
clean_sic(SIC, _info) = swap_space!(zeros(axes(_info.land_fraction)), SIC) ./ float_type_bcf(_info)(100.0)

issaturated(::PrescribedIceSimulation, q) = isnothing(q)

reinit!(sim::PrescribedIceSimulation) = reinit!(sim.integrator)