import SciMLBase
import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaCoupler: Checkpointer, FluxCalculator, Interfacer, Utilities

###
### Functions required by ClimaCoupler.jl for a SurfaceModelSimulation
###
"""
    EisenmanIceSimulation{P, Y, D, I}

Thermodynamic 0-layer, based on the Semtner 1976 model and later refined by
Eisenmen 2009 and Zhang et al 2021.

Note that Eisenman sea ice assumes gray radiation, no snow coverage, and
PartitionedStateFluxes for the surface flux calculation.
"""
struct EisenmanIceSimulation{P, Y, D, I} <: Interfacer.SeaIceModelSimulation
    params_ice::P
    Y_init::Y
    domain::D
    integrator::I
end
Interfacer.name(::EisenmanIceSimulation) = "EisenmanIceSimulation"

Base.@kwdef struct EisenmanIceParameters{FT <: AbstractFloat}
    z0m::FT = 1e-3                  # roughness length for momentum [m]
    z0b::FT = 1e-5                  # roughness length for buoyancy [m]
    C0_base::FT = 120               # ice base transfer coefficient [W / m2 / K]
    T_base::FT = 273.16             # ice base temperature [K]
    L_ice::FT = 3e8                 # latent heat coefficient for ice [J / m3]
    T_freeze::FT = 273.16           # temperature at freezing point [K]
    k_ice::FT = 2                   # thermal conductivity of ice [W / m / K]
    α::FT = 0.70                    # albedo
    σ::FT = 5.67e-8                 # Stefan Boltzmann constant [W / m^2 / K^4]
end

Base.@kwdef struct EisenmanOceanParameters{FT <: AbstractFloat}
    h::FT = 1                       # mixed layer depth [m]
    ρ::FT = 1020                    # density of the mixed layer [kg / m^3]
    c::FT = 4000                    # mixed layer heat capacity [J / kg / K ]
    z0m::FT = 1e-3                  # roughness length for momentum [m]
    z0b::FT = 1e-5                  # roughness length for buoyancy [m]
    α::FT = 0.1                     # albedo
end

"""
    EisenmanIceSimulation(::Type{FT}, tspan; space = nothing, area_fraction = nothing, thermo_params = nothing, stepper = CTS.RK4(), dt = 0.02, saveat = 1.0e10)

Initialize the Eisenman-Zhang sea ice model and simulation.
"""
function EisenmanIceSimulation(
    ::Type{FT},
    tspan;
    space = nothing,
    area_fraction = nothing,
    thermo_params = nothing,
    stepper = CTS.RK4(),
    dt = 0.02,
    saveat = [1.0e10],
) where {FT}

    params_ice = EisenmanIceParameters{FT}()
    params_ocean = EisenmanOceanParameters{FT}()
    params = (; p_i = params_ice, p_o = params_ocean)

    # initiate prognostic variables
    Y, Ya = state_init(params_ice, space)

    ode_algo = CTS.ExplicitAlgorithm(stepper)
    ode_function =
        CTS.ClimaODEFunction(T_exp! = ∑tendencies, dss! = (Y, p, t) -> CC.Spaces.weighted_dss!(Y, p.dss_buffer))

    cache = (;
        Ya = Ya,
        Δt = dt,
        params = params,
        area_fraction = area_fraction,
        ice_area_fraction = zeros(space),
        thermo_params = thermo_params,
        dss_buffer = CC.Spaces.create_dss_buffer(Y),
    )
    problem = SciMLBase.ODEProblem(ode_function, Y, Float64.(tspan), cache)
    integrator = SciMLBase.init(problem, ode_algo, dt = Float64(dt), saveat = Float64.(saveat), adaptive = false)

    sim = EisenmanIceSimulation(params, Y, space, integrator)
    return sim
end

# extensions required by Interfacer
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:air_density}) = sim.integrator.p.Ya.ρ_sfc
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:area_fraction}) = sim.integrator.p.area_fraction
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:beta}) = convert(eltype(sim.integrator.u), 1.0)
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:roughness_buoyancy}) =
    @. sim.integrator.p.params.p_i.z0b * (sim.integrator.p.ice_area_fraction) +
       sim.integrator.p.params.p_o.z0b .* (1 - sim.integrator.p.ice_area_fraction)
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:roughness_momentum}) =
    @. sim.integrator.p.params.p_i.z0m * (sim.integrator.p.ice_area_fraction) +
       sim.integrator.p.params.p_o.z0m .* (1 - sim.integrator.p.ice_area_fraction)
Interfacer.get_field(sim::EisenmanIceSimulation, ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}}) =
    @. sim.integrator.p.params.p_i.α * (sim.integrator.p.ice_area_fraction) +
       sim.integrator.p.params.p_o.α .* (1 - sim.integrator.p.ice_area_fraction)
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:surface_humidity}) = sim.integrator.u.q_sfc
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:surface_temperature}) = sim.integrator.u.T_sfc
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:water}) = nothing

"""
    Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:energy})

Extension of Interfacer.get_field to get the energy of the ocean.
It is the sum of the heat content of the mixed layer, the heat content of the ice, the heat flux from the ocean below ice.
"""
function Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:energy})
    p_i = sim.integrator.p.params.p_i
    p_o = sim.integrator.p.params.p_o
    C0_base = p_i.C0_base
    T_base = p_i.T_base
    L_ice = p_i.L_ice
    T_freeze = p_i.T_freeze
    k_ice = p_i.k_ice
    ocean_qflux = sim.integrator.p.Ya.ocean_qflux

    cache = sim.integrator.p
    Δt = cache.Δt
    e_base = cache.Ya.e_base
    ocean_qflux = cache.Ya.ocean_qflux

    FT = eltype(sim.integrator.u)

    hρc_ml = p_o.h * p_o.ρ * p_o.c

    e_ml = @. p_o.h * p_o.ρ * p_o.c * sim.integrator.u.T_ml # heat
    e_ice = @. p_i.L_ice * sim.integrator.u.h_ice # phase
    e_qflux = @. ocean_qflux * FT(sim.integrator.t)

    return @. e_ml + e_ice + e_qflux + e_base
end

function Interfacer.update_field!(sim::EisenmanIceSimulation, ::Val{:air_density}, field)
    parent(sim.integrator.p.Ya.ρ_sfc) .= parent(field)
end
function Interfacer.update_field!(sim::EisenmanIceSimulation, ::Val{:area_fraction}, field::CC.Fields.Field)
    sim.integrator.p.area_fraction .= field
end
function Interfacer.update_field!(sim::EisenmanIceSimulation, ::Val{:∂F_turb_energy∂T_sfc}, field)
    sim.integrator.p.Ya.∂F_turb_energy∂T_sfc .= field
end
function Interfacer.update_field!(sim::EisenmanIceSimulation, ::Val{:radiative_energy_flux_sfc}, field)
    parent(sim.integrator.p.Ya.F_rad) .= parent(field)
end
function Interfacer.update_field!(sim::EisenmanIceSimulation, ::Val{:turbulent_energy_flux}, field)
    parent(sim.integrator.p.Ya.F_turb) .= parent(field)
end

# extensions required by FieldExchanger
Interfacer.step!(sim::EisenmanIceSimulation, t) = Interfacer.step!(sim.integrator, t - sim.integrator.t, true)
Interfacer.reinit!(sim::EisenmanIceSimulation) = Interfacer.reinit!(sim.integrator)

# extensions required by FluxCalculator (partitioned fluxes)
function FluxCalculator.update_turbulent_fluxes!(sim::EisenmanIceSimulation, fields::NamedTuple)
    (; F_turb_energy) = fields
    @. sim.integrator.p.Ya.F_turb = F_turb_energy
end

"""
    Checkpointer.get_model_prog_state(sim::EisenmanIceSimulation)

Extension of Checkpointer.get_model_prog_state to get the model state.
"""
function Checkpointer.get_model_prog_state(sim::EisenmanIceSimulation)
    return sim.integrator.u
end

"""
    FluxCalculator.differentiate_turbulent_fluxes!(sim::EisenmanIceSimulation, args)

Extension of FluxCalculator.differentiate_turbulent_fluxes! from FluxCalculator to get the turbulent fluxes.
"""
FluxCalculator.differentiate_turbulent_fluxes!(sim::EisenmanIceSimulation, args) =
    FluxCalculator.differentiate_turbulent_fluxes!(sim::EisenmanIceSimulation, args..., ΔT_sfc = 0.1)

"""
    differentiate_turbulent_fluxes(sim::Interfacer.SurfaceModelSimulation, thermo_params, input_args, fluxes, δT_sfc = 0.1)

Differentiates the turbulent fluxes in the surface model simulation `sim` with respect to the surface temperature,
using δT_sfc as the perturbation.
"""
function FluxCalculator.differentiate_turbulent_fluxes!(
    sim::EisenmanIceSimulation,
    thermo_params,
    input_args,
    fluxes;
    δT_sfc = 0.1,
)
    (; thermo_state_int, surface_params, surface_scheme) = input_args
    thermo_state_sfc_dT = FluxCalculator.surface_thermo_state(sim, thermo_params, thermo_state_int, δT_sfc = δT_sfc)
    input_args = merge(input_args, (; thermo_state_sfc = thermo_state_sfc_dT))

    # set inputs based on whether the surface_scheme is `MoninObukhovScheme` or `BulkScheme`
    inputs = surface_inputs(surface_scheme, input_args)

    # calculate the surface fluxes
    _, _, F_shf_δT_sfc, F_lhf_δT_sfc, _ = get_surface_fluxes!(inputs, surface_params)

    (; F_shf, F_lhf) = fluxes

    # calculate the derivative
    ∂F_turb_energy∂T_sfc = @. (F_shf_δT_sfc + F_lhf_δT_sfc - F_shf - F_lhf) / δT_sfc

    Interfacer.update_field!(sim, Val(:∂F_turb_energy∂T_sfc), ∂F_turb_energy∂T_sfc)
end

###
### Eisenman-Zhang sea ice model-specific functions (not explicitly required by ClimaCoupler.jl)
###
"""
    state_init(p::EisenmanIceParameters, space)

Initialize the state vectors for the Eisenman-Zhang sea ice model.
"""
function state_init(p::EisenmanIceParameters, space::CC.Spaces.AbstractSpace)
    Y = CC.Fields.FieldVector(
        T_sfc = ones(space) .* p.T_freeze,
        h_ice = zeros(space),
        T_ml = ones(space) .* 277,
        q_sfc = CC.Fields.zeros(space),
    )
    Ya = CC.Fields.FieldVector(
        F_turb = CC.Fields.zeros(space),
        ∂F_turb_energy∂T_sfc = CC.Fields.zeros(space),
        F_rad = CC.Fields.zeros(space),
        e_base = CC.Fields.zeros(space),
        ocean_qflux = CC.Fields.zeros(space),
        ρ_sfc = CC.Fields.zeros(space),
    )
    return Y, Ya
end

"""
    get_∂F_rad_energy∂T_sfc(T_sfc, p)

Calculate the derivative of the radiative flux with respect to the surface temperature.
"""
function get_∂F_rad_energy∂T_sfc(T_sfc, p)
    FT = eltype(T_sfc)
    @. FT(4) * (FT(1) - p.α) * p.σ * T_sfc^3
end

"""
    solve_eisenman_model!(Y, Ya, p, Δt::FT)

Solve the Eisenman-Zhang sea ice model for one timestep.
"""
function solve_eisenman_model!(Y, Ya, p, thermo_params, Δt)

    # model parameter sets
    (; p_i, p_o) = p

    C0_base = p_i.C0_base
    T_base = p_i.T_base
    L_ice = p_i.L_ice
    T_freeze = p_i.T_freeze
    k_ice = p_i.k_ice

    hρc_ml = p_o.h * p_o.ρ * p_o.c

    # prognostic
    (; T_sfc, h_ice, T_ml, q_sfc) = Y

    # auxiliary
    (; F_turb, ∂F_turb_energy∂T_sfc, F_rad, ocean_qflux, e_base) = Ya

    # local
    F_atm = @. F_turb + F_rad
    ∂F_atmo∂T_sfc = get_∂F_rad_energy∂T_sfc.(T_sfc, Ref(p_i)) .+ ∂F_turb_energy∂T_sfc

    # ice thickness and mixed layer temperature changes due to atmosphereic and ocean fluxes
    ice_covered = parent(h_ice)[1] > 0

    FT = eltype(T_ml)
    if ice_covered # ice-covered
        F_base = @. C0_base * (T_ml - T_base)
        ΔT_ml = @. -(F_base - ocean_qflux) * FT(Δt) / (hρc_ml)
        Δh_ice = @. (F_atm - F_base - ocean_qflux) * FT(Δt) / L_ice
        @. e_base .+= F_base * FT(Δt)
    else # ice-free
        ΔT_ml = @. -(F_atm - ocean_qflux) * FT(Δt) / (hρc_ml)
        Δh_ice = 0
    end

    # T_ml is not allowed to be below freezing
    frazil_ice_formation = parent(T_ml .+ ΔT_ml)[1] < T_freeze
    if frazil_ice_formation
        # Note that ice formation, which requires T_m < T_freeze, implies that Δh_ice increases
        Δh_ice = @. Δh_ice - (T_ml + ΔT_ml - T_freeze) * (hρc_ml) / L_ice
        ΔT_ml = @. T_freeze - T_ml
    end

    # adjust ocean temperature if transition to ice-free
    transition_to_icefree = (parent(h_ice)[1] > 0) & (parent(h_ice .+ Δh_ice)[1] <= 0)
    if transition_to_icefree
        ΔT_ml = @. ΔT_ml - (h_ice + Δh_ice) * L_ice / (hρc_ml)
        Δh_ice = @. -h_ice
    end

    # solve for T_sfc
    remains_ice_covered = (parent(h_ice .+ Δh_ice)[1] > 0)
    if remains_ice_covered
        # if ice covered, solve implicity (for now one Newton iteration: ΔT_s = - F(T_s) / dF(T_s)/dT_s )
        h = @. h_ice + Δh_ice
        F_conductive = @. k_ice / h * (T_base - T_sfc)
        numerator = @. -F_atm + F_conductive
        denominator = @. k_ice / h + ∂F_atmo∂T_sfc
        ΔT_sfc = @. numerator / denominator
        surface_melting = (parent(T_sfc .+ ΔT_sfc)[1] > T_freeze)
        if surface_melting
            ΔT_sfc = @. T_freeze - T_sfc # NB: T_sfc not storing energy
        end
        # surface is ice-covered, so update T_sfc as ice surface temperature
        T_sfc .+= ΔT_sfc
        # update surface humidity
        @. q_sfc = TD.q_vap_saturation_generic.(thermo_params, T_sfc, Ya.ρ_sfc, TD.Ice())
    else # ice-free, so update T_sfc as mixed layer temperature
        T_sfc .= T_ml .+ ΔT_ml
        # update surface humidity
        @. q_sfc = TD.q_vap_saturation_generic.(thermo_params, T_sfc, Ya.ρ_sfc, TD.Liquid())
    end

    Y.T_ml .+= ΔT_ml
    Y.h_ice .+= Δh_ice
    Y.T_sfc .= T_sfc
    Y.q_sfc .= q_sfc

    return Y, Ya
end

"""
    ∑tendencies(dY, Y, cache, _)

Calculate the tendencies for the Eisenman-Zhang sea ice model.
"""
function ∑tendencies(dY, Y, cache, _)
    FT = eltype(dY)
    Δt = cache.Δt
    Ya = cache.Ya
    thermo_params = cache.thermo_params
    p = cache.params

    @. dY.T_ml = -Y.T_ml / Δt
    @. dY.h_ice = -Y.h_ice / Δt
    @. dY.T_sfc = -Y.T_sfc / Δt
    @. dY.q_sfc = -Y.q_sfc / Δt

    solve_eisenman_model!(Y, Ya, p, thermo_params, Δt)

    # Get dY/dt
    @. dY.T_ml += Y.T_ml / Δt
    @. dY.h_ice += Y.h_ice / Δt
    @. dY.T_sfc += Y.T_sfc / Δt
    @. dY.q_sfc = -Y.q_sfc / Δt

    # update ice area fraction (binary mask for now)
    cache.ice_area_fraction .= Utilities.binary_mask.(Y.h_ice)
end
