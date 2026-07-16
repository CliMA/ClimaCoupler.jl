import ClimaCore as CC
import ClimaTimeSteppers as CTS
import ClimaComms
import ..Checkpointer, ..FluxCalculator, ..Interfacer, ..Utilities

"""
    EisenmanIceSimulation{P, I}

Thermodynamic 0-layer sea ice model over a slab ocean mixed layer, based on the
Semtner (1976) model and later refined by Eisenman & Wettlaufer (2009) and
Zhang et al. (2021).

Prognostic variables are the ice thickness `h_ice`, the mixed layer temperature
`T_ml`, the surface temperature `T_sfc`, and the accumulated basal energy
`e_base` (used for energy bookkeeping). Ice cover is a binary mask
(ice-covered wherever `h_ice > 0`).

Note that the Eisenman sea ice model assumes no snow coverage.
"""
struct EisenmanIceSimulation{P, I} <: Interfacer.AbstractSeaIceSimulation
    params::P
    integrator::I
end

Base.@kwdef struct EisenmanIceParameters{FT <: AbstractFloat}
    z0m::FT         # roughness length for momentum [m]
    z0b::FT         # roughness length for buoyancy [m]
    C0_base::FT     # ice base transfer coefficient [W / m2 / K]
    T_base::FT      # ice base temperature [K]
    L_ice::FT       # latent heat coefficient for ice [J / m3]
    T_freeze::FT    # temperature at freezing point [K]
    k_ice::FT       # thermal conductivity of ice [W / m / K]
    α::FT           # albedo of ice
    ϵ::FT           # emissivity of ice
    σ::FT           # Stefan-Boltzmann constant [W / m2 / K4]
end

"""
    EisenmanIceParameters{FT}(coupled_param_dict; kwargs...)

Initialize the `EisenmanIceParameters` object with the coupled parameters.

# Arguments
- `coupled_param_dict`: a dictionary of coupled parameters (required)
- `z0m`: roughness length for momentum [m] (default: 1e-3)
- `z0b`: roughness length for buoyancy [m] (default: 1e-5)
- `C0_base`: ice base transfer coefficient [W / m2 / K] (default: 120)
- `T_base`: ice base temperature [K] (default: 273.16)
- `L_ice`: latent heat coefficient for ice [J / m3] (default: 3e8)
- `T_freeze`: temperature at freezing point [K] (default: 273.16)
- `k_ice`: thermal conductivity of ice [W / m / K] (default: 2)
- `α`: albedo of ice (default: 0.70)
- `ϵ`: emissivity of ice (default: 1)
"""
function EisenmanIceParameters{FT}(
    coupled_param_dict;
    z0m = FT(1e-3),
    z0b = FT(1e-5),
    C0_base = FT(120),
    T_base = FT(273.16),
    L_ice = FT(3e8),
    T_freeze = FT(273.16),
    k_ice = FT(2),
    α = FT(0.70),
    ϵ = FT(1),
) where {FT}
    return EisenmanIceParameters{FT}(;
        z0m,
        z0b,
        C0_base,
        T_base,
        L_ice,
        T_freeze,
        k_ice,
        α,
        ϵ,
        σ = coupled_param_dict["stefan_boltzmann_constant"],
    )
end

Base.@kwdef struct EisenmanOceanParameters{FT <: AbstractFloat}
    h::FT = 1                       # mixed layer depth [m]
    ρ::FT = 1020                    # density of the mixed layer [kg / m3]
    c::FT = 4000                    # mixed layer heat capacity [J / kg / K]
    z0m::FT = 1e-3                  # roughness length for momentum [m]
    z0b::FT = 1e-5                  # roughness length for buoyancy [m]
    α::FT = 0.1                     # albedo of the ice-free ocean
end

"""
    Interfacer.SeaIceSimulation(::Type{FT}, ::Val{:eisenman}; kwargs...)

Extension of the generic SeaIceSimulation constructor for the Eisenman-Zhang
sea ice model.
"""
function Interfacer.SeaIceSimulation(::Type{FT}, ::Val{:eisenman}; kwargs...) where {FT}
    return EisenmanIceSimulation(FT; kwargs...)
end

"""
    EisenmanIceSimulation(
        ::Type{FT};
        tspan,
        dt,
        saveat,
        boundary_space,
        coupled_param_dict,
        stepper = CTS.RK4(),
        extra_kwargs...,
    ) where {FT}

Initialize the Eisenman-Zhang sea ice model and simulation.
"""
function EisenmanIceSimulation(
    ::Type{FT};
    tspan,
    dt,
    saveat,
    boundary_space,
    coupled_param_dict,
    stepper = CTS.RK4(),
    extra_kwargs...,
) where {FT}
    params_ice = EisenmanIceParameters{FT}(coupled_param_dict)
    params_ocean = EisenmanOceanParameters{FT}()
    params = (; p_i = params_ice, p_o = params_ocean)

    # initiate prognostic variables and forcing fields
    Y, Ya = eisenman_state_init(params_ice, boundary_space)

    cache = (;
        Ya...,
        params = params,
        area_fraction = ones(boundary_space),
        Δt = dt,
        # scratch copy of the state, used to evaluate the discrete solve
        # without mutating the stepper's state (see `eisenman_ice_rhs!`)
        Y_scratch = similar(Y),
        # add dss_buffer to cache to avoid runtime dss allocation
        dss_buffer = Utilities.init_dss_buffer(Y),
    )

    ode_algo = CTS.ExplicitAlgorithm(stepper)
    ode_function = CTS.ClimaODEFunction(
        T_exp! = eisenman_ice_rhs!,
        dss! = (Y, p, t) -> Utilities.apply_dss!(Y, p.dss_buffer),
    )
    if dt isa Number
        dt = Float64(dt)
        tspan = Float64.(tspan)
        saveat = Float64.(saveat)
    end
    problem = CTS.ODEProblem(ode_function, Y, tspan, cache)
    integrator = CTS.init(problem, ode_algo; dt, saveat, adaptive = false)

    sim = EisenmanIceSimulation(params, integrator)

    # DSS state to ensure we have continuous fields
    dss_state!(sim)
    return sim
end

"""
    eisenman_state_init(p_i::EisenmanIceParameters, space)

Initialize the prognostic state `Y` and the atmospheric forcing fields `Ya`
for the Eisenman-Zhang sea ice model. The initial state is ice-free, with the
mixed layer at 277 K.
"""
function eisenman_state_init(
    p_i::EisenmanIceParameters{FT},
    space::CC.Spaces.AbstractSpace,
) where {FT}
    Y = CC.Fields.FieldVector(
        T_sfc = ones(space) .* p_i.T_freeze,
        h_ice = CC.Fields.zeros(space),
        T_ml = ones(space) .* FT(277),
        e_base = CC.Fields.zeros(space),
    )
    Ya = (;
        F_turb = CC.Fields.zeros(space),
        SW_d = CC.Fields.zeros(space),
        LW_d = CC.Fields.zeros(space),
        ocean_qflux = CC.Fields.zeros(space),
    )
    return Y, Ya
end

"""
    solve_eisenman_model!(Y, Ya, p, Δt)

Advance the Eisenman-Zhang sea ice model state `Y` by one timestep `Δt`, given
the atmospheric forcing in `Ya` (`F_turb`, `SW_d`, `LW_d`, `ocean_qflux`) and
the parameter set `p = (; p_i, p_o)`.

The net upward atmospheric flux is assembled from the coupler-provided
downwelling radiation and the turbulent flux:

    F_atm = F_turb + ϵ (σ T_sfc^4 - LW_d) - (1 - α) SW_d

In ice-covered conditions the ice thickness responds to the imbalance between
`F_atm`, the basal flux `C0_base (T_ml - T_base)`, and the prescribed ocean
q-flux, and the surface temperature is obtained from one Newton iteration on
the surface energy balance F_atm = F_conductive. In ice-free conditions the
mixed layer absorbs `F_atm` directly and `T_sfc = T_ml`. The mixed layer is
not allowed to cool below freezing (the energy deficit forms frazil ice), and
a transition to ice-free conditions returns the residual melt energy to the
mixed layer.
"""
function solve_eisenman_model!(Y, Ya, p, Δt)
    # model parameter sets
    (; p_i, p_o) = p
    (; C0_base, T_base, L_ice, T_freeze, k_ice, α, ϵ, σ) = p_i
    hρc_ml = p_o.h * p_o.ρ * p_o.c

    # prognostic
    (; T_sfc, h_ice, T_ml) = Y

    # atmospheric forcing
    (; F_turb, SW_d, LW_d, ocean_qflux) = Ya

    FT = eltype(T_ml)
    Δt = FT(Δt)

    # net upward atmospheric flux and its derivative w.r.t. surface temperature.
    # NOTE: the turbulent flux contribution ∂F_turb/∂T_sfc is dropped from the
    # derivative (the coupler no longer provides it since #1284); the Newton
    # update below therefore treats the turbulent flux explicitly.
    F_atm = @. F_turb + ϵ * (σ * T_sfc^4 - LW_d) - (1 - α) * SW_d
    ∂F_atm∂T_sfc = @. 4 * ϵ * σ * T_sfc^3

    ice_covered = @. h_ice > 0

    # ice thickness and mixed layer temperature changes due to atmospheric and
    # oceanic fluxes (basal flux only acts under ice; the mixed layer only
    # exchanges with the atmosphere when ice-free)
    F_base = @. C0_base * (T_ml - T_base) * ice_covered
    ΔT_ml = @. ifelse(
        ice_covered,
        -(F_base - ocean_qflux) * Δt / hρc_ml,
        -(F_atm - ocean_qflux) * Δt / hρc_ml,
    )
    Δh_ice = @. ifelse(ice_covered, (F_atm - F_base - ocean_qflux) * Δt / L_ice, zero(F_atm))
    @. Y.e_base += F_base * Δt

    # T_ml is not allowed to be below freezing; the energy deficit forms frazil ice
    frazil_ice_formation = @. (T_ml + ΔT_ml) < T_freeze
    @. Δh_ice = ifelse(
        frazil_ice_formation,
        Δh_ice - (T_ml + ΔT_ml - T_freeze) * hρc_ml / L_ice,
        Δh_ice,
    )
    @. ΔT_ml = ifelse(frazil_ice_formation, T_freeze - T_ml, ΔT_ml)

    # adjust ocean temperature if transition to ice-free
    transition_to_icefree = @. ice_covered & ((h_ice + Δh_ice) <= 0)
    @. ΔT_ml = ifelse(
        transition_to_icefree,
        ΔT_ml - (h_ice + Δh_ice) * L_ice / hρc_ml,
        ΔT_ml,
    )
    @. Δh_ice = ifelse(transition_to_icefree, -h_ice, Δh_ice)

    # solve for T_sfc
    remains_ice_covered = @. (h_ice + Δh_ice) > 0
    # if ice-covered, solve implicitly (one Newton iteration: δT_sfc = -F(T_sfc) / F'(T_sfc));
    # the surface temperature is capped at freezing (NB: T_sfc not storing energy).
    # `h_safe` avoids dividing by zero where the surface is ice-free (the value is unused there).
    h_safe = @. ifelse(remains_ice_covered, h_ice + Δh_ice, one(h_ice))
    δT_sfc = @. (-F_atm + k_ice / h_safe * (T_base - T_sfc)) /
       (k_ice / h_safe + ∂F_atm∂T_sfc)
    @. Y.T_sfc = ifelse(
        remains_ice_covered,
        min(T_sfc + δT_sfc, T_freeze),
        T_ml + ΔT_ml,
    )

    @. Y.T_ml += ΔT_ml
    @. Y.h_ice += Δh_ice

    return Y
end

"""
    eisenman_ice_rhs!(dY, Y, p, t)

Rhs method in the form required by `ClimaTimeSteppers`. The Eisenman-Zhang
model is a discrete per-timestep update rather than a continuous tendency, so
the update is evaluated on a scratch copy of the state (leaving `Y` untouched
for the stepper's stage evaluations) and returned as the equivalent tendency
`dY = (Y⁺ - Y) / Δt`.
"""
function eisenman_ice_rhs!(dY, Y, p, t)
    Δt = float(p.Δt)
    Y⁺ = p.Y_scratch
    Y⁺ .= Y
    solve_eisenman_model!(Y⁺, p, p.params, Δt)
    @. dY = (Y⁺ - Y) / Δt
end

# binary ice mask, used to blend ice and ocean surface properties
ice_mask(sim::EisenmanIceSimulation) =
    @. ifelse(sim.integrator.u.h_ice > 0, one(sim.integrator.u.h_ice), zero(sim.integrator.u.h_ice))

# extensions required by Interfacer
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:area_fraction}) =
    sim.integrator.p.area_fraction
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:ice_concentration}) =
    sim.integrator.p.area_fraction
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:emissivity}) =
    sim.integrator.p.params.p_i.ϵ
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:roughness_buoyancy}) =
    @. sim.integrator.p.params.p_i.z0b * $(ice_mask(sim)) +
       sim.integrator.p.params.p_o.z0b * (1 - $(ice_mask(sim)))
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:roughness_momentum}) =
    @. sim.integrator.p.params.p_i.z0m * $(ice_mask(sim)) +
       sim.integrator.p.params.p_o.z0m * (1 - $(ice_mask(sim)))
Interfacer.get_field(
    sim::EisenmanIceSimulation,
    ::Union{Val{:surface_direct_albedo}, Val{:surface_diffuse_albedo}},
) =
    @. sim.integrator.p.params.p_i.α * $(ice_mask(sim)) +
       sim.integrator.p.params.p_o.α * (1 - $(ice_mask(sim)))
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:roughness_model}) = :constant
Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:surface_temperature}) =
    sim.integrator.u.T_sfc

"""
    Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:energy})

Extension of Interfacer.get_field to get the energy of the surface column.
It is the sum of the heat content of the mixed layer, the latent heat of the
ice, the accumulated prescribed ocean q-flux, and the accumulated basal flux.
"""
function Interfacer.get_field(sim::EisenmanIceSimulation, ::Val{:energy})
    p_i = sim.integrator.p.params.p_i
    p_o = sim.integrator.p.params.p_o
    u = sim.integrator.u
    FT = eltype(u)
    hρc_ml = p_o.h * p_o.ρ * p_o.c
    t = FT(float(sim.integrator.t))
    return @. hρc_ml * u.T_ml +
       p_i.L_ice * u.h_ice +
       sim.integrator.p.ocean_qflux * t +
       u.e_base
end

function Interfacer.update_field!(
    sim::EisenmanIceSimulation,
    ::Val{:area_fraction},
    field::CC.Fields.Field,
)
    sim.integrator.p.area_fraction .= field
    return nothing
end
function Interfacer.update_field!(sim::EisenmanIceSimulation, ::Val{:SW_d}, field)
    Interfacer.remap!(sim.integrator.p.SW_d, field)
end
function Interfacer.update_field!(sim::EisenmanIceSimulation, ::Val{:LW_d}, field)
    Interfacer.remap!(sim.integrator.p.LW_d, field)
end
function Interfacer.update_field!(
    sim::EisenmanIceSimulation,
    ::Val{:turbulent_energy_flux},
    field,
)
    Interfacer.remap!(sim.integrator.p.F_turb, field)
end
Interfacer.update_field!(
    sim::EisenmanIceSimulation,
    ::Val{:turbulent_moisture_flux},
    field,
) = nothing

function FluxCalculator.update_turbulent_fluxes!(
    sim::EisenmanIceSimulation,
    fields::NamedTuple,
)
    Interfacer.update_field!(sim, Val(:turbulent_energy_flux), fields.F_lh .+ fields.F_sh)
    return nothing
end

"""
    Checkpointer.get_model_prog_state(sim::EisenmanIceSimulation)

Extension of Checkpointer.get_model_prog_state to get the model state.
"""
function Checkpointer.get_model_prog_state(sim::EisenmanIceSimulation)
    return sim.integrator.u
end

function Checkpointer.get_model_cache(sim::EisenmanIceSimulation)
    return sim.integrator.p
end

function Checkpointer.restore_cache!(sim::EisenmanIceSimulation, new_cache)
    old_cache = Checkpointer.get_model_cache(sim)
    for p in propertynames(old_cache)
        if getproperty(old_cache, p) isa CC.Fields.Field
            ArrayType = ClimaComms.array_type(getproperty(old_cache, p))
            parent(getproperty(old_cache, p)) .=
                ArrayType(parent(getproperty(new_cache, p)))
        end
    end
end

"""
    dss_state!(sim::EisenmanIceSimulation)

Perform DSS on the state of a component simulation, intended to be used
before the initial step of a run. This method acts on Eisenman-Zhang ice simulations.
"""
dss_state!(sim::EisenmanIceSimulation) =
    Utilities.apply_dss!(sim.integrator.u, sim.integrator.p.dss_buffer)
