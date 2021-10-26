import ClimateMachine.BalanceLaws:
# declaration
    vars_state,
# initialization
    init_state_prognostic!,
    init_state_auxiliary!,
    nodal_init_state_auxiliary!,
# rhs computation
    compute_gradient_argument!,
    compute_gradient_flux!,
    flux_first_order!,
    flux_second_order!,
    source!,
# boundary conditions
    boundary_conditions,
    boundary_state!

import ClimateMachine.DGMethods.NumericalFluxes:
    CentralNumericalFluxGradient,
    CentralNumericalFluxSecondOrder,
    NumericalFluxFirstOrder,
    NumericalFluxSecondOrder,
    RusanovNumericalFlux,
    numerical_boundary_flux_second_order!,
    numerical_flux_second_order!,
    numerical_boundary_flux_first_order!,
    numerical_flux_first_order!,
    normal_boundary_flux_second_order!

##
#     Declaration of variables
##
"""
    vars_state 
    - returns a NamedTuple of data types.
"""
function vars_state(model::ModelSetup, aux::Auxiliary, T)

    @vars begin
        x::T
        y::T
        z::T
        T_sfc::T  # land surface temperature
    end
end

function vars_state(::ModelSetup, ::Prognostic, T)
    @vars begin
        ρ::T
        ρu::SVector{3, T}
        ρθ::T
        F_ρθ_accum::T
    end
end

function vars_state(::ModelSetup, ::Gradient, T)
    @vars begin
        ∇ρ::T
        ∇u::SVector{3, T}
        ∇θ::T
    end
end

function vars_state(::ModelSetup, ::GradientFlux, T)
    @vars begin
        μ∇ρ::SVector{3, T}
        ν∇u::SMatrix{3, 3, T, 9}
        κ∇θ::SVector{3, T}
    end
end

##
#     Initialization of variables
##
"""
    init_state_xyz! 
    - sets up the initial fields within our state variables
    (e.g., prognostic, auxiliary, etc.), however it seems to not initialized
    the gradient flux variables by default.
"""

function nodal_init_state_auxiliary!(m::ModelSetup, state_auxiliary, tmp, geom)
    init_state_auxiliary!(m, m.physics.orientation, state_auxiliary, geom)

    state_auxiliary.T_sfc = 0
end

function init_state_auxiliary!(::ModelSetup, ::Union{NoOrientation, FlatOrientation}, state_auxiliary, geom)
    FT = eltype(state_auxiliary)
    _grav = FT(grav(param_set))
    @inbounds r = geom.coord[3]
    state_auxiliary.x = geom.coord[1]
    state_auxiliary.y = geom.coord[2]
    state_auxiliary.z = geom.coord[3]
end

function init_state_prognostic!(model::ModelSetup, state::Vars, aux::Vars, localgeo, t)
    x = aux.x
    y = aux.y
    z = aux.z

    parameters = model.parameters
    ic = model.initial_conditions

    state.ρ = ic.ρ(parameters, x, y, z)
    state.ρu = ic.ρu(parameters, x, y, z)
    state.ρθ = ic.ρθ(parameters, x, y, z)

    state.F_ρθ_accum = 0

    return nothing
end

"""
    LHS computations
"""
@inline function flux_first_order!(model::ModelSetup, flux::Grad, state::Vars, aux::Vars, t::Real, direction)
    flux.ρu += calc_pressure(model.physics.eos, state) * I
    calc_advective_flux!(flux, model.physics.advection, state, aux, t)

    return nothing
end

@inline function compute_gradient_argument!(model::ModelSetup, grad::Vars, state::Vars, aux::Vars, t::Real)
    calc_diffusive_flux_argument!(grad, model.physics.diffusion, state, aux, t)

    return nothing
end

@inline function compute_gradient_flux!(model::ModelSetup, gradflux::Vars, grad::Grad, state::Vars, aux::Vars, t::Real)
    calc_diffusive_flux!(gradflux, model.physics.diffusion, grad, state, aux, t)

    return nothing
end

@inline function flux_second_order!(
    model::ModelSetup,
    flux::Grad,
    state::Vars,
    gradflux::Vars,
    ::Vars,
    aux::Vars,
    t::Real,
)
    flux.ρ += gradflux.μ∇ρ
    flux.ρu += gradflux.ν∇u
    flux.ρθ += gradflux.κ∇θ

    return nothing
end

"""
    RHS computations
"""
@inline function source!(model::ModelSetup, source::Vars, state::Vars, gradflux::Vars, aux::Vars, t::Real, direction)
    source.F_ρθ_accum = calculate_land_sfc_fluxes(model, state, aux, t)

    return nothing
end
