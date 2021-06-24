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

Returns a NamedTuple of data types.
"""
function vars_state(model::SlabLandModelSetup, aux::Auxiliary, T)

    @vars begin
        x::T
        y::T
        z::T
        F_ρθ_prescribed::T # stores prescribed flux for secondary (land import)
    end
end

function vars_state(::SlabLandModelSetup, ::Prognostic, T)
    @vars begin
        T_sfc::T
    end
end

vars_state(::SlabLandModelSetup, ::Gradient, T) = @vars()

vars_state(::SlabLandModelSetup, ::GradientFlux, T) = @vars()


##
#     Initialization of variables
##
"""
    init_state_xyz!

Sets up the initial fields within our state variables
(e.g., prognostic, auxiliary, etc.), however it seems to not initialized
the gradient flux variables by default.
"""

function nodal_init_state_auxiliary!(
    m::SlabLandModelSetup,
    state_auxiliary,
    tmp,
    geom,
)
    init_state_auxiliary!(m, m.physics.orientation, state_auxiliary, geom)

    state_auxiliary.F_ρθ_prescribed = 0
end

function init_state_auxiliary!(
    ::SlabLandModelSetup,
    ::SphericalOrientation,
    state_auxiliary,
    geom,
)
    FT = eltype(state_auxiliary)
    _grav = FT(grav(param_set))
    r = norm(geom.coord)
    state_auxiliary.x = geom.coord[1]
    state_auxiliary.y = geom.coord[2]
    state_auxiliary.z = geom.coord[3]
    state_auxiliary.Φ = _grav * r
    state_auxiliary.∇Φ = _grav * geom.coord / r
end

function init_state_auxiliary!(
    ::SlabLandModelSetup,
    ::Union{NoOrientation,FlatOrientation},
    state_auxiliary,
    geom,
)
    FT = eltype(state_auxiliary)
    _grav = FT(grav(param_set))
    @inbounds r = geom.coord[3]
    state_auxiliary.x = geom.coord[1]
    state_auxiliary.y = geom.coord[2]
    state_auxiliary.z = geom.coord[3]
end

function init_state_prognostic!(model::SlabLandModelSetup, state::Vars, aux::Vars, localgeo, t)
    x = aux.x
    y = aux.y
    z = aux.z

    parameters = model.parameters
    ic = model.initial_conditions

    state.T_sfc = ic.T_sfc(parameters, x, y, z)

    return nothing
end

"""
    LHS computations
"""
flux_first_order!(model::SlabLandModelSetup, _...,) = nothing
compute_gradient_argument!(model::SlabLandModelSetup, _...,) = nothing
compute_gradient_flux!(model::SlabLandModelSetup, _...,) = nothing
flux_second_order!(model::SlabLandModelSetup, _...,) = nothing

"""
    RHS computations
"""
@inline function source!(
    model::SlabLandModelSetup,
    source::Vars,
    state::Vars,
    gradflux::Vars,
    aux::Vars,
    t::Real,
    direction,
)

    p = model.parameters
    G    = p.κ_s * (state.T_sfc - p.T_h) / p.h_s # simple soil physics
    source.T_sfc = - (aux.F_ρθ_prescribed + G) / (p.ρ_s * p.c_s * p.h_s)

    return nothing
end

"""
    Numerics and boundary conditions
    - no interior / exterior interface fluxes needed for the slab
"""
@inline boundary_conditions(model::SlabLandModelSetup) = model.boundary_conditions

numerical_flux_first_order!(::Nothing, model::Union{SlabLandModelSetup}, _...,) = nothing
numerical_flux_gradient!(nf, bc,  model::Union{SlabLandModelSetup}, _...,) = nothing
numerical_flux_second_order!(::Nothing, model::Union{SlabLandModelSetup}, _...,) = nothing

numerical_boundary_flux_first_order!(::Nothing, bc, model::Union{SlabLandModelSetup}, _...,) = nothing
numerical_boundary_flux_gradient!(::Nothing, bc, model::Union{SlabLandModelSetup}, _...,) = nothing
numerical_boundary_flux_second_order!(::Nothing, bc, model::Union{SlabLandModelSetup}, _...,) = nothing