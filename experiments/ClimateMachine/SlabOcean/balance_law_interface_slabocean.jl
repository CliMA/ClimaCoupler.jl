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

using Unitful

##
#     Declaration of variables
##
"""
    vars_state 
    - returns a NamedTuple of data types.
"""
function vars_state(model::SlabOceanModelSetup, aux::Auxiliary, T)

    @vars begin
        x::T
        y::T
        z::T
        F_ρe_prescribed::T # stores prescribed flux for secondary (land import)
    end
end

function vars_state(::SlabOceanModelSetup, ::Prognostic, T)
    @vars begin
        T_sfc::T
    end
end

vars_state(::SlabOceanModelSetup, ::Gradient, T) = @vars()

vars_state(::SlabOceanModelSetup, ::GradientFlux, T) = @vars()


##
#     Initialization of variables
##
"""
    init_state_xyz! 
    - sets up the initial fields within our state variables
    (e.g., prognostic, auxiliary, etc.), however it seems to not initialized
    the gradient flux variables by default.
"""

function nodal_init_state_auxiliary!(m::SlabOceanModelSetup, state_auxiliary, tmp, geom)
    init_state_auxiliary!(m, m.physics.orientation, state_auxiliary, geom)

    state_auxiliary.F_ρe_prescribed = 0
end

function init_state_auxiliary!(::SlabOceanModelSetup, ::SphericalOrientation, state_auxiliary, geom)
    FT = eltype(state_auxiliary)
    r = norm(geom.coord)
    state_auxiliary.x = geom.coord[1]
    state_auxiliary.y = geom.coord[2]
    state_auxiliary.z = geom.coord[3]

end

function init_state_auxiliary!(::SlabOceanModelSetup, ::Union{NoOrientation, FlatOrientation}, state_auxiliary, geom)
    FT = eltype(state_auxiliary)
    @inbounds r = geom.coord[3]
    state_auxiliary.x = geom.coord[1]
    state_auxiliary.y = geom.coord[2]
    state_auxiliary.z = geom.coord[3]
end

function init_state_prognostic!(model::SlabOceanModelSetup, state::Vars, aux::Vars, localgeo, t)
    x = aux.x
    y = aux.y
    z = aux.z

    parameters = model.physics.parameters
    ic = model.initial_conditions

    FT = eltype(state)
    state.T_sfc = FT(0)
    state.T_sfc = ic.T_sfc(parameters, x, y, z)

    return nothing
end

"""
    LHS computations
"""
flux_first_order!(model::SlabOceanModelSetup, _...) = nothing
compute_gradient_argument!(model::SlabOceanModelSetup, _...) = nothing
compute_gradient_flux!(model::SlabOceanModelSetup, _...) = nothing
flux_second_order!(model::SlabOceanModelSetup, _...) = nothing

"""
    RHS computations
"""
@inline function source!(
    model::SlabOceanModelSetup,
    source::Vars,
    state::Vars,
    gradflux::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    x = aux.x
    y = aux.y
    z = aux.z
    p = model.physics.parameters
    FT = eltype(aux)

    Q = FT(0)
    if boundary_mask(p, x, y, z)
        Q = idealized_qflux(p, x, y, z)
    end

    source.T_sfc = (aux.F_ρe_prescribed + Q) / (p.ρ_o * p.c_o * p.h_o)

    return nothing
end

@inline function idealized_qflux(p, x...)
    """
    Calculates additional tropical fluxes for a more realistic representation of surface tropical heating in idealized experiments
    """
    FT = eltype(p.Q_0)
    ϕ = lat(x...)
    Q = p.Q_0 * (FT(1) - 2 * ϕ^2 / p.L_w^2) * exp(-ϕ^2 / p.L_w^2) / cos(ϕ)

    ϕ_Qcutoff = FT(85 / 180 * pi) # this is not in the GFLD model but here needed to avoid polar discontinuities
    Q = abs(ϕ) > ϕ_Qcutoff ? FT(0) : Q
end

"""
    Numerics and boundary conditions
    - no interior / exterior interface fluxes needed for the slab
"""
@inline boundary_conditions(model::SlabOceanModelSetup) = model.boundary_conditions

numerical_flux_first_order!(::Nothing, model::Union{SlabOceanModelSetup}, _...) = nothing
numerical_flux_gradient!(nf, bc, model::Union{SlabOceanModelSetup}, _...) = nothing
function numerical_flux_second_order!(::Nothing, model::Union{SlabOceanModelSetup}, fluxᵀn, _...)
    fluxᵀn.T_sfc = fluxᵀn.T_sfc .* Float64(0) # this overwrites the internal numerical flux contributions to the normal flux
end

numerical_boundary_flux_first_order!(::Nothing, bc, model::Union{SlabOceanModelSetup}, _...) = nothing
numerical_boundary_flux_gradient!(::Nothing, bc, model::Union{SlabOceanModelSetup}, _...) = nothing
function numerical_boundary_flux_second_order!(::Nothing, bc, model::Union{SlabOceanModelSetup}, fluxᵀn, _...)
    fluxᵀn.T_sfc = fluxᵀn.T_sfc .* Float64(0) # this overwrites the internal numerical flux contributions to the normal flux
end
