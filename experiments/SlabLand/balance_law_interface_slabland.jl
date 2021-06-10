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
    - sets up the initial fields within our state variables
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
    #G    = p.κ_s * (state.T_sfc - p.T_h) / p.h_s # simple soil physics
    G = Float64(0)
    source.T_sfc = - (aux.F_ρθ_prescribed + G) / (p.ρ_s * p.c_s * p.h_s)

    return nothing
end

"""
    Boundary conditions
"""
@inline boundary_conditions(model::SlabLandModelSetup) = model.boundary_conditions

@inline function boundary_state!(
        numerical_flux,
        bc,
        model::SlabLandModelSetup,
        args...
    )
    # We need to apply boundary conditions for state variables. 
    # This depends on first, second, and high-order
    # operations, hence the diffusion model dependence.
    # TODO!: make work for higher-order diffusion
    diffusion = model.physics.diffusion # defaults to `nothing`
    calc_boundary_state!(numerical_flux, bc.T_sfc, model, diffusion, args...)
end

@inline vertical_unit_vector(::Orientation, aux) = aux.∇Φ / grav(param_set)
@inline vertical_unit_vector(::NoOrientation, aux) = @SVector [0, 0, 1]

# customized methods form specifying total normal fluxes
numerical_boundary_flux_second_order!(numerical_flux::Union{PenaltyNumFluxDiffusive}, bctype::NamedTuple{(:T_sfc,),Tuple{CoupledSecondaryBoundary}}, balance_law::SlabLandModelSetup, _...,) = nothing

function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive}, bctype::NamedTuple{(:T_sfc,),Tuple{Insulating}}, balance_law::SlabLandModelSetup, fluxᵀn::Vars{S}, _...,) where {S}
    FT = eltype(fluxᵀn)
    fluxᵀn.T_sfc = FT(0)
end

"""
function numerical_flux_second_order!(::PenaltyNumFluxDiffusive, 
  Penalty flux formulation of second order numerical flux. This formulation
  computes the CentralNumericalFluxSecondOrder term first (which is just the average
  of the + and - fluxes and an edge), and then adds a "penalty" flux that relaxes
  the edge state + and - toward each other.
"""
function numerical_flux_second_order!(
    ::PenaltyNumFluxDiffusive,
    bl::SlabLandModelSetup,
    fluxᵀn::Vars{S},
    n::SVector,
    state⁻::Vars{S},
    diff⁻::Vars{D},
    hyperdiff⁻::Vars{HD},
    aux⁻::Vars{A},
    state⁺::Vars{S},
    diff⁺::Vars{D},
    hyperdiff⁺::Vars{HD},
    aux⁺::Vars{A},
    t,
) where {S, HD, D, A}
    
    numerical_flux_second_order!(
        CentralNumericalFluxSecondOrder(),
        bl,
        fluxᵀn,
        n,
        state⁻,
        diff⁻,
        hyperdiff⁻,
        aux⁻,
        state⁺,
        diff⁺,
        hyperdiff⁺,
        aux⁺,
        t,
    )

    Fᵀn = parent(fluxᵀn)
    FT = eltype(Fᵀn)
    tau = FT(0.0) # decide if want to use penalty
    Fᵀn .+= tau * (parent(state⁻) - parent(state⁺))
end

numerical_flux_first_order!(::RoeNumericalFlux, model::Union{SlabLandModelSetup}, _...,) = nothing

"""
function preA(csolver)
    - saves and regrids the `EnergyFluxAtmos` couplerfield (i.e. regridded `mAtmos.state.F_ρθ_accum[mAtmos.boundary]`) to `F_ρθ_prescribed[mLand.boundary]` on the `mLand` grid
    csolver::CplSolver
"""
function preLand(csolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model

    mLand.odesolver.rhs!.state_auxiliary.F_ρθ_prescribed[mLand.boundary] .= 
        coupler_get(csolver.coupler, :EnergyFluxAtmos, mLand.grid, DateTime(0), u"J")
end

"""
function postA(csolver)
    - updates couplerfield `LandSurfaceTemerature` with `mLand.state.T_sfc[mLand.boundary]` regridded to the coupler grid, and updates the coupler time
    csolver::CplSolver
"""
function postLand(csolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    coupler_put!(csolver.coupler, :LandSurfaceTemerature, mLand.state.T_sfc[mLand.boundary], mLand.grid.numerical, DateTime(0), u"K")
end