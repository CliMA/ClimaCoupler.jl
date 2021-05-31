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
"""
    Declaration of state variables

    vars_state returns a NamedTuple of data types.
"""
function vars_state(model::SlabLandModelSetup, aux::Auxiliary, T)
    # orientation = model.physics.orientation

    @vars begin
        x::T
        y::T
        z::T
        Φ::T
        ∇Φ::SVector{3, T} # TODO: only needed for the linear model
        # orientation::vars_state(orientation, aux, T)
        F_ρθ_prescribed::T # stores prescribed flux for secondary (ocean import)
    end
end

function vars_state(::SlabLandModelSetup, ::Prognostic, T)
    @vars begin
        T_sfc::T
    end
end

vars_state(::SlabLandModelSetup, ::Gradient, T) = @vars()

vars_state(::SlabLandModelSetup, ::GradientFlux, T) = @vars()


"""
    Initialization of state variables

    init_state_xyz! sets up the initial fields within our state variables
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
    # init_state_auxiliary!(m, m.physics.ref_state, state_auxiliary, geom)

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


# function init_state_auxiliary!(
#     model::SlabLandModelSetup,
#     state_auxiliary::MPIStateArray,
#     grid,
#     direction,
# )
#     # domain = model.numerics.grid.domain
#     orientation = model.physics.orientation

#     # helper function for storing coordinates (used only once here)
#     function init_coords!(::SlabLandModelSetup, aux, geom)
#         aux.x = geom.coord[1]
#         aux.y = geom.coord[2]
#         aux.z = geom.coord[3]
    
#         return nothing
#     end

#     # store coordinates
#     init_state_auxiliary!(
#         model,
#         (model, aux, tmp, geom) -> init_coords!(model, aux, geom),
#         state_auxiliary,
#         grid,
#         direction,
#     )

#     # store orientation
#     init_state_auxiliary!(
#         model,
#         # (model, aux, tmp, geom) -> orientation_nodal_init_aux!(orientation, domain, aux, geom),
#         (model, aux, tmp, geom) -> orientation_nodal_init_aux!(orientation, aux, geom),
#         state_auxiliary,
#         grid,
#         direction,
#     )

#     # store vertical unit vector
#     orientation_gradient(model, orientation, state_auxiliary, grid, direction)
# end

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


"""
function preA(csolver)
    - saves and regrids the EnergyFluxB couplerfield (i.e. regridded mB.state.F_ρθ_accum[mB.boundary]) to ρθ_prescribed[mA.boundary] on the domainA grid
    csolver::CplSolver
"""
function preLand(csolver)
    mA = csolver.component_list.domainLand.component_model
    mB = csolver.component_list.domainAtmos.component_model
    # Set mean air-sea theta flux

    mA.odesolver.rhs!.state_auxiliary.F_ρθ_prescribed[mA.boundary] .= 
        coupler_get(csolver.coupler, :EnergyFluxAtmos, mA.grid, DateTime(0), u"J")
    # Set ocean boundary flux accumulator to 0. (this isn't used)

    # @info(
    #     "preocean",
    #     time = csolver.t,
    #     F_prescribed_max =
    #         maximum(mA.discretization.state_auxiliary.F_prescribed[mA.boundary]),
    #     F_prescribed_min =
    #         maximum(mA.discretization.state_auxiliary.F_prescribed[mA.boundary]),
    #     ocean_θ_surface_max = maximum(mA.state.θ[mA.boundary]),
    #     ocean_θ_surface_min = maximum(mA.state.θ[mA.boundary]),
    # )
end

"""
function postA(csolver)
    - updates couplerfield EnergyA with mA.state.ρθ[mA.boundary] regridded to the coupler grid, and updates the coupler time
    - updates EnergyA with mA.state.ρθ[mA.boundary], and the coupler time
    csolver::CplSolver
"""
function postLand(csolver)
    mA = csolver.component_list.domainLand.component_model
    mB = csolver.component_list.domainAtmos.component_model
    # @info(
    #     "postocean",
    #     time = csolver.t + csolver.dt,
    #     ocean_θ_surface_max = maximum(mA.state.θ[mA.boundary]),
    #     ocean_θ_surface_min = maximum(mA.state.θ[mA.boundary]),
    # )

    # Pass ocean exports to "coupler" namespace
    #  1. Ocean SST (value of θ at z=0)
    coupler_put!(csolver.coupler, :EnergyLand, mA.state.T_sfc[mA.boundary], mA.grid.numerical, DateTime(0), u"J")
end

#FluidBC{Impenetrable{FreeSlip},Insulating}

# function numerical_boundary_flux_second_order!(
#     numerical_flux::Union{PenaltyNumFluxDiffusive},
#     bctype,
#     balance_law::SlabLandModelSetup,
#     args... 
# ) where {S, D, A, HD}
#     #@show "hellooo!"

#     numerical_boundary_flux_second_order!( # this overwrites the total flux for T_sfc (using methods below)
#         numerical_flux,
#         bctype.T_sfc,
#         balance_law,
#         args...
#     )

# end

# customized methods form specifying total normal fluxes
function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive},
    bctype::NamedTuple{(:T_sfc,),Tuple{CoupledSecondaryBoundary}}, 
    balance_law::SlabLandModelSetup,
    fluxᵀn::Vars{S},
    normal_vector::SVector,
    state_prognostic⁻::Vars{S},
    state_gradient_flux⁻::Vars{D},
    state_hyperdiffusive⁻::Vars{HD},
    state_auxiliary⁻::Vars{A},
    state_prognostic⁺::Vars{S},
    state_gradient_flux⁺::Vars{D},
    state_hyperdiffusive⁺::Vars{HD},
    state_auxiliary⁺::Vars{A},
    t,
    state1⁻::Vars{S},
    diff1⁻::Vars{D},
    aux1⁻::Vars{A},
) where {S, D, A, HD}

    p = balance_law.parameters
    G    = p.κ_s * (state_prognostic⁻.T_sfc - p.T_h) / p.h_s # simple soil physics
    fluxᵀn.T_sfc = - state_auxiliary⁺.F_ρθ_prescribed + G

end
function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive}, bctype::NamedTuple{(:T_sfc,),Tuple{Insulating}}, balance_law::SlabLandModelSetup, fluxᵀn::Vars{S}, _...,) where {S}
    FT = eltype(fluxᵀn)
    fluxᵀn.T_sfc = FT(0)
end


# # customized flux for Insulating boundary
# numerical_boundary_flux_second_order!(
#     numerical_flux::Union{PenaltyNumFluxDiffusive},
#     bctype::FluidBC{Impenetrable{FreeSlip},Insulating},
#     balance_law::SlabLandModelSetup,
#     fluxᵀn::Vars{S},
#     normal_vector::SVector,
#     state_prognostic⁻::Vars{S},
#     state_gradient_flux⁻::Vars{D},
#     state_hyperdiffusive⁻::Vars{HD},
#     state_auxiliary⁻::Vars{A},
#     state_prognostic⁺::Vars{S},
#     state_gradient_flux⁺::Vars{D},
#     state_hyperdiffusive⁺::Vars{HD},
#     state_auxiliary⁺::Vars{A},
#     t,
#     state1⁻::Vars{S},
#     diff1⁻::Vars{D},
#     aux1⁻::Vars{A},
# ) where {S, D, HD, A} = normal_boundary_flux_second_order!(
#     numerical_flux,
#     bctype,
#     balance_law,
#     fluxᵀn,
#     normal_vector,
#     state_prognostic⁻,
#     state_gradient_flux⁻,
#     state_hyperdiffusive⁻,
#     state_auxiliary⁻,
#     state_prognostic⁺,
#     state_gradient_flux⁺,
#     state_hyperdiffusive⁺,
#     state_auxiliary⁺,
#     t,
#     state1⁻,
#     diff1⁻,
#     aux1⁻,
# )

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