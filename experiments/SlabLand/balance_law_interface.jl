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
function vars_state(model::ModelSetup, aux::Auxiliary, T)
    # orientation = model.physics.orientation

    @vars begin
        x::T
        y::T
        z::T
        Φ::T
        ∇Φ::SVector{3, T} # TODO: only needed for the linear model
        # orientation::vars_state(orientation, aux, T)
        T_sfc::T  # stores opposite face for primary (atmospheric import)
        F_ρθ_prescribed::T # stores prescribed flux for secondary (ocean import)
    end
end

function vars_state(::ModelSetup, ::Prognostic, T)
    @vars begin
        ρ::T
        ρu::SVector{3,T}
        ρθ::T
        F_ρθ_accum::T
    end
end

function vars_state(::ModelSetup, ::Gradient, T)
    @vars begin
        ∇ρ::T
        ∇u::SVector{3,T}
        ∇θ::T
    end
end

function vars_state(::ModelSetup, ::GradientFlux, T)
    @vars begin
        μ∇ρ::SVector{3,T}
        ν∇u::SMatrix{3,3,T,9}
        κ∇θ::SVector{3,T}
    end
end

"""
    Initialization of state variables

    init_state_xyz! sets up the initial fields within our state variables
    (e.g., prognostic, auxiliary, etc.), however it seems to not initialized
    the gradient flux variables by default.
"""

function nodal_init_state_auxiliary!(
    m::ModelSetup,
    state_auxiliary,
    tmp,
    geom,
)
    init_state_auxiliary!(m, m.physics.orientation, state_auxiliary, geom)
    # init_state_auxiliary!(m, m.physics.ref_state, state_auxiliary, geom)

    state_auxiliary.T_sfc = 0
    state_auxiliary.F_ρθ_prescribed = 0
end

function init_state_auxiliary!(
    ::ModelSetup,
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
    ::ModelSetup,
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
    state_auxiliary.Φ = _grav * r
    state_auxiliary.∇Φ = SVector{3, FT}(0, 0, _grav)
end


# function init_state_auxiliary!(
#     model::ModelSetup,
#     state_auxiliary::MPIStateArray,
#     grid,
#     direction,
# )
#     # domain = model.numerics.grid.domain
#     orientation = model.physics.orientation

#     # helper function for storing coordinates (used only once here)
#     function init_coords!(::ModelSetup, aux, geom)
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
@inline function flux_first_order!(
    model::ModelSetup,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    flux.ρu += calc_pressure(model.physics.eos, state) * I

    calc_advective_flux!(flux, model.physics.advection, state, aux, t)

    return nothing
end

@inline function compute_gradient_argument!(
    model::ModelSetup,
    grad::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    calc_diffusive_flux_argument!(grad, model.physics.diffusion, state, aux, t)

    return nothing
end

@inline function compute_gradient_flux!(
    model::ModelSetup,
    gradflux::Vars,
    grad::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
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
@inline function source!(
    model::ModelSetup,
    source::Vars,
    state::Vars,
    gradflux::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    coriolis = model.physics.coriolis
    gravity = model.physics.gravity
    orientation = model.physics.orientation

    calc_force!(source, coriolis, state, aux, orientation, t)
    calc_force!(source, gravity, state, aux, orientation, t)

    source.F_ρθ_accum = calculate_land_sfc_fluxes(model, state, aux) / model.parameters.c_p
    
    return nothing
end

"""
    calculate_land_sfc_fluxes(model::ModelSetup, state, aux)
- calculate furface fluxes using the bulk gradient diffusion theory
"""
function calculate_land_sfc_fluxes(model::ModelSetup, state, aux; MO_params = nothing) # should pass in the coupler state (also move to coupler), so can access states of both models derectly -e.g. callback?
    
    FT = eltype(state)
    p = model.parameters

    p_sfc = FT(100000)
    θ_a = state.ρθ / state.ρ 
    T_a = θ_a
    T_sfc = aux.T_sfc 
    θ_sfc = T_sfc  

    R_SW = (FT(1)-p.α) * p.τ * p.F_sol
    R_LW = p.ϵ * (p.σ * θ_sfc - p.F_a)
    SH   = p.c_p * p.g_a * (T_sfc - T_a) # g_a could be substituded by g_a=f(MO_params) (see Bonan P90), but first need to modify SurfaceFluxes.jl
      
    #LH   = p.lambda * p.g_w * (q_sat(T_sfc, p_sfc) - q_a) 

    F_tot = - (R_SW - R_LW - SH )#- LH 
end


"""
    Boundary conditions
"""
@inline boundary_conditions(model::ModelSetup) = model.boundary_conditions

@inline function boundary_state!(
        numerical_flux,
        bc::FluidBC,
        model::ModelSetup,
        args...
    )
    # We need to apply boundary conditions for state variables. 
    # This depends on first, second, and high-order
    # operations, hence the diffusion model dependence.
    # TODO!: make work for higher-order diffusion
    diffusion = model.physics.diffusion # defaults to `nothing`
    calc_boundary_state!(numerical_flux, bc.ρu, model, diffusion, args...)
    calc_boundary_state!(numerical_flux, bc.ρθ, model, diffusion, args...)
end

@inline vertical_unit_vector(::Orientation, aux) = aux.∇Φ / grav(param_set)
@inline vertical_unit_vector(::NoOrientation, aux) = @SVector [0, 0, 1]


"""
function preB(csolver)
    - saves and regrids the EnergyA couplerfield (i.e. regridded mA.state.ρθ[mA.boundary]) on coupler grid to mB.state.ρθ_secondary[mB.boundary] on the domainB grid
    csolver::CplSolver
"""
function preAtmos(csolver)
    mA = csolver.component_list.domainLand.component_model
    mB = csolver.component_list.domainAtmos.component_model
    # Set boundary SST used in atmos to SST of ocean surface at start of coupling cycle.
    mB.odesolver.rhs!.state_auxiliary.T_sfc[mB.boundary] .= 
        coupler_get(csolver.coupler, :EnergyLand, mB.grid.numerical, DateTime(0), u"J")
    # Set atmos boundary flux accumulator to 0.
    mB.state.F_ρθ_accum .= 0

    idx = varsindex(vars(mB.state), :ρθ)[1]
    idx_rho = varsindex(vars(mB.state), :ρ)[1]
    @info(
        "preatmos",
        time = csolver.t, #* "/" * mB.time.finish ,
        total_ρθA_ = weightedsum(mA.state, 1),
        total_ρθB = weightedsum(mB.state, idx) / weightedsum(mB.state, idx_rho),
        total_ρθ = weightedsum(mA.state, 1) + weightedsum(mB.state, idx) / weightedsum(mB.state, idx_rho),
        atmos_ρθ_surface_maxA = maximum(mA.state.T_sfc[mA.boundary]),
        ocean_ρθ_surface_maxB = maximum(mB.state.ρθ[mB.boundary]),
    )

    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.A[csolver.steps] = weightedsum(mA.state, 1)
    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.B[csolver.steps] = weightedsum(mB.state, idx) / weightedsum(mB.state, idx_rho)
end

"""
function postB(csolver)
    - updates couplerfield EnergyFluxB with mB.state.F_ρθ_accum[mB.boundary] regridded to the coupler grid, and updates the coupler time
    csolver::CplSolver
"""
function postAtmos(csolver)
    mA = csolver.component_list.domainLand.component_model
    mB = csolver.component_list.domainAtmos.component_model
    # Pass atmos exports to "coupler" namespace
    # 1. Save mean θ flux at the Atmos boundary during the coupling period
    coupler_put!(csolver.coupler, :EnergyFluxAtmos, mB.state.F_ρθ_accum[mB.boundary] ./ csolver.dt,
        mB.grid.numerical, DateTime(0), u"J")

    # @info(
    #     "postatmos",
    #     time = time = csolver.t + csolver.dt,
    #     total_θ_atmos = weightedsum(mB.state, 1),
    #     total_θ_ocean = weightedsum(mA.state, 1),
    #     total_F_accum = mean(mB.state.F_accum[mB.boundary]) * 1e6 * 1e6,
    #     total_θ =
    #         weightedsum(mB.state, 1) +
    #         weightedsum(mA.state, 1) +
    #         mean(mB.state.F_ρθ_accum[mB.boundary]) * 1e6 * 1e6,
    #     F_accum_max = maximum(mB.state.F_accum[mB.boundary]),
    #     F_avg_max = maximum(mB.state.F_accum[mB.boundary] ./ csolver.dt),
    #     atmos_θ_surface_max = maximum(mB.state.θ[mB.boundary]),
    #     ocean_θ_surface_max = maximum(mA.state.θ[mA.boundary]),
    # )
end


#FluidBC{Impenetrable{FreeSlip},Insulating}

function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive},
    bctype::FluidBC,
    balance_law::ModelSetup,
    args... 
) where {S, D, A, HD}
    #@show "hellooo!"
    normal_boundary_flux_second_order!( # for ρu BCs (to be implemented in boundary_state!)
        numerical_flux,
        bctype,
        balance_law,
        args...
    )

    numerical_boundary_flux_second_order!( # this overwrites the total flux for ρθ (using methods below)
        numerical_flux,
        bctype.ρθ,
        balance_law,
        args...
    )

end

# customized methods form specifying total normal fluxes
function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive},
    bctype::CoupledPrimaryBoundary,
    balance_law::ModelSetup,
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

    fluxᵀn.ρθ = calculate_land_sfc_fluxes(balance_law, state_prognostic⁻, state_auxiliary⁻) / balance_law.parameters.c_p # W/m^2

end

function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive}, bctype::Insulating, balance_law::ModelSetup, fluxᵀn::Vars{S}, _...,) where {S}
    FT = eltype(fluxᵀn)
    fluxᵀn.ρθ = FT(0)
end


# # customized flux for Insulating boundary
# numerical_boundary_flux_second_order!(
#     numerical_flux::Union{PenaltyNumFluxDiffusive},
#     bctype::FluidBC{Impenetrable{FreeSlip},Insulating},
#     balance_law::ModelSetup,
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
    bl::ModelSetup,
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