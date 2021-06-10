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
    m::ModelSetup,
    state_auxiliary,
    tmp,
    geom,
)
    init_state_auxiliary!(m, m.physics.orientation, state_auxiliary, geom)

    state_auxiliary.T_sfc = 0
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
end

function init_state_prognostic!(model::ModelSetup, state::Vars, aux::Vars, localgeo, t)
    x = aux.x
    y = aux.y
    z = aux.z

    parameters = model.parameters
    ic = model.initial_conditions

    state.ρ = ic.ρ(parameters, x, y, z)
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
    source.F_ρθ_accum = calculate_land_sfc_fluxes(model, state, aux, t) 
    
    return nothing
end

"""
    calculate_land_sfc_fluxes(model::ModelSetup, state, aux, t)
- calculate furface fluxes using the bulk gradient diffusion theory
"""
function calculate_land_sfc_fluxes(model::ModelSetup, state, aux, t; MO_params = nothing) # should pass in the coupler state (also move to coupler), so can access states of both models derectly -e.g. callback?
    
    FT = eltype(state)
    p = model.parameters

    p_sfc = p.pₒ
    θ_a = state.ρθ / state.ρ 
    T_a = θ_a
    T_sfc = aux.T_sfc 
    θ_sfc = T_sfc  

    R_SW = (FT(1)-p.α) * p.τ * p.F_sol * (FT(1) .+ sin(t * 2π / p.diurnal_period) )
    R_LW = p.ϵ * (p.σ * θ_sfc .^ 4 - p.F_a)
    SH   = state.ρ * p.cp_d * p.g_a * (T_sfc - T_a) # g_a could be substituded by g_a=f(MO_params) (see Bonan P90), but first need to modify SurfaceFluxes.jl
      
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
    calc_boundary_state!(numerical_flux, bc.ρθ, model, diffusion, args...)
end

function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive},
    bctype::FluidBC,
    balance_law::ModelSetup,
    args... 
) where {S, D, A, HD}

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
    F_tot = calculate_land_sfc_fluxes(balance_law, state_prognostic⁻, state_auxiliary⁻, t)  # W/m^2
    fluxᵀn.ρθ = - F_tot  / balance_law.parameters.cp_d

end

function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive}, bctype::Insulating, balance_law::ModelSetup, fluxᵀn::Vars{S}, _...,) where {S}
    FT = eltype(fluxᵀn)
    fluxᵀn.ρθ = FT(0)
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
    tau = FT(0.0) # decide if want to use this penalty
    Fᵀn .+= tau * (parent(state⁻) - parent(state⁺))
end

"""
function preAtmos(csolver)
    - saves and regrids the `LandSurfaceTemerature` couplerfield (i.e. regridded `mLand.state.T_sfc[mLand.boundary]``) on coupler grid to `mAtmos.aux.T_sfc[mAtmos.boundary]` on the `mAtmos` grid
    - resets the accumulated flux `F_ρθ_accum` in `mAtmos` to 0
    csolver::CplSolver
"""
function preAtmos(csolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    # Set boundary T_sfc used in atmos at the start of the coupling cycle.
    mAtmos.odesolver.rhs!.state_auxiliary.T_sfc[mAtmos.boundary] .= 
        coupler_get(csolver.coupler, :LandSurfaceTemerature, mAtmos.grid.numerical, DateTime(0), u"K")
    # Set atmos boundary flux accumulator to 0.
    mAtmos.state.F_ρθ_accum .= 0

    # get indicies of the required variables from the `mAtmos` MPIStateArray
    idx = varsindex(vars(mAtmos.state), :ρθ)[1]
    idx_rho = varsindex(vars(mAtmos.state), :ρ)[1]

    # Calculate, print and log domain-averaged energy
    FT = eltype(mAtmos.state)
    p = cpl_solver.component_list.domainAtmos.component_model.model.parameters # maybe the parameter list could be saved in the coupler? (this would requre all the compute kernels below to import the coupler)
    nel = mAtmos.grid.resolution.elements.vertical
    po = mAtmos.grid.resolution.polynomial_order.vertical

    E_Land = weightedsum(mLand.state, 1) .* p.ρ_s .* p.h_s .* p.c_s  # J / m^2
    E_Atmos = weightedsum(mAtmos.state, idx) .* p.cp_d .* p.zmax ./  nel ./ (po+1) ./ (po+2) # J / m^2 

    @info(
        "preatmos",
        time = string(csolver.t) * "/" * string(mAtmos.time.finish),
        total_energyA = E_Land,
        total_energyB = E_Atmos,
        total_energy = E_Atmos + E_Land,
        land_T_sfc = maximum(mLand.state.T_sfc[mLand.boundary]),
        atmos_ρθ_sfc = maximum(mAtmos.state.ρθ[mAtmos.boundary]),
    )

    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.A[csolver.steps] = E_Land
    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.B[csolver.steps] = E_Atmos

end

"""
function postAtmos(csolver)
    - updates couplerfield `EnergyFluxAtmos` with mAtmos.state.F_ρθ_accum[mAtmos.boundary] regridded to the coupler grid, and updates the coupler time
    csolver::CplSolver
"""
function postAtmos(csolver)
    mLand = csolver.component_list.domainLand.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    # Pass atmos exports to "coupler" namespace
    # 1. Save mean θ flux at the Atmos boundary during the coupling period
    coupler_put!(csolver.coupler, :EnergyFluxAtmos, mAtmos.state.F_ρθ_accum[mAtmos.boundary] ./ csolver.dt,
        mAtmos.grid.numerical, DateTime(0), u"J")
end