"""
    module CplTestingBL
 Defines an equation set (balance law) for testing coupling

 Defines kernels to evaluates RHS of

    ```
    ∂θ
    --  = - ∇ • ( u θ + κ ∇ θ ) = - ∇ • F
    ∂t
    ```

    Where
    
     - `θ` is the tracer (e.g. potential temperature)
     - `u` is the initial advection velocity (constant)
     - `κ` is the diffusivity tensor

 subject to specified prescribed gradient or flux bc's

 - [`CplMainBL`](@ref)
     balance law struct created by this module

"""
#module CplMainBL

export CplMainBL
export PenaltyNumFluxDiffusive
export ExteriorBoundary
export CoupledPrimaryBoundary, CoupledSecondaryBoundary

using ClimateMachine.BalanceLaws:
    Auxiliary, BalanceLaw, Gradient, GradientFlux, Prognostic

import ClimateMachine.BalanceLaws:
    boundary_conditions,
    boundary_state!,
    compute_gradient_argument!,
    compute_gradient_flux!,
    flux_first_order!,
    flux_second_order!,
    init_state_prognostic!,
    nodal_init_state_auxiliary!,
    source!,
    vars_state,
    wavespeed

using ClimateMachine.Mesh.Geometry: LocalGeometry
using ClimateMachine.MPIStateArrays

using ClimateMachine.DGMethods.NumericalFluxes:
    CentralNumericalFluxGradient,
    CentralNumericalFluxSecondOrder,
    NumericalFluxFirstOrder,
    NumericalFluxSecondOrder,
    RusanovNumericalFlux

import ClimateMachine.DGMethods.NumericalFluxes:
    numerical_boundary_flux_second_order!, numerical_flux_second_order!, numerical_boundary_flux_first_order!, numerical_flux_first_order!


using ClimateMachine.VariableTemplates

using LinearAlgebra
using StaticArrays

"""
    CplMainBL <: BalanceLaw
    - used in both the atmos and ocean component models
"""
struct CplMainBL{BLP, BCS, PS} <: BalanceLaw
    bl_prop::BLP  # initial condition, ...
    boundaryconditions::BCS
    param_set::PS
end

l_type = CplMainBL

"""
Declare prognostic state variables
    - θ
    - F_accum: "Shadow" variable used to capture boundary fluxes that we want to accumulate over 
    a single timestep and export to coupling as time integrals. We use a shadow variable
    because we want to integrate over whatever timestepper is being used.
    Eventually we should have the ability to potentially use a 2d field here.
    The shadow variable needs to be zeroed at the start of each
    coupling cycle for a component.
    - u: advective velocity - currently constant, but kept in prognostic variables to enable
    easier use of overintegration  
"""
function vars_state(bl::l_type, ::Prognostic, FT)
    @vars begin
        θ::FT
        F_accum::FT 
        u::SVector{3, FT}
    end
end

"""
Declare Aaxiliary state variables
  - for array index and real world coordinates and for
  θ value at reference time used to compute κ.
"""
function vars_state(bl::l_type, ::Auxiliary, FT)
    @vars begin
        npt::Int    # no. nodes
        elnum::Int  # no. elems

        xc::FT      # Cartesian x
        yc::FT      # Cartesian y
        zc::FT      # Cartesian z
        
        θⁱⁿⁱᵗ::FT   # unused in default setup
        θ_secondary::FT  # stores opposite face for primary (atmospheric import)
        F_prescribed::FT # stores prescribed flux for secondary (ocean import)
    end
end

"""
  Gradient computation stage input (and output) variable symbols
"""
function vars_state(bl::l_type, ::Gradient, FT)
    @vars begin
        ∇θ::FT 
        ∇θⁱⁿⁱᵗ::FT # unused in default setup
    end
end

function vars_state(bl::l_type, ::GradientFlux, FT)
    @vars begin
        κ∇θ::SVector{3, FT}
    end
end

"""
  Initialize prognostic state variables
"""
function init_state_prognostic!(
    bl::l_type,
    Q::Vars,
    A::Vars,
    geom::LocalGeometry,
    FT,
)
    npt = getproperty(geom, :n)
    elnum = getproperty(geom, :e)
    x = A.xc 
    y = A.yc
    z = A.zc

    Q.θ = bl.bl_prop.init_theta(npt, elnum, x, y, z)
    Q.F_accum = 0
    Q.u = bl.bl_prop.init_u(npt, elnum, x, y, z)
    nothing
end

"""
  Initialize auxiliary state variables
"""
function nodal_init_state_auxiliary!(
    bl::l_type,
    A::Vars,
    tmp::Vars,
    geom::LocalGeometry,
    _...,
)
    npt = getproperty(geom, :n)
    elnum = getproperty(geom, :e)
    x = geom.coord[1]
    y = geom.coord[2]
    z = geom.coord[3]

    A.npt, A.elnum, A.xc, A.yc, A.zc =
        bl.bl_prop.init_aux_geom(npt, elnum, x, y, z)
    
    A.θⁱⁿⁱᵗ = 0
    A.θ_secondary = 0
    A.F_prescribed = 0
    nothing
end


"""
Compute Kernels:
"""
#====

Atmos

----

Land


====#


"""
  Set source terms 
  - for prognostic state external sources
  - for recording boundary flux terms into shadow variables for export to coupler
"""
function source!(bl::l_type, S::Vars, Q::Vars, G::Vars, A::Vars, _...)
    #S.θ=bl.bl_prop.source_theta(Q.θ,A.npt,A.elnum,A.xc,A.yc,A.zc,A.θ_secondary)
    # Record boundary condition fluxes as needed by adding to shadow
    # prognostic variable

    S.F_accum = (Q.θ - A.θ_secondary) * bl.bl_prop.coupling_lambda()
    nothing
end

"""
  Set values to have gradients computed.
"""
function compute_gradient_argument!(bl::l_type, G::Vars, Q::Vars, A::Vars, t)
    G.∇θ = Q.θ
    G.∇θⁱⁿⁱᵗ = A.θⁱⁿⁱᵗ
    nothing
end

"""
  Compute diffusivity tensor times computed gradient to give net gradient flux.
"""
function compute_gradient_flux!(
    bl::l_type,
    GF::Vars,
    G::Grad,
    Q::Vars,
    A::Vars,
    t,
)
    # "Non-linear" form (for time stepped)
    ### κ¹,κ²,κ³=bl.bl_prop.calc_kappa_diff(G.∇θ,A.npt,A.elnum,A.xc,A.yc,A.zc)
    # "Linear" form (for implicit)
    κ¹, κ², κ³ =
        bl.bl_prop.calc_kappa_diff(G.∇θⁱⁿⁱᵗ, A.npt, A.elnum, A.xc, A.yc, A.zc)
    # Maybe I should pass both G.∇θ and G.∇θⁱⁿⁱᵗ?

    # Messy diffusion tensor rotation - should find code for doing this succinctly!
    # we want κ¹, κ² to act in lat and lon directions respectively i.e. tangential 
    # to shell, we want κ³ to act normal to shell i.e. z or r direction.
    ## don't use regular norm(), its not great on performance
    mynorm(x::Float64,y::Float64,z::Float64)  = ( x^2 + y^2 + z^2 ) ^ 0.5
    ## excessively type stable forms!
    r̂ⁿᵒʳᵐ(x::Float64,y::Float64,z::Float64) = mynorm(x,y,z) ≈ 0 ? 1 : mynorm(x, y, z)^(-1)
    ϕ̂ⁿᵒʳᵐ(x::Float64,y::Float64,z::Float64) = mynorm(x,y,Float64(0)) ≈ 0 ? 1 : ( mynorm(x, y, z) * mynorm(x, y, Float64(0)) )^(-1)
    λ̂ⁿᵒʳᵐ(x::Float64,y::Float64,z::Float64) = mynorm(x,y,Float64(0)) ≈ 0 ? 1 :   mynorm(x, y, Float64(0))^(-1)
    r̂(x::Float64,y::Float64,z::Float64) = r̂ⁿᵒʳᵐ(x,y,z) * @SVector([x, y, z])
    ϕ̂(x::Float64,y::Float64,z::Float64) = ϕ̂ⁿᵒʳᵐ(x,y,z) * @SVector [x*z, y*z, -(x^2 + y^2)]
    λ̂(x::Float64,y::Float64,z::Float64) = λ̂ⁿᵒʳᵐ(x,y,z) * @SVector [-y, x, 0]

    # Sphere surface bits 
    x=A.xc
    y=A.yc
    z=A.zc
    Gλ = G.∇θ'*λ̂(x,y,z)
    Fλ = κ¹*Gλ*λ̂(x,y,z)
    Gϕ = G.∇θ'*ϕ̂(x,y,z)
    Fϕ = κ²*Gϕ*ϕ̂(x,y,z)
    Gr = G.∇θ'*r̂(x,y,z)
    Fr = κ³*Gr*r̂(x,y,z)
    GF.κ∇θ = Diagonal(@SVector([1, 1, 1])) * (Fλ + Fϕ + Fr)

    # Original cartesian form
    #GF.κ∇θ = Diagonal(@SVector([κ¹, κ², κ³])) * G.∇θ   
    #GF.κ∇θ = Diagonal(@SVector([κ¹, κ¹, κ¹])) * G.∇θ   
    nothing
end

"""
  Pass flux components for second order term into update kernel.
"""
function flux_second_order!(
    bl::l_type,
    F::Grad,
    Q::Vars,
    GF::Vars,
    H::Vars,
    A::Vars,
    t,
)
    F.θ += GF.κ∇θ
    nothing
end

# Boundary conditions
"""
  Define boundary condition flags/types to iterate over, for now keep it simple.
- 3 AbstractCouplerBoundary types:
    1. ExteriorBoundary
    2. CoupledPrimaryBoundary
    3. CoupledSecondaryBoundary
- PenaltyNumFluxDiffusive: additional NumericalFluxSecondOrder to account for an additional penalty term
"""
struct PenaltyNumFluxDiffusive <: NumericalFluxSecondOrder end

abstract type AbstractCouplerBoundary end

# ## 1. ExteriorBoundary
# flux is 0 across the boundary
struct ExteriorBoundary <: AbstractCouplerBoundary end

"""
  Zero normal gradient boundary condition.
"""

function boundary_state!(
    nF::Union{CentralNumericalFluxGradient},
    bc::ExteriorBoundary,
    bl::l_type,
    Q⁺::Vars,
    A⁺::Vars,
    n,
    Q⁻::Vars,
    A⁻::Vars,
    t,
    _...,
)
    Q⁺.θ = Q⁻.θ
    nothing
end
function numerical_boundary_flux_second_order!(
    numerical_flux::Union{NumericalFluxSecondOrder},
    bctype::ExteriorBoundary,
    balance_law::l_type,
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
    
    fluxᵀn.θ += 0
end


# ## 2. CoupledPrimaryBoundary
# # compute flux based on opposite face
# # also need to accumulate net flux across boundary
 struct CoupledPrimaryBoundary <: AbstractCouplerBoundary end

function boundary_state!(
    nF::Union{CentralNumericalFluxGradient},
    bc::CoupledPrimaryBoundary,
    bl::l_type,
    Q⁺::Vars,
    A⁺::Vars,
    n,
    Q⁻::Vars,
    A⁻::Vars,
    t,
    _...,
)
    
    Q⁺.θ = Q⁻.θ # Q⁺.θ=A⁺.θ_secondary
    
    nothing
end

function numerical_boundary_flux_second_order!(
    numerical_flux::Union{NumericalFluxSecondOrder},
    bctype::CoupledPrimaryBoundary,
    balance_law::l_type,
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
    fluxᵀn.θ =
        (state_prognostic⁻.θ - state_auxiliary⁺.θ_secondary) *
        balance_law.bl_prop.coupling_lambda() # W/m^2

end


# θ: J / m^3
# ∇θ: J / m^4
# κ: (m^2/s)
# κ∇θ: J/(m^2 s) = W/m^2 = J/m^3 * m/s

#  - clean up and write primer
#  - imports and exports
#  - add boundary auxiliary interface
#  - add boundary flux accumulation interface
#  - vector quantity
#  - run on sphere to check normal terms are all correct
#  - single stack integration
#  - Held-Suarez with wind stress and heat

## 3. CoupledSecondaryBoundary
# use prescribed flux computed in primary
struct CoupledSecondaryBoundary  <: AbstractCouplerBoundary end
function boundary_state!(
    nF::Union{CentralNumericalFluxGradient},
    bc::CoupledSecondaryBoundary,
    bl::l_type,
    Q⁺::Vars,
    A⁺::Vars,
    n,
    Q⁻::Vars,
    A⁻::Vars,
    t,
    _...,
)
    Q⁺.θ = Q⁻.θ
    nothing
end
function numerical_boundary_flux_second_order!(
    numerical_flux::Union{NumericalFluxSecondOrder},
    bctype::CoupledSecondaryBoundary,
    balance_law::l_type,
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
    fluxᵀn.θ = -state_auxiliary⁺.F_prescribed
end

function wavespeed(bl::l_type, _...)
    # Used in Rusanov term.
    # Only active if there is a flux first order term?
    bl.bl_prop.get_wavespeed()
end

"""
  Penalty flux formulation of second order numerical flux. This formulation
  computes the CentralNumericalFluxSecondOrder term first (which is just the average
  of the + and - fluxes and an edge), and then adds a "penalty" flux that relaxes
  the edge state + and - toward each other.
"""
function numerical_flux_second_order!(
    ::PenaltyNumFluxDiffusive,
    bl::l_type,
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
    tau = bl.bl_prop.get_penalty_tau()
    Fᵀn .+= tau * (parent(state⁻) - parent(state⁺))
end

"""
  First order flux for advection equation

  NumericalFluxFirstOrder is an abstract type that currently generalizes
  RusanovNumericalFlux, CentralNumericalFluxFirstOrder, RoeNumericalFlux, HLLCNumericalFlux.
"""
function flux_first_order!(
    bl::l_type,
    F::Grad,
    Q::Vars,
    A::Vars,
    t::Real,
    directions,
)
    F.θ += Q.u * Q.θ  
    nothing
end
function boundary_state!(
    nf::Union{NumericalFluxFirstOrder},
    bc::Union{CoupledPrimaryBoundary,CoupledSecondaryBoundary,ExteriorBoundary},
    atmos::l_type,
    Q⁺::Vars,
    A⁺::Vars,
    n,
    Q⁻::Vars,
    A⁻::Vars,
    t,
    state_int⁻,
    aux_int⁻,
)
    Q⁺.θ = Q⁻.θ
    nothing
end

function boundary_conditions(bl::l_type, _...)
    bl.boundaryconditions
end

# #this one is for the general adv only test
# function boundary_state!(
#     nF,
#     bc,
#     bl::l_type,
#     Q⁺,
#     state_gradient_flux⁺,
#     state_hyperdiffusive⁺,
#     state_auxiliary⁺,
#     normal_vector,
#     Q⁻,
#     state_gradient_flux⁻,
#     state_hyperdiffusive⁻,
#     state_auxiliary⁻,
#     t,
#     state1⁻,
#     diff1⁻,
#     aux1⁻,
# )
#     Q⁺.θ = Q⁻.θ # Q⁺.θ=A⁺.θ_secondary
    
#     nothing
# end

"""
  Set a default set of properties and their default values
  - init_aux_geom   :: function to initialize geometric terms stored in aux.
  - init_theta      :: function to set initial θ values.
  - source_theta    :: function to add a source term to θ.
  - calc_kappa_diff :: function to set diffusion coeffiecient(s).
  - get_wavespeed   :: function to return a wavespeed for Rusanov computations (there aren't any in this model)
  - get_penalty_tau :: function to set timescale on which to bring state+ and state- together
  - theta_shadow_boundary_flux :: function to set boundary flux into shadow variable for passing to coupler
"""
function prop_defaults()
    bl_prop = NamedTuple()

    function init_aux_geom(npt, elnum, x, y, z)
        return npt, elnum, x, y, z
    end
    # init_aux_geom(_...)=(return 0., 0., 0., 0., 0.)
    bl_prop = (bl_prop..., init_aux_geom = init_aux_geom)

    init_theta(_...) = (return 0.0)
    bl_prop = (bl_prop..., init_theta = init_theta)

    source_theta(_...) = (return 0.0)
    bl_prop = (bl_prop..., source_theta = source_theta)

    calc_kappa_diff(_...) = (return 0.0, 0.0, 0.0)
    bl_prop = (bl_prop..., calc_kappa_diff = calc_kappa_diff)

    get_wavespeed(_...) = (return 0.0)
    bl_prop = (bl_prop..., get_wavespeed = get_wavespeed)

    get_penalty_tau(_...) = (return 1.0)
    bl_prop = (bl_prop..., get_penalty_tau = get_penalty_tau)

    theta_shadow_boundary_flux(_...) = (return 0.0)
    bl_prop =
        (bl_prop..., theta_shadow_boundary_flux = theta_shadow_boundary_flux)

    coupling_lambda(_...) = (return 0.0)
    bl_prop = (bl_prop..., coupling_lambda = coupling_lambda)

    init_u(_...) = (return 0.0)
    bl_prop = (bl_prop..., init_u = init_u)

    bl_prop = (bl_prop..., LAW = CplMainBL)
end

"""
  Helper functions to communicate between components before and after timestepping
  - preatmos(csolver)
  - postatmos(csolver)
  - preocean(csolver)
  - postocean(csolver)
"""

function preatmos(csolver)
    mA = csolver.component_list.atmosphere.component_model
    mO = csolver.component_list.ocean.component_model
    # Set boundary SST used in atmos to SST of ocean surface at start of coupling cycle.
    mA.discretization.state_auxiliary.θ_secondary[mA.boundary] .= 
        coupler_get(csolver.coupler, :Ocean_SST, mA.grid, DateTime(0), u"°C")
    # Set atmos boundary flux accumulator to 0.
    mA.state.F_accum .= 0

    @info(
        "preatmos",
        time = csolver.t,
        total_θ_atmos = weightedsum(mA.state, 1),
        total_θ_ocean = weightedsum(mO.state, 1),
        total_θ = weightedsum(mA.state, 1) + weightedsum(mO.state, 1),
        atmos_θ_surface_max = maximum(mA.state.θ[mA.boundary]),
        ocean_θ_surface_max = maximum(mO.state.θ[mO.boundary]),
    )
end

function postatmos(csolver)
    mA = csolver.component_list.atmosphere.component_model
    mO = csolver.component_list.ocean.component_model
    # Pass atmos exports to "coupler" namespace
    # 1. Save mean θ flux at the Atmos boundary during the coupling period
    coupler_put!(csolver.coupler, :Atmos_MeanAirSeaθFlux, mA.state.F_accum[mA.boundary] ./ csolver.dt,
        mA.grid, DateTime(0), u"°C")

    @info(
        "postatmos",
        time = time = csolver.t + csolver.dt,
        total_θ_atmos = weightedsum(mA.state, 1),
        total_θ_ocean = weightedsum(mO.state, 1),
        total_F_accum = mean(mA.state.F_accum[mA.boundary]) * 1e6 * 1e6,
        total_θ =
            weightedsum(mA.state, 1) +
            weightedsum(mO.state, 1) +
            mean(mA.state.F_accum[mA.boundary]) * 1e6 * 1e6,
        F_accum_max = maximum(mA.state.F_accum[mA.boundary]),
        F_avg_max = maximum(mA.state.F_accum[mA.boundary] ./ csolver.dt),
        atmos_θ_surface_max = maximum(mA.state.θ[mA.boundary]),
        ocean_θ_surface_max = maximum(mO.state.θ[mO.boundary]),
    )
end

function preocean(csolver)
    mA = csolver.component_list.atmosphere.component_model
    mO = csolver.component_list.ocean.component_model
    # Set mean air-sea theta flux
    mO.discretization.state_auxiliary.F_prescribed[mO.boundary] .= 
        coupler_get(csolver.coupler, :Atmos_MeanAirSeaθFlux, mO.grid, DateTime(0), u"°C")
    # Set ocean boundary flux accumulator to 0. (this isn't used)
    mO.state.F_accum .= 0

    @info(
        "preocean",
        time = csolver.t,
        F_prescribed_max =
            maximum(mO.discretization.state_auxiliary.F_prescribed[mO.boundary]),
        F_prescribed_min =
            maximum(mO.discretization.state_auxiliary.F_prescribed[mO.boundary]),
        ocean_θ_surface_max = maximum(mO.state.θ[mO.boundary]),
        ocean_θ_surface_min = maximum(mO.state.θ[mO.boundary]),
    )
end

function postocean(csolver)
    mA = csolver.component_list.atmosphere.component_model
    mO = csolver.component_list.ocean.component_model
    @info(
        "postocean",
        time = csolver.t + csolver.dt,
        ocean_θ_surface_max = maximum(mO.state.θ[mO.boundary]),
        ocean_θ_surface_min = maximum(mO.state.θ[mO.boundary]),
    )

    # Pass ocean exports to "coupler" namespace
    #  1. Ocean SST (value of θ at z=0)
    coupler_put!(csolver.coupler, :Ocean_SST, mO.state.θ[mO.boundary], mO.grid, DateTime(0), u"°C")
end


#end