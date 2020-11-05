# Model A
#
using ClimateMachine.DGMethods.NumericalFluxes:
    NumericalFluxFirstOrder, NumericalFluxSecondOrder, NumericalFluxGradient

using ClimateMachine.DGMethods.NumericalFluxes:
    CentralNumericalFluxGradient, CentralNumericalFluxSecondOrder

using LinearAlgebra: I, dot, Diagonal

using ClimateMachine.BalanceLaws:
   BalanceLaw, Prognostic, Auxiliary, Gradient, GradientFlux

using ClimateMachine.VariableTemplates

import ClimateMachine.DGMethods:
     init_state_auxiliary!, update_auxiliary_state!, update_auxiliary_state_gradient!, vars_state, VerticalDirection, boundary_state!, compute_gradient_flux!, init_state_prognostic!, flux_first_order!, flux_second_order!, source!, wavespeed, compute_gradient_argument!

using StaticArrays

"""
 AModel{M} <: BalanceLaw
"""

# Create a new linear model instance
abstract type AbstractAModel <: BalanceLaw end
struct AModel{FT} <: AbstractAModel 
 function AModel{FT}(
  ;
 ) where {FT<: AbstractFloat}
    return new{FT}(
    )
 end
end

function ADGModel(
   bl::AModel,
   grid,
   nfnondiff,
   nfdiff,
   gnf;
   kwargs...,
   )

   modeldata=()

   return DGModel(bl,grid,nfnondiff,nfdiff,gnf;kwargs...,modeldata=modeldata,)
end

"""
 Set model state variables and operators
"""

# State variable and initial value, just one for now, θ
##
vars_state(m::AModel, ::Prognostic, FT) = @vars(θ::FT)

init_state_prognostic!( args...) = ( init_state_cp!( args... ) )
function init_state_cp!( m::IVDCModel, Q::Vars, A::Vars, coords, t,)
  @inbounds begin
    Q.θ = 10.
  end
  return nothing
end
##

##
vars_state(m::IVDCModel, ::Auxiliary, FT) = @vars(θ_init::FT)

function init_state_auxiliary!(m::IVDCModel,A::Vars, _...)
  @inbounds begin
    A.θ_init = -0
  end
  return nothing
end
##
#

# Variables and operations used in differentiating first derivatives
##
vars_state(m::IVDCModel, ::Gradient, FT) = @vars(∇θ::FT, ∇θ_init::FT,)
@inline function compute_gradient_argument!(
    m::IVDCModel,
    G::Vars,
    Q::Vars,
    A,
    t,
)
    G.∇θ = Q.θ
    G.∇θ_init = A.θ_init

    return nothing
end
##
#

# Variables and operations used in differentiating second derivatives
##
vars_state(m::IVDCModel, ::GradientFlux, FT) = @vars(κ∇θ::SVector{3, FT})
@inline function compute_gradient_flux!(
    m::IVDCModel,
    D::Vars,
    G::Grad,
    Q::Vars,
    A::Vars,
    t,
)
    κ = diffusivity_tensor(m, G.∇θ_init[3])
    D.κ∇θ = -κ * G.∇θ 
    return nothing
end
##

##
## Set vertical diffusivity profile based on vertical hydrography profile
@inline function diffusivity_tensor(m::IVDCModel, ∂θ∂z)
    κᶻ = m.κᶻ
    κᶜ = m.κᶜ
    ∂θ∂z < 0 ? κ = (@SVector [0, 0, κᶜ]) : κ = (@SVector [0, 0, κᶻ])
    return Diagonal(-κ)
end
##
#

# Function to apply I to state variable

##
@inline function source!(
    m::IVDCModel,
    S::Vars,
    Q::Vars,
    D::Vars,
    A::Vars,
    t,
    direction,
)
    ivdc_dt = m.dt
    @inbounds begin
     S.θ = Q.θ/ivdc_dt
    end

    return nothing
end
##
#

# Numerical fluxes and boundaries

##
function flux_first_order!(::IVDCModel, _...) end

function flux_second_order!(
    ::IVDCModel,
    F::Grad,
    Q::Vars,
    D::Vars,
    H::Vars,
    A::Vars,
    t,
)
    F.θ += D.κ∇θ
end

function wavespeed(m::IVDCModel, n⁻, _...)
    C = abs(SVector(m.cʰ, m.cʰ, m.cᶻ)' * n⁻)
    return C
end
##

function boundary_state!(
    nf::Union{NumericalFluxFirstOrder, NumericalFluxGradient, CentralNumericalFluxGradient},
    m::IVDCModel,
    Q⁺,
    A⁺,
    n,
    Q⁻,
    A⁻,
    bctype,
    t,
    _...,
)
    Q⁺.θ = Q⁻.θ
    return nothing
end

###    From -  function numerical_boundary_flux_gradient! , DGMethods/NumericalFluxes.jl
###    boundary_state!(
###        numerical_flux,
###        balance_law,
###        state_conservative⁺,
###        state_auxiliary⁺,
###        normal_vector,
###        state_conservative⁻,
###        state_auxiliary⁻,
###        bctype,
###        t,
###        state1⁻,
###        aux1⁻,
###    )

function boundary_state!(
    nf::Union{NumericalFluxSecondOrder,CentralNumericalFluxSecondOrder},
    m::IVDCModel,
    Q⁺,
    D⁺,
    A⁺,
    n⁻,
    Q⁻,
    D⁻,
    A⁻,
    bctype,
    t,
    _...,
)
    Q⁺.θ = Q⁻.θ
    D⁺.κ∇θ = n⁻ * -0
    # D⁺.κ∇θ = Q/(rho*cp)   ( at top boundary for heat flux )
    return nothing
end

###    boundary_state!(
###        numerical_flux,
###        balance_law,
###        state_conservative⁺,
###        state_gradient_flux⁺,
###        state_auxiliary⁺,
###        normal_vector,
###        state_conservative⁻,
###        state_gradient_flux⁻,
###        state_auxiliary⁻,
###        bctype,
###        t,
###        state1⁻,
###        diff1⁻,
###        aux1⁻,
###    )

###    boundary_flux_second_order!(
###        numerical_flux,
###        balance_law,
###        Grad{S}(flux),
###        state_conservative⁺,
###        state_gradient_flux⁺,
###        state_hyperdiffusive⁺,
###        state_auxiliary⁺,
###        normal_vector,
###        state_conservative⁻,
###        state_gradient_flux⁻,
###        state_hyperdiffusive⁻,
###        state_auxiliary⁻,
###        bctype,
###        t,
###        state1⁻,
###        diff1⁻,
###        aux1⁻,
###    )

