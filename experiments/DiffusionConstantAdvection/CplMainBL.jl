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


# θ: J / m^3
# ∇θ: J / m^4
# κ: (m^2/s)
# κ∇θ: J/(m^2 s) = W/m^2 = J/m^3 * m/s

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
    RusanovNumericalFlux,
    numerical_boundary_flux_second_order!, 
    numerical_flux_second_order!, 
    numerical_boundary_flux_first_order!, 
    numerical_flux_first_order!


using ClimateMachine.VariableTemplates

using LinearAlgebra
using StaticArrays
using ClimateMachine.Orientations
"""
    CplMainBL <: BalanceLaw
    - used in both the atmos and ocean component models
"""
struct CplMainBL{BLP, BCS, PS, O} <: BalanceLaw
    bl_prop::BLP  # initial condition, ...
    boundaryconditions::BCS
    param_set::PS
    orientation::O
end

l_type = CplMainBL

"""
function vars_state(bl::l_type, ::Prognostic, FT)
- Declare prognostic state variables
    - `θ`: tracer
    - `F_accum`: "Shadow" variable used to capture boundary fluxes that we want to accumulate over 
    a single timestep and export to coupling as time integrals. We use a shadow variable
    because we want to integrate over whatever timestepper is being used.
    Eventually we should have the ability to potentially use a 2d field here.
    The shadow variable needs to be zeroed at the start of each
    coupling cycle for a component.
    - `u`: advective velocity - currently constant, but kept in prognostic variables to enable easy overintegration  
"""
function vars_state(bl::l_type, ::Prognostic, FT)
    @vars begin
        θ::FT
        F_accum::FT 
        u::SVector{3, FT}
    end
end

"""
function vars_state(bl::l_type, st::Auxiliary, FT)
- Declare Aaxiliary state variables

  `npt`::Int    # no. nodes
  `elnum`::Int  # no. elems

  `xc`::FT      # Cartesian x
  `yc`::FT      # Cartesian y
  `zc`::FT      # Cartesian z
  
  `θⁱⁿⁱᵗ`::FT   # unused in default setup
  `θ_secondary`::FT  # stores opposite face for primary (atmospheric import)
  `F_prescribed`::FT # stores prescribed flux for secondary (ocean import)
  `orientation`::vars_state(bl.orientation, st, FT)
"""
function vars_state(bl::l_type, st::Auxiliary, FT)
    @vars begin
        npt::Int    # no. nodes
        elnum::Int  # no. elems

        xc::FT      # Cartesian x
        yc::FT      # Cartesian y
        zc::FT      # Cartesian z
        
        θⁱⁿⁱᵗ::FT   # unused in default setup
        θ_secondary::FT  # stores opposite face for primary (atmospheric import)
        F_prescribed::FT # stores prescribed flux for secondary (ocean import)
        orientation::vars_state(bl.orientation, st, FT)
    end
end

vars_state(::Orientation, ::Auxiliary, FT) = @vars(Φ::FT, ∇Φ::SVector{3, FT})

"""
function vars_state(bl::l_type, ::Gradient, FT)
  - Pre-gradient computation variables 
  `∇θ`::FT 
  `∇θⁱⁿⁱᵗ`::FT # unused in default setup
"""
function vars_state(bl::l_type, ::Gradient, FT)
    @vars begin
        ∇θ::FT 
        ∇θⁱⁿⁱᵗ::FT # unused in default setup
    end
end

"""
function vars_state(bl::l_type, ::GradientFlux, FT)
  - Post-gradient computation variable
  `κ∇θ`::SVector{3, FT}
"""
function vars_state(bl::l_type, ::GradientFlux, FT)
    @vars begin
        κ∇θ::SVector{3, FT}
    end
end

"""
function init_state_prognostic!
  Point-wise initialization of prognostic state variables
  `bl`::l_type, balance law
  `Q`::Vars, state variables
  `A`::Vars, auxilliary variables
  `geom`::LocalGeometry, 
  `FT`,
"""

function init_state_prognostic!(
    bl::l_type,
    Q::Vars,
    A::Vars,
    geom::LocalGeometry,
    FT,
)
    npt = A.npt
    elnum = A.elnum
    x = A.xc 
    y = A.yc
    z = A.zc

    Q.θ = bl.bl_prop.init_theta(npt, elnum, x, y, z)
    Q.u = bl.bl_prop.init_u(npt, elnum, x, y, z)
        
    Q.F_accum = 0
    nothing
end

"""
function nodal_init_state_auxiliary!
    Point-wise initialization of auxiliary state variables
    `bl`::l_type, balance law
    `A`::Vars, auxilliary variables
    `tmp`::Vars,
    `geom`::LocalGeometry,

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

    # This is necesary for ∇Φ, which is used to get the `vertical_unit_vector`
    FT = eltype(A)
    A.orientation.Φ = grav(param_set) * (norm(geom.coord) - planet_radius(param_set))

    nothing
end

""""
function init_state_auxiliary!
    - array-wise initialization of auxiliary state variables
    `model`::l_type, balance law
    `state_auxiliary`::MPIStateArray,
    `grid`,
    `direction`,
"""
function init_state_auxiliary!(
    model::l_type,
    state_auxiliary::MPIStateArray,
    grid,
    direction,
)

    init_state_auxiliary!(
        bl,
        (bl, A, tmp, geom) -> nodal_init_state_auxiliary!(bl,A,tmp,geom),
        state_auxiliary,
        grid,
        direction,
    )    

    # update ∇Φ in state_auxiliary.orientation.∇Φ (necessary for the unit vector notmal to the vertical surfces of the reference element)
    auxiliary_field_gradient!(
        model,
        state_auxiliary,
        ("orientation.∇Φ",),
        state_auxiliary,
        ("orientation.Φ",),
        grid,
        direction,
    )
end

"""
function source!
  - for prognostic state external sources
  - for recording boundary flux terms into shadow variables for export to coupler
  bl::l_type, balance law
  S::Vars, source variables
  Q::Vars, prognostic variables
  G::Vars, pre-gradient variables
  A::Vars, auxiliary variables 
  """
function source!(bl::l_type, S::Vars, Q::Vars, G::Vars, A::Vars, _...)
    #S.θ=bl.bl_prop.source_theta(Q.θ,A.npt,A.elnum,A.xc,A.yc,A.zc,A.θ_secondary)
    # Record boundary condition fluxes as needed by adding to shadow
    # prognostic variable

    S.F_accum = (Q.θ - A.θ_secondary) * bl.bl_prop.coupling_lambda()     
    nothing
end

"""
function compute_gradient_argument!
  Set values to have gradients computed.
  `bl`::l_type, balance law
  `G`::Vars, pre-gradient variables
  `Q`::Vars, prognostic variables 
  `A`::Vars, auxiliary variables 
  `t`, time
"""
function compute_gradient_argument!(bl::l_type, G::Vars, Q::Vars, A::Vars, t)
    G.∇θ = Q.θ
    G.∇θⁱⁿⁱᵗ = A.θⁱⁿⁱᵗ
    nothing
end

"""
function compute_gradient_flux!
  Compute diffusivity tensor times computed gradient to give net gradient flux.
  bl::l_type, balance law
  GF::Vars, post-gradient (gradient flux) variables
  G::Grad, pre-gradient variables
  Q::Vars, prognostic variables
  A::Vars, auxiliary variables 
  t,
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
    F =
        bl.bl_prop.calc_diff_flux(G.∇θ, A.npt, A.elnum, A.xc, A.yc, A.zc)

    GF.κ∇θ = F
    nothing
end

"""
function flux_second_order!(
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
abstract type AbstractCouplerBoundary end
  Define boundary condition flags/types to iterate over, for now keep it simple.
- 3 AbstractCouplerBoundary types:
    1. ExteriorBoundary
    2. CoupledPrimaryBoundary
    3. CoupledSecondaryBoundary

"""
abstract type AbstractCouplerBoundary end

"""
struct PenaltyNumFluxDiffusive <: NumericalFluxSecondOrder end
- PenaltyNumFluxDiffusive: additional NumericalFluxSecondOrder to account for an additional penalty term
"""
struct PenaltyNumFluxDiffusive <: NumericalFluxSecondOrder end

# ## 1. ExteriorBoundary
"""
struct ExteriorBoundary <: AbstractCouplerBoundary end
  Zero normal gradient boundary condition.
"""
struct ExteriorBoundary <: AbstractCouplerBoundary end

# ## 2. CoupledPrimaryBoundary
"""
struct CoupledPrimaryBoundary <: AbstractCouplerBoundary end
# # compute flux based on opposite face
# # also need to accumulate net flux across boundary
"""
struct CoupledPrimaryBoundary <: AbstractCouplerBoundary end

## 3. CoupledSecondaryBoundary
"""
struct CoupledSecondaryBoundary <: AbstractCouplerBoundary end
# # use prescribed flux computed in primary
"""
struct CoupledSecondaryBoundary  <: AbstractCouplerBoundary end

# first order BC (Neumann, Insulating = no bdry flux)
function boundary_state!(
    nF::Union{NumericalFluxFirstOrder,NumericalFluxGradient},
    bc::Union{ExteriorBoundary,CoupledPrimaryBoundary,CoupledSecondaryBoundary},
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
    return nothing
end

# 2nd order BCs (Neumann, Insulating = no bdry flux) - not called if using customized numerical_boundary_flux_second_order! below
# function boundary_state!(
#     nf,
#     bc::Union{ExteriorBoundary,CoupledPrimaryBoundary,CoupledSecondaryBoundary},
#     m::l_type,
#     state⁺::Vars,
#     diff⁺::Vars,
#     hyperdiff⁺::Vars,
#     aux⁺::Vars,
#     n⁻,
#     state⁻::Vars,
#     diff⁻::Vars,
#     hyperdiff⁻::Vars,
#     aux⁻::Vars,
#     t,
#     _...,
# )
#     # Apply Neumann BCs
#     FT = eltype(state⁺)
#     diff⁺.κ∇θ = n⁻ * FT(0.0)
# end

# customized flux for ExteriorBoundary
function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive},
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

    fluxᵀn.θ = 0.0

end

# customized flux for CoupledPrimaryBoundary
function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive},
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
        balance_law.bl_prop.coupling_lambda() # W/m^2   # T / s

end

# customized flux for CoupledSecondaryBoundary
function numerical_boundary_flux_second_order!(
    numerical_flux::Union{PenaltyNumFluxDiffusive},
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

    fluxᵀn.θ = -state_auxiliary⁺.F_prescribed # T m /s^2 

end

function wavespeed(bl::l_type, _...)
    # Used in Rusanov term.
    # Only active if there is a flux first order term?
    bl.bl_prop.get_wavespeed()
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

    # Remove surplus flux normal to the reference element surface
    #k̂ = vertical_unit_vector(bl, A)
    #F.θ += (SDiagonal(1, 1, 1) - k̂ * k̂')*Q.u * Q.θ  
    nothing
end


function boundary_conditions(bl::l_type, _...)
    bl.boundaryconditions
end

"""
  Set a default set of properties and their default values
  - init_aux_geom   :: function to initialize geometric terms stored in aux.
  - init_theta      :: function to set initial θ values.
  - source_theta    :: function to add a source term to θ.
  - calc_diff_flux  :: function to set calculation of diffusive fluxes
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

    calc_diff_flux(_...) = (return 0.0, 0.0, 0.0)
    bl_prop = (bl_prop..., calc_diff_flux = calc_diff_flux)

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

# Helper functions to communicate between components before and after timestepping

"""
function preatmos(csolver)
    - saves and regrids the OceanSST field (mO.state.θ[mO.boundary]) to θ_secondary[mA.boundary] on the atmos grid
    csolver::CplSolver
"""
function preatmos(csolver)
    mA = csolver.component_list.atmosphere.component_model
    mO = csolver.component_list.ocean.component_model
    # Set boundary SST used in atmos to SST of ocean surface at start of coupling cycle.
    mA.discretization.state_auxiliary.θ_secondary[mA.boundary] .= 
        CouplerMachine.coupler_get(csolver.coupler, :Ocean_SST, mA.grid, DateTime(0), u"°C")
    # Set atmos boundary flux accumulator to 0.
    mA.state.F_accum .= 0

    @info(
        "preatmos",
        endtime = simulation.simtime[2],
        time = csolver.t,
        total_θ_atmos = weightedsum(mA.state, 1),
        total_θ_ocean = weightedsum(mO.state, 1),
        total_θ = weightedsum(mA.state, 1) + weightedsum(mO.state, 1),
        atmos_θ_surface_max = maximum(mA.state.θ[mA.boundary]),
        ocean_θ_surface_max = maximum(mO.state.θ[mO.boundary]),
    )

    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.A[csolver.steps] = weightedsum(mA.state, 1)
    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.O[csolver.steps] = weightedsum(mO.state, 1)
end

"""
function preatmos(csolver)
    - updates Atmos_MeanAirSeaθFlux with mA.state.F_accum[mA.boundary], and the coupler time
    csolver::CplSolver
"""
function postatmos(csolver)
    mA = csolver.component_list.atmosphere.component_model
    mO = csolver.component_list.ocean.component_model
    # Pass atmos exports to "coupler" namespace
    # 1. Save mean θ flux at the Atmos boundary during the coupling period
    CouplerMachine.coupler_put!(csolver.coupler, :Atmos_MeanAirSeaθFlux, mA.state.F_accum[mA.boundary] ./ csolver.dt,
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

"""
function preocean(csolver)
    - saves and regrids the Atmos_MeanAirSeaθFlux field (mA.state.F_accum[mA.boundary]) to F_prescribed[mO.boundary] on the ocean grid
    csolver::CplSolver
"""
function preocean(csolver)
    mA = csolver.component_list.atmosphere.component_model
    mO = csolver.component_list.ocean.component_model
    # Set mean air-sea theta flux
    mO.discretization.state_auxiliary.F_prescribed[mO.boundary] .= 
        CouplerMachine.coupler_get(csolver.coupler, :Atmos_MeanAirSeaθFlux, mO.grid, DateTime(0), u"°C")
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

"""
function postocean(csolver)
    - updates Ocean_SST with mO.state.θ[mO.boundary], and the coupler time
    csolver::CplSolver
"""
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
    CouplerMachine.coupler_put!(csolver.coupler, :Ocean_SST, mO.state.θ[mO.boundary], mO.grid, DateTime(0), u"°C")
end


#end