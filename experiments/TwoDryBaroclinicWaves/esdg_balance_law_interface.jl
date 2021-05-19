import ClimateMachine.BalanceLaws:
    # declaration
    vars_state,
    # initialization
    nodal_init_state_auxiliary!,
    init_state_prognostic!,
    init_state_auxiliary!,
    # rhs computation
    compute_gradient_argument!,
    compute_gradient_flux!,
    flux_first_order!,
    flux_second_order!,
    source!,
    # boundary conditions
    boundary_conditions,
    boundary_state!


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

import ClimateMachine.DGMethods.NumericalFluxes.numerical_boundary_flux_first_order!

struct DryReferenceState{TP}
    temperature_profile::TP
end

using Unitful
"""

    Defines an equation set (balance law) for testing coupling

 Defines kernels to evaluates RHS of

    ```

DryAtmosModel: 

    ∂ρe
    -- + ∇ • ( ρu * (ρe + p) / ρ ) = 0
    ∂t    

    ∂ρu
    -- + ∇ • ( p * I + ρu ⊗ ρu / ρ) = - 2Ω × ρu
    ∂t    

    ∂ρ
    -- + ∇ • ( ρu ) = 0
    ∂t    


DryAtmosLinearModel:

    ∂ρe
    -- + ∇ • ( (ρeᵣ + pᵣ) / ρᵣ * ρu ) = - state.ρu' * ∇Φ
    ∂t    

    ∂ρu
    -- + ∇ • ( -0) = 0
    ∂t    

    ∂ρ
    -- + ∇ • ( ρu ) = 0
    ∂t    


(NB: `split_explicit_implicit == false` by default)

Boundary conditions:

External:
- first order:
    state⁺.ρ = state⁻.ρ
    state⁺.ρu -= 2 * dot(state⁻.ρu, n) .* SVector(n)
    state⁺.ρe = state⁻.ρe ( ≡ fluxᵀn.ρe = 0)
    aux⁺.Φ = aux⁻.Φ

PrimaryCoupledBoundary:
    - state⁺.ρ = state⁻.ρ
    - state⁺.ρu -= 2 * dot(state⁻.ρu, n) .* SVector(n)
    - fluxᵀn.ρe = (ρu * (ρe + p) / ρ )ᵀ ⋅ n = λ(ρe - ρe_secondary) 
    - aux⁺.Φ = aux⁻.Φ

SecondaryCoupledBoundary:
    - state⁺.ρ = state⁻.ρ
    - state⁺.ρu -= 2 * dot(state⁻.ρu, n) .* SVector(n)
    - fluxᵀn.ρe = (ρu * (ρe + p) / ρ )ᵀ ⋅ n = F_prescribed = ∫ F_accum dt / Δt_coupler
    - aux⁺.Φ = aux⁻.Φ


"""


"""
    Declaration of state variables

    vars_state returns a NamedTuple of data types.
"""
function vars_state(m::Union{DryAtmosModel,DryAtmosLinearModel}, st::Auxiliary, FT)
    @vars begin
        x::FT
        y::FT
        z::FT
        Φ::FT
        ∇Φ::SVector{3, FT} # TODO: only needed for the linear model
        ref_state::vars_state(m, m.physics.ref_state, st, FT)
        ρe_secondary::FT  # stores opposite face for primary (atmospheric import)
        F_ρe_prescribed::FT # stores prescribed flux for secondary (ocean import)
    end
end

vars_state(::Union{DryAtmosModel,DryAtmosLinearModel}, ::DryReferenceState, ::Auxiliary, FT) =
    @vars(T::FT, p::FT, ρ::FT, ρe::FT)
vars_state(::Union{DryAtmosModel,DryAtmosLinearModel}, ::NoReferenceState, ::Auxiliary, FT) = @vars()

function vars_state(::Union{DryAtmosModel,DryAtmosLinearModel}, ::Prognostic, FT)
    @vars begin
        ρ::FT
        ρu::SVector{3, FT}
        ρe::FT
        F_ρe_accum::FT
    end
end

function vars_state(::DryAtmosModel, ::Entropy, FT)
    @vars begin
        ρ::FT
        ρu::SVector{3, FT}
        ρe::FT
        Φ::FT
    end
end

"""
    Initialization of state variables

    init_state_xyz! sets up the initial fields within our state variables
    (e.g., prognostic, auxiliary, etc.), however it seems to not initialized
    the gradient flux variables by default.
"""
function init_state_prognostic!(
        model::Union{DryAtmosModel,DryAtmosLinearModel},
        state::Vars,
        aux::Vars,
        localgeo,
        t
    )
    x = aux.x
    y = aux.y
    z = aux.z

    parameters = model.parameters
    ic = model.initial_conditions

    if !isnothing(ic)
        state.ρ  = ic.ρ(parameters, x, y, z)
        state.ρu = ic.ρu(parameters, x, y, z)
        state.ρe = ic.ρe(parameters, x, y, z)
        state.F_ρe_accum = 0
    end

    return nothing
end

function nodal_init_state_auxiliary!(
    m::Union{DryAtmosModel,DryAtmosLinearModel},
    state_auxiliary,
    tmp,
    geom,
)
    init_state_auxiliary!(m, m.physics.orientation, state_auxiliary, geom)
    init_state_auxiliary!(m, m.physics.ref_state, state_auxiliary, geom)

    state_auxiliary.ρe_secondary = 0
    state_auxiliary.F_ρe_prescribed = 0
end

function init_state_auxiliary!(
    ::Union{DryAtmosModel,DryAtmosLinearModel},
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
    ::Union{DryAtmosModel,DryAtmosLinearModel},
    ::FlatOrientation,
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

function init_state_auxiliary!(
    ::Union{DryAtmosModel,DryAtmosLinearModel},
    ::NoReferenceState,
    state_auxiliary,
    geom,
) end

function init_state_auxiliary!(
    m::Union{DryAtmosModel,DryAtmosLinearModel},
    ref_state::DryReferenceState,
    state_auxiliary,
    geom,
)
    FT = eltype(state_auxiliary)
    z = altitude(m, m.physics.orientation, geom)
    T, p = ref_state.temperature_profile(param_set, z)

    _R_d::FT = R_d(param_set)
    ρ = p / (_R_d * T)
    Φ = state_auxiliary.Φ
    ρu = SVector{3, FT}(0, 0, 0)

    state_auxiliary.ref_state.T = T
    state_auxiliary.ref_state.p = p
    state_auxiliary.ref_state.ρ = ρ
    state_auxiliary.ref_state.ρe = totalenergy(ρ, ρu, p, Φ)
end

"""
    LHS computations
"""
@inline function flux_first_order!(
    model::Union{DryAtmosModel,DryAtmosLinearModel},
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
    direction,
)

    lhs = model.physics.lhs
    ntuple(Val(length(lhs))) do s
        Base.@_inline_meta
        calc_flux!(flux, lhs[s], state, aux, t)
    end
end

"""
    RHS computations
"""
function source!(m::DryAtmosModel, source, state_prognostic, state_auxiliary, _...)
    sources = m.physics.sources

    ntuple(Val(length(sources))) do s
        Base.@_inline_meta
        calc_force!(source, sources[s], state_prognostic, state_auxiliary)

    end

    source.F_ρe_accum = (state_prognostic.ρe - state_auxiliary.ρe_secondary) * m.parameters.λ_coupler
end

function source!(
    m::DryAtmosLinearModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
    ::NTuple{1, Dir},
) where {Dir <: Direction}
    sources = m.physics.sources

    ntuple(Val(length(sources))) do s
        Base.@_inline_meta
        calc_force!(source, sources[s], state, aux)
    end
end

"""
    Boundary conditions
"""

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

boundary_conditions(model::Union{DryAtmosModel,DryAtmosLinearModel}) = model.boundary_conditions

function boundary_state!(
    ::NumericalFluxFirstOrder,
    bctype,
    ::Union{DryAtmosModel,DryAtmosLinearModel},
    state⁺,
    aux⁺,
    n,
    state⁻,
    aux⁻,
    _...,
)
    state⁺.ρ = state⁻.ρ
    state⁺.ρu -= 2 * dot(state⁻.ρu, n) .* SVector(n)
    state⁺.ρe = state⁻.ρe
    aux⁺.Φ = aux⁻.Φ
end

function boundary_state!(
    nf::NumericalFluxSecondOrder,
    bc,
    lm::Union{DryAtmosModel,DryAtmosLinearModel},
    args...,
)
    nothing
end

"""
    Utils
"""
function vertical_unit_vector(::Union{DryAtmosModel,DryAtmosLinearModel}, aux::Vars)
    FT = eltype(aux)
    aux.∇Φ / FT(grav(param_set))
end

function altitude(::Union{DryAtmosModel,DryAtmosLinearModel}, ::SphericalOrientation, geom)
    FT = eltype(geom)
    _planet_radius::FT = planet_radius(param_set)
    norm(geom.coord) - _planet_radius
end

function altitude(::Union{DryAtmosModel,DryAtmosLinearModel}, ::FlatOrientation, geom)
    @inbounds geom.coord[3]
end

# Helper functions to communicate between components before and after timestepping

"""
function preB(csolver)
    - saves and regrids the EnergyA couplerfield (i.e. regridded mA.state.ρe[mA.boundary]) on coupler grid to mB.state.ρe_secondary[mB.boundary] on the domainB grid
    csolver::CplSolver
"""
function preB(csolver)
    mA = csolver.component_list.domainA.component_model
    mB = csolver.component_list.domainB.component_model
    # Set boundary SST used in atmos to SST of ocean surface at start of coupling cycle.
    mB.odesolver.rhs!.state_auxiliary.ρe_secondary[mB.boundary] .= 
        CouplerMachine.coupler_get(csolver.coupler, :EnergyA, mB.grid.numerical, DateTime(0), u"J")
    # Set atmos boundary flux accumulator to 0.
    mB.state.F_ρe_accum .= 0

    @info(
        "preatmos",
        time = csolver.t, #* "/" * mB.time.finish ,
        total_θ_atmos = weightedsum(mA.state, 1),
        total_θ_ocean = weightedsum(mA.state, 1),
        total_θ = weightedsum(mA.state, 1) + weightedsum(mA.state, 1),
        atmos_θ_surface_max = maximum(mA.state.ρe[mA.boundary]),
        ocean_θ_surface_max = maximum(mA.state.ρe[mA.boundary]),
    )

    # isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.A[csolver.steps] = weightedsum(mA.state, 1)
    # isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.B[csolver.steps] = weightedsum(mB.state, 1)
end

"""
function postB(csolver)
    - updates couplerfield EnergyFluxB with mB.state.F_ρe_accum[mB.boundary] regridded to the coupler grid, and updates the coupler time
    csolver::CplSolver
"""
function postB(csolver)
    mA = csolver.component_list.domainA.component_model
    mB = csolver.component_list.domainB.component_model
    # Pass atmos exports to "coupler" namespace
    # 1. Save mean θ flux at the Atmos boundary during the coupling period
    CouplerMachine.coupler_put!(csolver.coupler, :EnergyFluxB, mB.state.F_ρe_accum[mB.boundary] ./ csolver.dt,
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
    #         mean(mB.state.F_ρe_accum[mB.boundary]) * 1e6 * 1e6,
    #     F_accum_max = maximum(mB.state.F_accum[mB.boundary]),
    #     F_avg_max = maximum(mB.state.F_accum[mB.boundary] ./ csolver.dt),
    #     atmos_θ_surface_max = maximum(mB.state.θ[mB.boundary]),
    #     ocean_θ_surface_max = maximum(mA.state.θ[mA.boundary]),
    # )
end

"""
function preA(csolver)
    - saves and regrids the EnergyFluxB couplerfield (i.e. regridded mB.state.F_ρe_accum[mB.boundary]) to ρe_prescribed[mA.boundary] on the domainA grid
    csolver::CplSolver
"""
function preA(csolver)
    mA = csolver.component_list.domainA.component_model
    mB = csolver.component_list.domainB.component_model
    # Set mean air-sea theta flux

    mA.odesolver.rhs!.state_auxiliary.F_ρe_prescribed[mA.boundary] .= 
        CouplerMachine.coupler_get(csolver.coupler, :EnergyFluxB, mA.grid, DateTime(0), u"J")
    # Set ocean boundary flux accumulator to 0. (this isn't used)
    mA.state.F_ρe_accum .= 0

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
    - updates couplerfield EnergyA with mA.state.ρe[mA.boundary] regridded to the coupler grid, and updates the coupler time
    - updates EnergyA with mA.state.ρe[mA.boundary], and the coupler time
    csolver::CplSolver
"""
function postA(csolver)
    mA = csolver.component_list.domainA.component_model
    mB = csolver.component_list.domainB.component_model
    # @info(
    #     "postocean",
    #     time = csolver.t + csolver.dt,
    #     ocean_θ_surface_max = maximum(mA.state.θ[mA.boundary]),
    #     ocean_θ_surface_min = maximum(mA.state.θ[mA.boundary]),
    # )

    # Pass ocean exports to "coupler" namespace
    #  1. Ocean SST (value of θ at z=0)
    CouplerMachine.coupler_put!(csolver.coupler, :EnergyA, mA.state.ρe[mA.boundary], mA.grid.numerical, DateTime(0), u"J")
end





#__________________

function numerical_boundary_flux_first_order!(
    numerical_flux::NumericalFluxFirstOrder,
    bctype::ExteriorBoundary,
    balance_law::BalanceLaw,
    fluxᵀn::Vars{S},
    normal_vector::SVector,
    state_prognostic⁻::Vars{S},
    state_auxiliary⁻::Vars{A},
    state_prognostic⁺::Vars{S},
    state_auxiliary⁺::Vars{A},
    t,
    direction,
    state1⁻::Vars{S},
    aux1⁻::Vars{A},
) where {S, A}

    boundary_state!(
        numerical_flux,
        bctype,
        balance_law,
        state_prognostic⁺,
        state_auxiliary⁺,
        normal_vector,
        state_prognostic⁻,
        state_auxiliary⁻,
        t,
        state1⁻,
        aux1⁻,
    )

    numerical_flux_first_order!(
        numerical_flux,
        balance_law,
        fluxᵀn,
        normal_vector,
        state_prognostic⁻,
        state_auxiliary⁻,
        state_prognostic⁺,
        state_auxiliary⁺,
        t,
        direction,
    )

end

function numerical_boundary_flux_first_order!(
    numerical_flux::NumericalFluxFirstOrder,
    bctype::CoupledPrimaryBoundary,
    balance_law::BalanceLaw,
    fluxᵀn::Vars{S},
    normal_vector::SVector,
    state_prognostic⁻::Vars{S},
    state_auxiliary⁻::Vars{A},
    state_prognostic⁺::Vars{S},
    state_auxiliary⁺::Vars{A},
    t,
    direction,
    state1⁻::Vars{S},
    aux1⁻::Vars{A},
) where {S, A}

    boundary_state!(
        numerical_flux,
        bctype,
        balance_law,
        state_prognostic⁺,
        state_auxiliary⁺,
        normal_vector,
        state_prognostic⁻,
        state_auxiliary⁻,
        t,
        state1⁻,
        aux1⁻,
    )


    numerical_flux_first_order!(
        numerical_flux,
        balance_law,
        fluxᵀn,
        normal_vector,
        state_prognostic⁻,
        state_auxiliary⁻,
        state_prognostic⁺,
        state_auxiliary⁺,
        t,
        direction,
    )

    fluxᵀn.ρe =
        (state_prognostic⁻.ρe - state_auxiliary⁺.ρe_secondary) * balance_law.parameters.λ_coupler
end
function numerical_boundary_flux_first_order!(
    numerical_flux::NumericalFluxFirstOrder,
    bctype::CoupledPrimaryBoundary,
    balance_law::BalanceLaw,
    fluxᵀn::Vars{S},
    normal_vector::SVector,
    state_prognostic⁻::Vars{S},
    state_auxiliary⁻::Vars{A},
    state_prognostic⁺::Vars{S},
    state_auxiliary⁺::Vars{A},
    t,
    direction,
    state1⁻::Vars{S},
    aux1⁻::Vars{A},
) where {S, A}

    boundary_state!(
        numerical_flux,
        bctype,
        balance_law,
        state_prognostic⁺,
        state_auxiliary⁺,
        normal_vector,
        state_prognostic⁻,
        state_auxiliary⁻,
        t,
        state1⁻,
        aux1⁻,
    )


    numerical_flux_first_order!(
        numerical_flux,
        balance_law,
        fluxᵀn,
        normal_vector,
        state_prognostic⁻,
        state_auxiliary⁻,
        state_prognostic⁺,
        state_auxiliary⁺,
        t,
        direction,
    )
    
    fluxᵀn.ρe = -state_auxiliary⁺.F_ρe_prescribed
end

