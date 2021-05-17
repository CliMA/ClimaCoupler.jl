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

struct DryReferenceState{TP}
    temperature_profile::TP
end

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
    - saves and regrids the EnergyFluxA field (mA.state.ρe[mA.boundary]) to mB.state.ρe_accum[mB.boundary] on the atmos grid
    csolver::CplSolver
"""
function preB(csolver)
    mA = csolver.component_list.domainA.component_model
    mB = csolver.component_list.domainB.component_model
    # Set boundary SST used in atmos to SST of ocean surface at start of coupling cycle.
    mB.discretization.state_auxiliary.ρe_secondary[mB.boundary] .= 
        CouplerMachine.get(csolver.coupler, :EnergyA, mB.grid, DateTime(0), u"J")
    # Set atmos boundary flux accumulator to 0.
    mB.state.ρe_accum .= 0

    @info(
        "preatmos",
        endtime = simulation.simtime[2],
        time = csolver.t,
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
function preatmos(csolver)
    - updates Atmos_MeanAirSeaθFlux with mA.state.F_accum[mA.boundary], and the coupler time
    csolver::CplSolver
"""
function postB(csolver)
    mA = csolver.component_list.domainA.component_model
    mB = csolver.component_list.domainB.component_model
    # Pass atmos exports to "coupler" namespace
    # 1. Save mean θ flux at the Atmos boundary during the coupling period
    CouplerMachine.put!(csolver.coupler, :EnergyFluxB, mB.state.ρe_accum[mB.boundary] ./ csolver.dt,
        mB.grid, DateTime(0), u"J")

    # @info(
    #     "postatmos",
    #     time = time = csolver.t + csolver.dt,
    #     total_θ_atmos = weightedsum(mB.state, 1),
    #     total_θ_ocean = weightedsum(mA.state, 1),
    #     total_F_accum = mean(mB.state.F_accum[mB.boundary]) * 1e6 * 1e6,
    #     total_θ =
    #         weightedsum(mB.state, 1) +
    #         weightedsum(mA.state, 1) +
    #         mean(mB.state.ρe_accum[mB.boundary]) * 1e6 * 1e6,
    #     F_accum_max = maximum(mB.state.F_accum[mB.boundary]),
    #     F_avg_max = maximum(mB.state.F_accum[mB.boundary] ./ csolver.dt),
    #     atmos_θ_surface_max = maximum(mB.state.θ[mB.boundary]),
    #     ocean_θ_surface_max = maximum(mA.state.θ[mA.boundary]),
    # )
end

"""
function preocean(csolver)
    - saves and regrids the EnergyFluxB field (mB.state.ρe_accum[mB.boundary]) to ρe_accum[mA.boundary] on the lower-domain grid
    csolver::CplSolver
"""
function preA(csolver)
    mA = csolver.component_list.domainA.component_model
    mB = csolver.component_list.domainB.component_model
    # Set mean air-sea theta flux
    mA.discretization.state_auxiliary.ρe_accum[mA.boundary] .= 
        CouplerMachine.get(csolver.coupler, :EnergyFluxB, mA.grid, DateTime(0), u"J")
    # Set ocean boundary flux accumulator to 0. (this isn't used)
    mA.state.ρe_accum .= 0

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
function postocean(csolver)
    - updates Ocean_SST with mA.state.θ[mA.boundary], and the coupler time
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
    CouplerMachine.put!(csolver.coupler, :EnergyFluxA, mA.state.ρe[mA.boundary], mA.grid, DateTime(0), u"J")
end
