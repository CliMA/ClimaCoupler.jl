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
import ClimateMachine.NumericalFluxes:
    numerical_boundary_flux_first_order!

struct DryReferenceState{TP}
    temperature_profile::TP
end

"""
    Declaration of state variables
    vars_state returns a NamedTuple of data types.
"""
function vars_state(m::DryAtmosModel, st::Auxiliary, FT)
    @vars begin
        x::FT
        y::FT
        z::FT
        Φ::FT
        ∇Φ::SVector{3, FT} # TODO: only needed for the linear model
        ref_state::vars_state(m, m.physics.ref_state, st, FT)
        T_sfc::FT  # Ocean surface temperature
    end
end

vars_state(::DryAtmosModel, ::DryReferenceState, ::Auxiliary, FT) =
    @vars(T::FT, p::FT, ρ::FT, ρu::SVector{3, FT}, ρe::FT, ρq::FT)
vars_state(::DryAtmosModel, ::NoReferenceState, ::Auxiliary, FT) = @vars()

function vars_state(::DryAtmosModel, ::Prognostic, FT)
    @vars begin
        ρ::FT
        ρu::SVector{3, FT}
        ρe::FT
        ρq::FT
        F_ρe_accum::FT # accumulated energy flux
    end
end

"""
    Initialization of state variables
    init_state_xyz! sets up the initial fields within our state variables
    (e.g., prognostic, auxiliary, etc.), however it seems to not initialized
    the gradient flux variables by default.
"""
function init_state_prognostic!(
        model::DryAtmosModel,
        state::Vars,
        aux::Vars,
        localgeo,
        t
    )
    x = aux.x
    y = aux.y
    z = aux.z

    parameters = model.physics.parameters
    ic = model.initial_conditions

    # TODO!: Set to 0 by default or assign IC
    if !isnothing(ic)
        state.ρ  = ic.ρ(parameters, x, y, z)
        state.ρu = ic.ρu(parameters, x, y, z)
        state.ρe = ic.ρe(parameters, x, y, z)
        state.ρq = ic.ρq(parameters, x, y, z)
    end

    state.F_ρe_accum = 0    

    return nothing
end

function nodal_init_state_auxiliary!(
    model::DryAtmosModel,
    state_auxiliary,
    tmp,
    geom,
)
    init_state_auxiliary!(model, model.physics.orientation, state_auxiliary, geom)
    init_state_auxiliary!(model, model.physics.ref_state, state_auxiliary, geom)

    state_auxiliary.T_sfc = model.physics.parameters.T_h
end

function init_state_auxiliary!(
    model::DryAtmosModel,
    ::SphericalOrientation,
    state_auxiliary,
    geom,
)
    g = model.physics.parameters.g

    r = norm(geom.coord)
    state_auxiliary.x = geom.coord[1]
    state_auxiliary.y = geom.coord[2]
    state_auxiliary.z = geom.coord[3]
    state_auxiliary.Φ = g * r
    state_auxiliary.∇Φ = g * geom.coord / r
end

function init_state_auxiliary!(
    model::DryAtmosModel,
    ::FlatOrientation,
    state_auxiliary,
    geom,
)
    g = model.physics.parameters.g

    FT = eltype(state_auxiliary)
    
    r = geom.coord[3]
    state_auxiliary.x = geom.coord[1]
    state_auxiliary.y = geom.coord[2]
    state_auxiliary.z = geom.coord[3]
    state_auxiliary.Φ = g * r
    state_auxiliary.∇Φ = SVector{3, FT}(0, 0, g)
end

function init_state_auxiliary!(
    ::DryAtmosModel,
    ::NoReferenceState,
    state_auxiliary,
    geom,
) end

function init_state_auxiliary!(
    model::DryAtmosModel,
    ref_state::DryReferenceState,
    state_auxiliary,
    geom,
)
    orientation = model.physics.orientation   
    R_d         = model.physics.parameters.R_d
    γ           = model.physics.parameters.γ
    Φ           = state_auxiliary.Φ

    FT = eltype(state_auxiliary)

    # Calculation of a dry reference state
    z = altitude(model, orientation, geom)
    T, p = ref_state.temperature_profile(model.physics.parameters, z)
    ρ  = p / R_d / T
    ρu = SVector{3, FT}(0, 0, 0)
    ρe = p / (γ - 1) + dot(ρu, ρu) / 2ρ + ρ * Φ
    ρq = FT(0)

    state_auxiliary.ref_state.T  = T
    state_auxiliary.ref_state.p  = p
    state_auxiliary.ref_state.ρ  = ρ
    state_auxiliary.ref_state.ρu = ρu
    state_auxiliary.ref_state.ρe = ρe
    state_auxiliary.ref_state.ρq = ρq    
end

"""
    LHS computations
"""
@inline function flux_first_order!(
    model::DryAtmosModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
    direction,
)
    lhs = model.physics.lhs
    physics = model.physics
    
    ntuple(Val(length(lhs))) do s
        Base.@_inline_meta
        calc_component!(flux, lhs[s], state, aux, physics)
    end
end

"""
    RHS computations
"""
function source!(
    model::DryAtmosModel, 
    source, 
    state_prognostic, 
    state_auxiliary, 
    _...
) # TODO: please enable t import here :)
    sources = model.physics.sources
    physics = model.physics

    ntuple(Val(length(sources))) do s
        Base.@_inline_meta
        calc_component!(source, sources[s], state_prognostic, state_auxiliary, physics)
    end

end

"""
    Boundary conditions with defaults
"""
boundary_conditions(model::DryAtmosModel) = model.boundary_conditions

function boundary_state!(_...)
    nothing
end

"""
    Utils
"""
function altitude(model::DryAtmosModel, ::SphericalOrientation, geom)
    return norm(geom.coord) - model.physics.parameters.a
end

function altitude(::DryAtmosModel, ::FlatOrientation, geom)
    @inbounds geom.coord[3]
end

"""
function preAtmos(csolver)
    - saves and regrids the `SeaSurfaceTemerature` couplerfield (i.e. regridded `mOcean.state.T_sfc[mOcean.boundary]``) on coupler grid to `mAtmos.aux.T_sfc[mAtmos.boundary]` on the `mAtmos` grid
    - resets the accumulated flux `F_ρe_accum` in `mAtmos` to 0
    csolver::CplSolver
"""
function preAtmos(csolver)
    mOcean = csolver.component_list.domainOcean.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    # Set boundary T_sfc used in atmos at the start of the coupling cycle.
    mAtmos.odesolver.rhs!.state_auxiliary.T_sfc[mAtmos.boundary] .= 
        coupler_get(csolver.coupler, :SeaSurfaceTemerature, mAtmos.grid.numerical, DateTime(0), u"K")
    # Set atmos boundary flux accumulator to 0.
    mAtmos.state.F_ρe_accum .= 0

    # get indicies of the required variables from the `mAtmos` MPIStateArray
    idx = varsindex(vars(mAtmos.state), :ρe)[1]
    idx_rho = varsindex(vars(mAtmos.state), :ρ)[1]

    # Calculate, print and log domain-averaged energy
    FT = eltype(mAtmos.state)

    model_type = cpl_solver.component_list.domainAtmos.component_model.model[1] # TODO: bl should be more accessible
    p = model_type.model.physics.parameters # maybe the parameter list could be saved in the coupler? (this would requre all the compute kernels below to import the coupler)
    nel = mAtmos.grid.resolution.elements.vertical
    po = mAtmos.grid.resolution.polynomial_order.vertical

    E_Ocean = weightedsum(mOcean.state, 1) .* p.ρ_o .* p.h_o .* p.c_o  # J / m^2
    E_Atmos = weightedsum(mAtmos.state, idx) .* p.cp_d .* p.H ./  nel ./ (po+1) ./ (po+2) # J / m^2 

    @info(
        "preatmos",
        time = string(csolver.t) * "/" * string(mAtmos.time.finish),
        total_energyA = E_Ocean,
        total_energyB = E_Atmos,
        total_energy = E_Atmos + E_Ocean,
        Ocean_T_sfc = maximum(mOcean.state.T_sfc[mOcean.boundary]),
        atmos_ρe_sfc = maximum(mAtmos.state.ρe[mAtmos.boundary]),
    )

    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.A[csolver.steps] = E_Ocean
    isnothing(csolver.fluxlog) ? nothing : csolver.fluxlog.B[csolver.steps] = E_Atmos

end

"""
function postAtmos(csolver)
    - updates couplerfield `EnergyFluxAtmos` with mAtmos.state.F_ρe_accum[mAtmos.boundary] regridded to the coupler grid, and updates the coupler time
    csolver::CplSolver
"""
function postAtmos(csolver)
    
    mOcean = csolver.component_list.domainOcean.component_model
    mAtmos = csolver.component_list.domainAtmos.component_model
    # @show mAtmos.state.ρe
    # Pass atmos exports to "coupler" namespace
    # 1. Save mean θ flux at the Atmos boundary during the coupling period

    @info(
        "postatmos",
        time = string(csolver.t) * "/" * string(mAtmos.time.finish),
        Ocean_T_sfc = maximum(mOcean.state.T_sfc[mOcean.boundary]),
        atmos_ρe_sfc = maximum(mAtmos.state.ρe[mAtmos.boundary]),
    )

    coupler_put!(csolver.coupler, :EnergyFluxAtmos, mAtmos.state.F_ρe_accum[mAtmos.boundary] ./ csolver.dt,
        mAtmos.grid.numerical, DateTime(0), u"J")
    #@show mAtmos.state.F_ρe_accum[simAtmos.boundary]
    end