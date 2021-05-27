using ClimateMachine, MPI
using ClimateMachine.DGMethods.NumericalFluxes
using ClimateMachine.DGMethods

using ClimateMachine.ODESolvers
using ClimateMachine.ODESolvers: solve!

using ClimateMachine.Atmos: SphericalOrientation, latitude, longitude, AtmosFilterPerturbations

using CLIMAParameters
using CLIMAParameters.Planet: MSLP, R_d, day, grav, Omega, planet_radius

using ClimateMachine.Mesh.Filters

using ClimateMachine.BalanceLaws:
    BalanceLaw,
    AbstractStateType,
    Prognostic,
    Auxiliary,
    Gradient,
    GradientFlux,
    GradientLaplacian,
    Hyperdiffusive,
    UpwardIntegrals,
    DownwardIntegrals,
    vars_state,
    number_states

using CouplerMachine

using Unitful
using Dates: DateTime
using Statistics
using StaticArrays

struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

ClimateMachine.init()
FT = Float64

# Shared functions
include("domains.jl")
include("abstractions.jl")
include("callbacks.jl")

# Main balance law and its components
include("CplMainBL.jl")

nstepsA = 10
nstepsO = 5

#  Background atmos and ocean diffusivities
const κᵃʰ = FT(1e4) * 0.0
const κᵃᶻ = FT(1e-1)
const κᵒʰ = FT(1e3) * 0.0
const κᵒᶻ = FT(1e-4)
const τ_airsea = FT(60 * 86400)
const L_airsea = FT(500)
const λ_airsea = FT(L_airsea / τ_airsea)
coupling_lambda() = (λ_airsea)
const u_max = FT(1e-5) # max. advective velocity in radians

function main(::Type{FT}) where {FT}
    
    # Domain
    # confusing name - better might be to use something like DeepSphericalShellDomain directly?
    ΩO = DeepSphericalShellDomain(radius = FT(planet_radius(param_set)) - FT(4e3), height = FT(4e3))
    ΩA = DeepSphericalShellDomain(radius = FT(planet_radius(param_set)) , height = FT(4e3))

    # Grid
    nelem = (;horizontal = 8, vertical = 4)
    polynomialorder = (;horizontal = 5, vertical = 5)
    overintegrationorder = (;horizontal = 1, vertical = 1)
    
    gridA = DiscontinuousSpectralElementGrid(ΩA, nelem, polynomialorder)
    gridO = DiscontinuousSpectralElementGrid(ΩO, nelem, polynomialorder)

    # Numerics-specific options
    numerics = (NFfirstorder = CentralNumericalFluxFirstOrder(), NFsecondorder = PenaltyNumFluxDiffusive(), overint_params = (overintegrationorder, polynomialorder) ) #, NFsecondorder = CentralNumericalFluxSecondOrder() )#PenaltyNumFluxDiffusive() )#, overint_params = (overintegrationorder, polynomialorder) ) 

    # Timestepping
    Δt_ = calculate_dt(gridA, wavespeed = u_max*(ΩA.radius), diffusivity = maximum([κᵃʰ, κᵃᶻ]) ) 
    
    t_time, end_time = ( 0  , 20Δt_ )

    # Collect spatial info, timestepping, balance law and DGmodel for the two components
    boundary_mask( xc, yc, zc ) = @. ( xc^2 + yc^2 + zc^2 )^0.5 ≈ planet_radius(param_set)
    
    # 1. Atmos component
    
    ## Prop atmos functions to override defaults
    atmos_structure(λ, ϕ, r) = FT(30)#30.0 + 10.0 * cos(ϕ) * sin(5λ)
    atmos_θⁱⁿⁱᵗ(npt, el, x, y, z) = atmos_structure( lon(x,y,z), lat(x,y,z), rad(x,y,z) )                # Set atmosphere initial state function
    #atmos_θ_shadowflux(θᵃ, θᵒ, npt, el, xc, yc, zc) = FT(0.0)
    atmos_θ_shadowflux(θᵃ, θᵒ, npt, el, xc, yc, zc) = is_surface(xc,yc,zc) ? (1.0 / τ_airsea) * (θᵃ - θᵒ) : 0.0 # Set atmosphere shadow boundary flux function
    atmos_calc_kappa_diff(_...) = κᵃʰ, κᵃʰ, κᵃᶻ               # Set atmos diffusion coeffs
    atmos_source_θ(θᵃ, npt, el, xc, yc, zc, θᵒ) = FT(0.0)     # Set atmos source!
    atmos_get_penalty_tau(_...) = FT(3.0 * 0.0)               # Set penalty term tau (for debugging)

    ## Set atmos advective velocity (constant in time) and convert to Cartesian
    uˡᵒⁿ(λ, ϕ, r) = u_max * r * cos(ϕ)
    atmos_uⁱⁿⁱᵗ(npt, el, x, y, z) = (     0 * r̂(x,y,z) 
                                        + 0 * ϕ̂(x,y,z)
                                        + uˡᵒⁿ(lon(x,y,z), lat(x,y,z), rad(x,y,z)) * λ̂(x,y,z) ) 

    ## Collect atmos props
    bl_propA = prop_defaults()
    bl_propA = (;bl_propA..., 
                init_theta = atmos_θⁱⁿⁱᵗ, 
                theta_shadow_boundary_flux = atmos_θ_shadowflux, 
                calc_kappa_diff = atmos_calc_kappa_diff,
                source_theta = atmos_source_θ,
                get_penalty_tau = atmos_get_penalty_tau,
                coupling_lambda = coupling_lambda,
                init_u = atmos_uⁱⁿⁱᵗ
                )
    ## Setup atmos component model
    mA = CplModel(;
        grid = gridA,
        equations = CplMainBL(
            bl_propA,
            (CoupledPrimaryBoundary(), ExteriorBoundary()),
            param_set,
        ),
        boundary_z = boundary_mask,
        nsteps = nstepsA,
        dt = Δt_/ nstepsA,
        timestepper = LSRK54CarpenterKennedy,
        numerics...,
    )

    # 2. Ocean component
    ## Prop ocean functions to override defaults
    tropical_heating_1(λ, ϕ, r) = 30.0 + 10.0 * cos(ϕ) * sin(5λ)
    tropical_heating_2(λ, ϕ, r) = 30.0 + 10.0 * cos(ϕ) + 1 * sin(5λ) * cos(ϕ)
    tropical_heating(λ, ϕ, r) = tropical_heating_1(λ, ϕ, r)
    ocean_θⁱⁿⁱᵗ(npt, el, x, y, z) = tropical_heating( lon(x,y,z), lat(x,y,z), rad(x,y,z) )                    # Set ocean initial state function
    ocean_calc_kappa_diff(_...) = κᵒʰ, κᵒʰ, κᵒᶻ               # Set ocean diffusion coeffs
    ocean_source_θ(θᵃ, npt, el, xc, yc, zc, θᵒ) = FT(0.0)     # Set ocean source!
    ocean_get_penalty_tau(_...) = FT(0.15 * 0.0)               # Set penalty term tau (for debugging)
    ocean_uⁱⁿⁱᵗ(xc, yc, zc, npt, el) = SVector(FT(0.0), FT(0.0), FT(0.0)) # Set ocean advective velocity

    ## Collect ocean props
    bl_propO = prop_defaults()
    bl_propO = (;bl_propO..., 
                init_theta = ocean_θⁱⁿⁱᵗ, 
                calc_kappa_diff = ocean_calc_kappa_diff,
                source_theta = ocean_source_θ,
                get_penalty_tau = ocean_get_penalty_tau,
                coupling_lambda = coupling_lambda,
                init_u = ocean_uⁱⁿⁱᵗ,
                )  
    ## Setup ocean component model
    mO = CplModel(;
        grid = gridO,
        equations = CplMainBL(
            bl_propO,
            (ExteriorBoundary(), CoupledSecondaryBoundary()),
            param_set,
        ),
        boundary_z = boundary_mask,
        nsteps = nstepsO,
        dt = Δt_ / nstepsO,
        timestepper = LSRK54CarpenterKennedy,
        numerics...,
    )

    # Create a Coupler State object for holding import/export fields.
    # Try using Dict here - not sure if that will be OK with GPU
    coupler = CplState()
    coupler_register!(coupler, :Ocean_SST, deepcopy(mO.state.θ[mO.boundary]), mO.grid, DateTime(0), u"°C")
    coupler_register!(coupler, :Atmos_MeanAirSeaθFlux, deepcopy(mA.state.F_accum[mA.boundary]), mA.grid, DateTime(0), u"°C")
    

    # Instantiate a coupled timestepper that steps forward the components and
    # implements mapings between components export bondary states and
    # other components imports.

    compA = (pre_step = preatmos, component_model = mA, post_step = postatmos)
    compO = (pre_step = preocean, component_model = mO, post_step = postocean)
    component_list = (atmosphere = compA, ocean = compO)
    cpl_solver = CplSolver(
        component_list = component_list,
        coupler = coupler,
        coupling_dt = Δt_,
        t0 = 0.0,
    )
    
    # For now applying callbacks only to atmos.
    callbacks = (
        ExponentialFiltering(),
        VTKOutput((
            iteration = string(5Δt_)*"ssecs" ,
            overdir ="output",
            overwrite = true,
            number_sample_points = 0
            )...,),   
    )

    simulation = (;
        coupled_odesolver = cpl_solver,
        odesolver = cpl_solver.component_list.atmosphere.component_model.odesolver,
        state = cpl_solver.component_list.atmosphere.component_model.state,
        dgmodel =  cpl_solver.component_list.atmosphere.component_model.discretization,
        callbacks = callbacks,
        simtime = (t_time, end_time),
        name = "Coupler_UnitTest_ctrl",
        )

    return simulation
end


function run(cpl_solver, numberofsteps, cbvector)
    # Run the model
    solve!(
        nothing,
        cpl_solver;
        numberofsteps = numberofsteps,
        callbacks = cbvector,
    )
end

simulation = main(Float64);
nsteps = Int(simulation.simtime[2] / simulation.coupled_odesolver.dt)
cbvector = create_callbacks(simulation, simulation.odesolver)
println("Initialized. Running...")
@time run(simulation.coupled_odesolver, nsteps, cbvector)
