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

# Shared functions
include("domains.jl")
include("abstractions.jl")
include("callbacks.jl")

# Main balance law and its components
include("CplMainBL.jl")

FT = Float64
nstepsA = 1 # steps per coupling cycle (atmos)
nstepsO = 1 # steps per coupling cycle (ocean)
totalsteps = 400 # total simulation coupled cycle steps

# Background atmos and ocean horizontal and vertical diffusivities
const κᵃʰ, κᵃᶻ = ( FT(0.0) , FT(1e-1) )
const κᵒʰ, κᵒᶻ = ( FT(0.0) , FT(1e-4) )
const νᵃʰ, νᵃᶻ = ( FT(0.0) , FT(1e-1) )
const νᵒʰ, νᵒᶻ = ( FT(0.0) , FT(1e-4) )

# Collects tracer valus/fluxes for conservation checks
mutable struct LogThatFlux{F}
    A::F
    O::F
end

# Coupling coefficient
const τ_airsea,  L_airsea = ( FT(60 * 86400), FT(500) ) #( FT(1), FT(0) )#( FT(60 * 86400), FT(500) ) 
coupling_lambda() = (L_airsea / τ_airsea)



# Max. advective velocity in radians
const u_max = FT(1e-5) 

"""
function main(::Type{FT}) where {FT}
    sets up the experiment and initializes the run

    FT: Float64
"""
function main(::Type{FT}) where {FT}
    
    # Domain
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
    Δt_ = calculate_dt(gridA, wavespeed = u_max*(ΩA.radius), diffusivity = maximum([κᵃʰ, κᵃᶻ, κᵒʰ, κᵒᶻ]), dif_direction = VerticalDirection() ) # 
    
    t_time, end_time = ( 0  , totalsteps * Δt_ )

    # Collect spatial info, timestepping, balance law and DGmodel for the two components
    epss = sqrt(eps(FT))
    boundary_mask( param_set, xc, yc, zc ) = @. abs(( xc^2 + yc^2 + zc^2 )^0.5 - planet_radius(param_set)) < epss
    
    # 1. Atmos component
    
    ## Prop atmos functions to override defaults
    get_atmos_diff() = [κᵃʰ, κᵃʰ, κᵃᶻ]
    get_atmos_visc() = [νᵃʰ, νᵃʰ, νᵃᶻ]
    atmos_structure(λ, ϕ, r) = FT(30) + FT(10.0) * r / ΩA.radius + 10.0 * cos(ϕ) * sin(5λ)
    atmos_θⁱⁿⁱᵗ(npt, el, x, y, z) = atmos_structure( lon(x,y,z), lat(x,y,z), rad(x,y,z) )                # Set atmosphere initial state function
    atmos_calc_diff_flux(∇θ, npt, el, x, y, z, κᵃʰ, κᵃᶻ) = - (
                κᵃʰ * ∇θ'* λ̂_quick(x,y,z) * λ̂_quick(x,y,z) +  # zonal
                κᵃʰ * ∇θ'* ϕ̂_quick(x,y,z) * ϕ̂_quick(x,y,z) +  # meridional
                κᵃᶻ * ∇θ'* r̂_quick(x,y,z) * r̂_quick(x,y,z) ) # vertical

    atmos_calc_visc_flux(∇u, npt, el, x, y, z, νᵃʰ, νᵃᶻ) = [ atmos_calc_diff_flux(∇u[:,1], npt, el, x, y, z, νᵃʰ, νᵃᶻ) atmos_calc_diff_flux(∇u[:,2], npt, el, x, y, z, νᵃʰ, νᵃᶻ) atmos_calc_diff_flux(∇u[:,3], npt, el, x, y, z, νᵃʰ, νᵃᶻ) ]

    atmos_get_penalty_tau(_...) = FT(3.0 * 0.0)               # Set penalty term tau (for debugging)
    #atmos_θ_shadowflux(θᵃ, θᵒ, npt, el, xc, yc, zc) = is_surface(xc,yc,zc) ? (1.0 / τ_airsea) * (θᵃ - θᵒ) : 0.0 # Set atmosphere shadow boundary flux function
    #atmos_source_θ(θᵃ, npt, el, xc, yc, zc, θᵒ) = FT(0.0)     # Set atmos source!

    ## Set atmos advective velocity (constant in time) and convert to Cartesian
    uˡᵒⁿ(λ, ϕ, r) = u_max * r * cos(ϕ)

    atmos_uⁱⁿⁱᵗ_(npt, el, x, y, z) = (    FT(0) * r̂(x,y,z) 
                                        + FT(0) * ϕ̂(x,y,z)
                                        + uˡᵒⁿ(lon(x,y,z), lat(x,y,z), rad(x,y,z)) * λ̂(x,y,z) ) 
    
    # atmos_uⁱⁿⁱᵗ_(npt, el, x, y, z) = (    FT(0) * r̂(x,y,z) 
    #                                     + FT(0) * ϕ̂(x,y,z)
    #                                     + FT(0) * λ̂(x,y,z) ) 

    atmos_uⁱⁿⁱᵗ(npt, el, x, y, z) = (rad(x,y,z)  < ΩA.radius + epss) || (rad(x,y,z) > ΩA.radius + ΩA.height - epss) ? SVector(FT(0.0), FT(0.0), FT(0.0)) : atmos_uⁱⁿⁱᵗ_(npt, el, x, y, z) # constant (in rads) u, but 0 at boundaries

    ## Collect atmos props
    bl_propA = prop_defaults()
    bl_propA = (;bl_propA..., 
                init_theta = atmos_θⁱⁿⁱᵗ, 
                get_diff = get_atmos_diff,
                get_visc = get_atmos_visc,
                calc_diff_flux = atmos_calc_diff_flux,
                calc_visc_flux = atmos_calc_visc_flux,
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
            SphericalOrientation(),
        ),
        boundary_z = boundary_mask,
        nsteps = nstepsA,
        dt = Δt_/ nstepsA,
        timestepper = LSRK54CarpenterKennedy,
        numerics...,
    )

    # 2. Ocean component
    ## Prop ocean functions to override defaults
    get_ocean_diff() = [κᵒʰ, κᵒʰ, κᵒᶻ]
    get_ocean_visc() = [νᵒʰ, νᵒʰ, νᵒᶻ]
    tropical_heating_1(λ, ϕ, r) = FT(0.0) #+ FT(10.0) * cos(ϕ) * sin(5λ)
    tropical_heating_2(λ, ϕ, r) = 30.0 + 10.0 * cos(ϕ) + 1 * sin(5λ) * cos(ϕ)
    tropical_heating(λ, ϕ, r) = tropical_heating_1(λ, ϕ, r)
    ocean_θⁱⁿⁱᵗ(npt, el, x, y, z) = tropical_heating( lon(x,y,z), lat(x,y,z), rad(x,y,z) )                    # Set ocean initial state function
    ocean_calc_diff_flux(∇θ, npt, el, x, y, z, κᵒʰ, κᵒᶻ) = - (
                κᵒʰ * ∇θ' * λ̂_quick(x,y,z) * λ̂_quick(x,y,z) +  # zonal
                κᵒʰ * ∇θ' * ϕ̂_quick(x,y,z) * ϕ̂_quick(x,y,z) +  # meridional
                κᵒᶻ * ∇θ' * r̂_quick(x,y,z) * r̂_quick(x,y,z) ) # vertical
    ocean_calc_visc_flux(∇u, npt, el, x, y, z, νᵃʰ, νᵃᶻ) = [ ocean_calc_diff_flux(∇u[:,1], npt, el, x, y, z, νᵃʰ, νᵃᶻ) ocean_calc_diff_flux(∇u[:,2], npt, el, x, y, z, νᵃʰ, νᵃᶻ) ocean_calc_diff_flux(∇u[:,3], npt, el, x, y, z, νᵃʰ, νᵃᶻ) ]
    ocean_get_penalty_tau(_...) = FT(0.15 * 0.0)               # Set penalty term tau (for debugging)
    ocean_uⁱⁿⁱᵗ(xc, yc, zc, npt, el) = SVector(FT(0.0), FT(0.0), FT(0.0)) # Set ocean advective velocity
    #ocean_source_θ(θᵃ, npt, el, xc, yc, zc, θᵒ) = FT(0.0)     # Set ocean source!
    
    ## Collect ocean props
    bl_propO = prop_defaults()
    bl_propO = (;bl_propO..., 
                init_theta = ocean_θⁱⁿⁱᵗ, 
                get_diff = get_ocean_diff,
                get_visc = get_ocean_visc,
                calc_diff_flux = ocean_calc_diff_flux,
                calc_visc_flux = ocean_calc_visc_flux,
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
            SphericalOrientation(),
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
    register_cpl_field!(coupler, :Ocean_SST, deepcopy(mO.state.θ[mO.boundary]), mO.grid, DateTime(0), u"°C")
    register_cpl_field!(coupler, :Atmos_MeanAirSeaθFlux, deepcopy(mA.state.F_accum[mA.boundary]), mA.grid, DateTime(0), u"°C")

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
        fluxlog = LogThatFlux( zeros(totalsteps) , zeros(totalsteps) )
    )
    
    # For now applying callbacks only to atmos.
    callbacks = (
        #ExponentialFiltering(),
        # VTKOutput((
        #     iteration = string(200Δt_)*"ssecs" ,
        #     overdir ="output",
        #     overwrite = true,
        #     number_sample_points = 0
        #     )...,),   
    )

    simulation = (;
        coupled_odesolver = cpl_solver,
        odesolver = cpl_solver.component_list.atmosphere.component_model.stepper,
        state = cpl_solver.component_list.atmosphere.component_model.state,
        dgmodel =  cpl_solver.component_list.atmosphere.component_model.discretization,
        callbacks = callbacks,
        simtime = (t_time, end_time),
        name = "Coupler_SphereUnitTest_ctrl",
        )

    return simulation
end

"""
function run(cpl_solver, numberofsteps, cbvector)
    runs the simulation

    cpl_solver: CplSolver
    numberofsteps: number of coupling cycles
    cbvector: vector of callbacks

"""
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
nsteps = totalsteps
cbvector = create_callbacks(simulation, simulation.odesolver)
println("Initialized. Running...")
@time run(simulation.coupled_odesolver, nsteps, cbvector)

# Check conservation
using Plots
fluxA = simulation.coupled_odesolver.fluxlog.A
fluxO = simulation.coupled_odesolver.fluxlog.O
fluxT = fluxA .+ fluxO
time = collect(1:1:totalsteps)
rel_error = [ ((fluxT .- fluxT[1]) / fluxT[1]) ]
plot(time .* simulation.coupled_odesolver.dt,rel_error)

plot(time .* simulation.coupled_odesolver.dt,[fluxA .- fluxA[1],fluxO .- fluxO[1],fluxT .- fluxT[1]])