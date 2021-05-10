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
totalsteps = 1000 # total simulation steps

#  Background atmos and ocean horizontal and vertical diffusivities
const κᵃʰ = FT(0.0) #FT(1e4) 
const κᵃᶻ = FT(1e-1)
const κᵒʰ = FT(0.0)
const κᵒᶻ = FT(1e-4)

# Collects tracer valus/fluxes for conservation checks
mutable struct LogThatFlux{F}
    A::F
    O::F
end

# Coupling coefficient
τ_airsea,  L_airsea = ( FT(60 * 86400), FT(500) )
coupling_lambda() = (L_airsea / τ_airsea)

# Max. advective velocity in m/s
const u_max = FT(1e-5 * 6371000) 

function main(::Type{FT}) where {FT}
    
    # Domain (approximately emulating the dimensions of one of the faces of the cubed sphere)
    # Grid
    nelem = (;horizontal = 8, vertical = 4)
    polynomialorder = (;horizontal = 5, vertical = 5)
    overintegrationorder = (;horizontal = 1, vertical = 1)

    ΩA = RectangularDomain(;
        Ne = (nelem.horizontal, nelem.horizontal, nelem.vertical),
        Np = polynomialorder.horizontal,
        x = (0, 1e7),
        y = (0, 1e7),
        z = (0, 4e3),
        periodicity = (true, true, false),
    )
    ΩO = RectangularDomain(;
        Ne = (nelem.horizontal, nelem.horizontal, nelem.vertical),
        Np = polynomialorder.horizontal,
        x = (0, 1e7),
        y = (0, 1e7),
        z = (-4e3, 0),
        periodicity = (true, true, false),
    )

    ## Grid
    btags = ((0,0),(0,0),(1,2))
    gridA = DiscontinuousSpectralElementGrid(ΩA; boundary_tags = btags)
    gridO = DiscontinuousSpectralElementGrid(ΩO; boundary_tags = btags)

    # Numerics-specific options
    numerics = (NFfirstorder = CentralNumericalFluxFirstOrder(), NFsecondorder = PenaltyNumFluxDiffusive(), overint_params = (overintegrationorder, polynomialorder) ) #, NFsecondorder = CentralNumericalFluxSecondOrder() )#PenaltyNumFluxDiffusive() )#, overint_params = (overintegrationorder, polynomialorder) ) 

    # Time step 
    Δt_ = calculate_dt(gridA, wavespeed = u_max, diffusivity = maximum([κᵃʰ, κᵃᶻ]) ) 
    t_time, end_time = ( 0  , totalsteps * Δt_ )    
    
    # 1. Atmos component
    
    ## Prop atmos functions to override defaults
    atmos_θⁱⁿⁱᵗ(npt, el, x, y, z) = FT(30)             # Set atmosphere initial state function
    atmos_calc_diff_flux(∇θ, npt, el, x, y, z) = Diagonal(@SVector([κᵃʰ, κᵃʰ, κᵃᶻ])) * ∇θ               # Set atmos diffusion coeffs
    atmos_get_penalty_tau(_...) = FT(3.0 * 0.0)               # Set penalty term tau (for debugging)
    #atmos_θ_shadowflux(θᵃ, θᵒ, npt, el, xc, yc, zc) = is_surface(xc,yc,zc) ? (1.0 / τ_airsea) * (θᵃ - θᵒ) : 0.0 # Set atmosphere shadow boundary flux function
    #atmos_source_θ(θᵃ, npt, el, xc, yc, zc, θᵒ) = FT(0.0)     # Set atmos source!

    ## Set atmos advective velocity (constant in time) and convert to Cartesian
    epss = sqrt(eps(FT))
    #atmos_uⁱⁿⁱᵗ(npt, el, x, y, z) = (z  < ΩA.z[1]+epss) || (z > ΩA.z[2]-epss) ? SVector(FT(0.0), FT(0.0), FT(0.0)) : SVector(u_max, FT(0.0), FT(0.0)) # constant u, but 0 at boundaries
    atmos_uⁱⁿⁱᵗ(npt, el, x, y, z) = SVector(FT(0.0), FT(0.0), FT(0.0))  #(z  < 0.00001) || (z > 3999.999) ? SVector(FT(0.0), FT(0.0), FT(0.0)) : SVector(u_max, FT(0.0), FT(0.0)) # constant u, but 0 at boundaries
    
    ## Collect atmos props
    bl_propA = prop_defaults()
    bl_propA = (;bl_propA..., 
                init_theta = atmos_θⁱⁿⁱᵗ, 
                calc_diff_flux = atmos_calc_diff_flux,
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
            FlatOrientation(),
        ),
        nsteps = nstepsA,
        dt = Δt_/ nstepsA,
        timestepper = LSRK54CarpenterKennedy,
        numerics...,
    )

    # 2. Ocean component
    ## Prop ocean functions to override defaults
    ocean_θⁱⁿⁱᵗ(npt, el, x, y, z) =  FT(0.0) #+ FT(10.0) * cos(π*x/1e6)                # Set ocean initial state function
    ocean_calc_diff_flux(∇θ, npt, el, x, y, z) = Diagonal(@SVector([κᵒʰ, κᵒʰ, κᵒᶻ])) * ∇θ              # Set ocean diffusion coeffs
    ocean_get_penalty_tau(_...) = FT(0.15 * 0.0)               # Set penalty term tau (for debugging)
    ocean_uⁱⁿⁱᵗ(xc, yc, zc, npt, el) = SVector(FT(0.0), FT(0.0), FT(0.0)) # Set ocean advective velocity
    #ocean_source_θ(θᵃ, npt, el, xc, yc, zc, θᵒ) = FT(0.0)     # Set ocean source!
    
    ## Collect ocean props
    bl_propO = prop_defaults()
    bl_propO = (;bl_propO..., 
                init_theta = ocean_θⁱⁿⁱᵗ, 
                calc_diff_flux = ocean_calc_diff_flux,
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
            FlatOrientation(),
        ),
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
        # ExponentialFiltering(),
        # VTKOutput((
        #     iteration = string(100Δt_)*"ssecs" ,
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
        name = "Coupler_BoxUnitTest_ctrl",
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
nsteps = totalsteps
cbvector = create_callbacks(simulation, simulation.odesolver)
println("Initialized. Running...")
@time run(simulation.coupled_odesolver, nsteps, cbvector)

# Check conservation
fluxA = simulation.coupled_odesolver.fluxlog.A
fluxO = simulation.coupled_odesolver.fluxlog.O
fluxT = fluxA .+ fluxO
#plot(collect(1:1:totalsteps),[fluxA .- fluxA[1], fluxO .- fluxO[1], fluxT .- fluxT[1] ])

using Plots
time = collect(1:1:totalsteps)[1:10:end]
rel_error = [ ((fluxT .- fluxT[1]) / fluxT[1])[1:10:end] ]
plot(time,rel_error)
