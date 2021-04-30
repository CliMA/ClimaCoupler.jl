# # Diffusion Equation in a Vertical Column
# ## Setup
# - A vertical diffusive column split over two Cartesian components (atmos & ocean columns)
# - Boundary conditions used to set diffusive flux in other components have different vertical discretizations and timesteps
# - The coupler maps fields between the two components

# ## Equations
# The RHS evaluation follows:
# ```math
# \frac{\partial \theta}{\partial t} = - \nabla \cdot [\kappa(\phi_{init}) \nabla \phi]
# ```
# where
#  -  $\theta$ is the tracer (e.g. potential temperature)
#  -  $\kappa$ is the diffusivity tensor

# # Import packages
using ClimateMachine
using MPI
using Statistics

## To couple
using CouplerMachine
using Unitful, Dates

## To create meshes & grids
using ClimateMachine.Ocean.Domains
using ClimateMachine.Grids
import ClimateMachine.DGMethods.NumericalFluxes: NumericalFluxSecondOrder

## To setup some callbacks
using ClimateMachine.GenericCallbacks

## To invoke timestepper
using ClimateMachine.ODESolvers
using ClimateMachine.ODESolvers: solve!
using ClimateMachine.MPIStateArrays: weightedsum

ClimateMachine.init()
const FT = Float64;

if !(:CplTestingBL in names(Main))
    include("CplTestingBL.jl") # allows re-run of script without restarting julia
end
using .CplTestingBL

# # Set simulation parameters

couple_dt = 3600.0 # timestep at which coupled components sync
nstepsA = 10 # atmos steps per coupled timestep
nstepsO = 5 # ocean steps per coupled timestep

##  Haney like relaxation time scale and a length scale (Haney, 1971).
##  Air-sea exchange vigor is governed by length/time-scale.
const τ_airsea = FT(60 * 86400)
const L_airsea = FT(500)
const λ_airsea = FT(L_airsea / τ_airsea)
function coupling_lambda()
    return (λ_airsea)
end;

##  Background atmos and ocean diffusivities
const κᵃʰ = FT(1e4) * 0.0
const κᵃᶻ = FT(1e-1)
const κᵒʰ = FT(1e3) * 0.0
const κᵒᶻ = FT(1e-4)

# # Set up coupled model
# Define component models and initialize the coupler

function main(::Type{FT}) where {FT}
    ## Domain
    Np = 4
    ΩA = RectangularDomain(
        Ne = (10, 10, 5),
        Np = Np,
        x = (0, 1e6),
        y = (0, 1e6),
        z = (0, 1e5),
        periodicity = (true, true, false),
    )
    ΩO = RectangularDomain(
        Ne = (10, 10, 4),
        Np = Np,
        x = (0, 1e6),
        y = (0, 1e6),
        z = (-4e3, 0),
        periodicity = (true, true, false),
    )

    ## Grid
    btags = ((0,0),(0,0),(1,2))
    gridA = DiscontinuousSpectralElementGrid(ΩA; boundary_tags = btags)
    gridO = DiscontinuousSpectralElementGrid(ΩO; boundary_tags = btags)

    ## Numerics-specific options
    numerics = (NFsecondorder = CplTestingBL.PenaltyNumFluxDiffusive(),)

    ## Callbacks (TODO)
    callbacks = ()

    ## Collect spatial info, timestepping, balance law and DGmodel for the two components

    ## 1. Atmos component
    mA = CplModel(;
        grid = gridA,
        equations = CplTestBL(
            bl_propA,
            (CoupledPrimaryBoundary(), ExteriorBoundary()),
        ),
        nsteps = nstepsA,
        dt = couple_dt / nstepsA,
        numerics...,
    )

    ## 2. Ocean component
    mO = CplModel(;
        grid = gridO,
        equations = CplTestBL(
            bl_propO,
            (ExteriorBoundary(), CoupledSecondaryBoundary()),
        ),
        nsteps = nstepsO,
        dt = couple_dt / nstepsO,
        numerics...,
    )
#+
# Create the coupler object for holding import/export fields and performs mappings
# and instantiate the coupled timestepper:
    coupler = CplState()
    register_cpl_field!(coupler, :Ocean_SST, deepcopy(mO.state.θ[mO.boundary]), mO.grid, DateTime(0), u"°C")
    register_cpl_field!(coupler, :Atmos_MeanAirSeaθFlux, deepcopy(mA.state.F_accum[mA.boundary]), mA.grid, DateTime(0), u"°C")

    compA = (pre_step = preatmos, component_model = mA, post_step = postatmos)
    compO = (pre_step = preocean, component_model = mO, post_step = postocean)
    component_list = (atmosphere = compA, ocean = compO)
    cpl_solver = CplSolver(
        component_list = component_list,
        coupler = coupler,
        coupling_dt = couple_dt,
        t0 = 0.0,
    )

    return cpl_solver, callbacks
end

function run(cpl_solver, numberofsteps, cbvector)
    solve!(
        nothing,
        cpl_solver;
        numberofsteps = numberofsteps,
        callbacks = cbvector,
    )
end

# # Define `pre_step` and `post_step` functions
# Each component model must define `pre_step` and `post_step` functions.
# In the `pre_step`, a component imports necessary boundary state and flux data from the coupler.
# In the `post_step`, a component exports boundary data to the coupler to be later received by other components.

function get_components(csolver)
    mA = csolver.component_list.atmosphere.component_model
    mO = csolver.component_list.ocean.component_model
    return mA, mO
end

function preatmos(csolver)
    mA, mO = get_components(csolver)
    
    ## Set boundary SST used in atmos to SST of ocean surface at start of coupling cycle.
    mA.discretization.state_auxiliary.θ_secondary[mA.boundary] .= 
        CouplerMachine.get(csolver.coupler, :Ocean_SST, mA.grid, DateTime(0), u"°C")
    ## Set atmos boundary flux accumulator to 0.
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
    mA, mO = get_components(csolver)

    ## Pass atmos exports to "coupler" namespace
    ## 1. Save mean θ flux at the Atmos boundary during the coupling period
    CouplerMachine.put!(csolver.coupler, :Atmos_MeanAirSeaθFlux, mA.state.F_accum[mA.boundary] ./ csolver.dt,
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
    mA, mO = get_components(csolver)

    ## Set mean air-sea theta flux
    mO.discretization.state_auxiliary.F_prescribed[mO.boundary] .= 
        CouplerMachine.get(csolver.coupler, :Atmos_MeanAirSeaθFlux, mO.grid, DateTime(0), u"°C")
    ## Set ocean boundary flux accumulator to 0. (this isn't used)
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
    mA, mO = get_components(csolver)
    @info(
        "postocean",
        time = csolver.t + csolver.dt,
        ocean_θ_surface_max = maximum(mO.state.θ[mO.boundary]),
        ocean_θ_surface_min = maximum(mO.state.θ[mO.boundary]),
    )

    ## Pass ocean exports to "coupler" namespace
    ##  1. Ocean SST (value of θ at z=0)
    CouplerMachine.put!(csolver.coupler, :Ocean_SST, mO.state.θ[mO.boundary], mO.grid, DateTime(0), u"°C")
end

# # Specify balance law

## Set atmosphere initial state function
function atmos_init_theta(xc, yc, zc, npt, el)
    return 30.0
end
## Set atmosphere shadow boundary flux function
function atmos_theta_shadow_boundary_flux(θᵃ, θᵒ, npt, el, xc, yc, zc)
    if zc == 0.0
        tflux = (1.0 / τ_airsea) * (θᵃ - θᵒ)
    else
        tflux = 0.0
    end
    return tflux
end
## Set atmsophere diffusion coeffs
function atmos_calc_kappa_diff(_...)
    return κᵃʰ, κᵃʰ, κᵃᶻ
end
## Set atmos source!
function atmos_source_theta(θᵃ, npt, el, xc, yc, zc, θᵒ)
    tsource = 0.0
    if zc == 0.0
        ## tsource = -(1. / τ_airsea)*( θᵃ-θᵒ )
    end
    return tsource
end
## Set penalty term tau (for debugging)
function atmos_get_penalty_tau(_...)
    return FT(3.0 * 0.0)
end
## Create atmos component
bl_propA = CplTestingBL.prop_defaults()

bl_propA = (;bl_propA..., init_theta = atmos_init_theta, 
            theta_shadow_boundary_flux = atmos_theta_shadow_boundary_flux)
bl_propA = (bl_propA..., init_theta = atmos_init_theta)
bl_propA =
    (bl_propA..., theta_shadow_boundary_flux = atmos_theta_shadow_boundary_flux)
bl_propA = (bl_propA..., calc_kappa_diff = atmos_calc_kappa_diff)
bl_propA = (bl_propA..., source_theta = atmos_source_theta)
bl_propA = (bl_propA..., get_penalty_tau = atmos_get_penalty_tau)
bl_propA = (bl_propA..., coupling_lambda = coupling_lambda)

## Set initial temperature profile
function ocean_init_theta(xc, yc, zc, npt, el)
    return 20.0
end
## Set boundary source imported from atmos
function ocean_source_theta(θ, npt, el, xc, yc, zc, air_sea_flux_import)
    sval = 0.0
    if zc == 0.0
        ## sval=air_sea_flux_import
    end
    return sval
end
## Set ocean diffusion coeffs
function ocean_calc_kappa_diff(_...)
    ## return κᵒʰ,κᵒʰ,κᵒᶻ*FT(100.)
    return κᵒʰ, κᵒʰ, κᵒᶻ # m^2 s^-1
end
## Set penalty term tau (for debugging)
function ocean_get_penalty_tau(_...)
    return FT(0.15 * 0.0)
end
## Create ocean component
bl_propO = CplTestingBL.prop_defaults()
bl_propO = (bl_propO..., init_theta = ocean_init_theta)
bl_propO = (bl_propO..., source_theta = ocean_source_theta)
bl_propO = (bl_propO..., calc_kappa_diff = ocean_calc_kappa_diff)
bl_propO = (bl_propO..., get_penalty_tau = ocean_get_penalty_tau)
bl_propO = (bl_propO..., coupling_lambda = coupling_lambda)

# # Run simulation
simulation, cbvector = main(Float64);
nsteps = 10
println("Initialized. Running...")
run(simulation, nsteps, cbvector)
