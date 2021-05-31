#!/usr/bin/env julia --project
include("utilities/boilerplate.jl")
using CouplerMachine

########
# Set up parameters and initial conditions
########
include("parameters_initialconditions.jl")

FT = Float64


# Collects tracer valus/fluxes for conservation checks
mutable struct LogThatFlux{F}
    A::F
    B::F
end

########
# Set up domain
########

# A: low atmos level (land)
domainA = SingleStack(
    xmax = FT(parameters.xmax), 
    ymax = FT(parameters.ymax), 
    zmin = -FT(parameters.zmax),
    zmax = FT(0),
)
# B: high atmos level (atmosphere)
domainB = SingleStack(
    xmax = FT(parameters.xmax), 
    ymax = FT(parameters.ymax), 
    zmin = FT(0),
    zmax = FT(parameters.zmax),
)

gridA = DiscretizedDomain(
    domainA;
    elements = (vertical = 20, horizontal = 1),
    polynomial_order = (vertical = 5, horizontal = 5),
    overintegration_order = (vertical = 1, horizontal = 1),
)

gridB = DiscretizedDomain(
    domainB;
    elements = (vertical = 20, horizontal = 1),
    polynomial_order = (vertical = 5, horizontal = 5),
    overintegration_order = (vertical = 1, horizontal = 1),
)

########
# Set up model physics
######## 

physics = Physics(
    orientation = FlatOrientation(),
    diffusion   = ConstantViscosity( FT(0), FT(0), FT(1e3) ),
    eos         = DryIdealGas{Float64}(R = parameters.R_d, pₒ = parameters.pₒ, γ = 1 / (1 - parameters.κ)),
)


########
# Set up boundary conditions
########

bcsA = (
    bottom = (ρu = Impenetrable(FreeSlip()), ρθ = Insulating()),
    top =    (ρu = Impenetrable(FreeSlip()), ρθ = CoupledSecondaryBoundary()),
)
bcsB = (
    bottom = (ρu = Impenetrable(FreeSlip()), ρθ = CoupledPrimaryBoundary()),
    top =    (ρu = Impenetrable(FreeSlip()), ρθ = Insulating()),
)

########
# Set up model
########

modelA = ModelSetup(
    physics = physics,
    boundary_conditions = bcsA,
    initial_conditions = (ρ = ρ₀A, ρu = ρu⃗₀A, ρθ = ρθ₀A),
    numerics = (flux = RoeNumericalFlux(),flux_second_order = PenaltyNumFluxDiffusive(), direction = EveryDirection()),
    parameters = parameters,
)

modelB = ModelSetup(
    physics = physics,
    boundary_conditions = bcsB,
    initial_conditions =(ρ = ρ₀B, ρu = ρu⃗₀B, ρθ = ρθ₀B),
    numerics = (flux = RoeNumericalFlux(),flux_second_order = PenaltyNumFluxDiffusive(), direction = EveryDirection()),
    parameters = parameters,
)

########
# Set up time steppers (could be done automatically in simulation)
########
Δt  = min_node_distance(gridA.numerical) / parameters.cₛ * 0.25
total_steps = 50
start_time = 0
end_time = Δt * total_steps#30 * 24 * 3600
method = SSPRK22Heuns
callbacks = (
  Info(),
  CFL(),
  VTKState(
    iteration = Int(floor(6*3600/Δt)), 
    filepath = "./output/SingleStack"),
)

########
# Set up simulation
########
nstepsA = 1
nstepsB = 1

epss = sqrt(eps(Float64))
boundary_mask( param_set, xc, yc, zc ) = @. abs(zc) < epss
    
simA = CplSimulation(
    modelA;
    grid = gridA,
    timestepper = (method = method, timestep = Δt / nstepsA),
    time        = (start = start_time, finish = end_time),
    nsteps      = nstepsA,
    boundary_z = boundary_mask,
    callbacks   = callbacks,
)
simB = CplSimulation(
    modelB;
    grid = gridB,
    timestepper = (method = method, timestep = Δt / nstepsB),
    time        = (start = start_time, finish = end_time),
    nsteps      = nstepsB,
    boundary_z = boundary_mask,
    callbacks   = callbacks,
)

## Create a Coupler State object for holding imort/export fields.
coupler = CplState()
register_cpl_field!(coupler, :EnergyA, deepcopy(simA.state.ρθ[simA.boundary]), simA.grid, DateTime(0), u"J") # value on top of domainA for calculating upward flux into domainB
register_cpl_field!(coupler, :EnergyFluxB, deepcopy(simB.state.F_ρθ_accum[simB.boundary]), simB.grid, DateTime(0), u"J") # downward flux

compA = (pre_step = preA, component_model = simA, post_step = postA)
compB = (pre_step = preB, component_model = simB, post_step = postB)
component_list = (domainA = compA, domainB = compB)


cpl_solver = CplSolver(
    component_list = component_list,
    coupler = coupler,
    coupling_dt = Δt,
    t0 = FT(start_time),
    fluxlog = LogThatFlux( zeros(total_steps) , zeros(total_steps) )
)

########
# Run the simulation
########
numberofsteps = Int( round((end_time - start_time) / Δt))
evolve!(cpl_solver, numberofsteps)

########
# Check conservation
########
using Plots
fluxA = cpl_solver.fluxlog.A
fluxB = cpl_solver.fluxlog.B
fluxT = fluxA .+ fluxB
time = collect(1:1:total_steps)
rel_error = [ ((fluxT .- fluxT[1]) / fluxT[1]) ]
plot(time .* cpl_solver.dt, rel_error)

using Statistics
plot(time .* cpl_solver.dt, [fluxA .- fluxA[1] fluxB .- fluxB[1]] )