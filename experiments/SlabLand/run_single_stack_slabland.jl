#!/usr/bin/env julia --project
include("utilities/boilerplate.jl")
include("balance_law_interface_slabland.jl")

using CouplerMachine

########
# Set up parameters and initial conditions
########
include("parameters_initialconditions_slabland.jl")


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
domainLand = SingleStack(
    xmax = FT(parameters.xmax), 
    ymax = FT(parameters.ymax), 
    zmin = -FT(parameters.zmax),
    zmax = FT(0),
)
# B: high atmos level (atmosphere)
domainAtmos = SingleStack(
    xmax = FT(parameters.xmax), 
    ymax = FT(parameters.ymax), 
    zmin = FT(0),
    zmax = FT(parameters.zmax),
)

gridLand = DiscretizedDomain(
    domainLand;
    elements = (vertical = 20, horizontal = 1),
    polynomial_order = (vertical = 5, horizontal = 5),
    overintegration_order = (vertical = 1, horizontal = 1),
)

gridAtmos = DiscretizedDomain(
    domainAtmos;
    elements = (vertical = 20, horizontal = 1),
    polynomial_order = (vertical = 5, horizontal = 5),
    overintegration_order = (vertical = 1, horizontal = 1),
)

########
# Set up model physics
######## 



physicsLand = Physics(
    orientation = FlatOrientation(),
    #diffusion   = ConstantViscosity( FT(0), FT(0), FT(1e3) ),
    eos         = DryIdealGas{Float64}(R = parameters.R_d, pₒ = parameters.pₒ, γ = 1 / (1 - parameters.κ)),
)

physicsAtmos = Physics(
    orientation = FlatOrientation(),
    diffusion   = ConstantViscosity( FT(0), FT(0), FT(1e-3) ),
    eos         = DryIdealGas{Float64}(R = parameters.R_d, pₒ = parameters.pₒ, γ = 1 / (1 - parameters.κ)),
)




########
# Set up boundary conditions
########

bcsLand = (
    bottom =  (T_sfc = Insulating(),),
    top =     (T_sfc = CoupledSecondaryBoundary(),),
)
bcsAtmos = (
    bottom = (ρu = Impenetrable(FreeSlip()), ρθ = CoupledPrimaryBoundary()),
    top =    (ρu = Impenetrable(FreeSlip()), ρθ = Insulating()),
)

########
# Set up model
########

modelLand = SlabLandModelSetup(
    physics = physicsLand,
    boundary_conditions = bcsLand,
    initial_conditions = (T_sfc = T_sfc₀,),
    numerics = (flux = RoeNumericalFlux(),flux_second_order = PenaltyNumFluxDiffusive(), direction = EveryDirection()),
    parameters = parameters,
)

modelAtmos = ModelSetup(
    physics = physicsAtmos,
    boundary_conditions = bcsAtmos,
    initial_conditions =(ρ = ρ₀B, ρu = ρu⃗₀B, ρθ = ρθ₀B),
    numerics = (flux = RoeNumericalFlux(),flux_second_order = PenaltyNumFluxDiffusive(), direction = EveryDirection()),
    parameters = parameters,
)

########
# Set up time steppers (could be done automatically in simulation)
########
Δt  = min_node_distance(gridAtmos.numerical) / parameters.cₛ * 0.25
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
    
simLand = CplSimulation(
    modelLand;
    grid = gridLand,
    timestepper = (method = method, timestep = Δt / nstepsA),
    time        = (start = start_time, finish = end_time),
    nsteps      = nstepsA,
    boundary_z = boundary_mask,
    callbacks   = callbacks,
)
simAtmos = CplSimulation(
    modelAtmos;
    grid = gridAtmos,
    timestepper = (method = method, timestep = Δt / nstepsB),
    time        = (start = start_time, finish = end_time),
    nsteps      = nstepsB,
    boundary_z = boundary_mask,
    callbacks   = callbacks,
)

## Create a Coupler State object for holding imort/export fields.
coupler = CplState()
register_cpl_field!(coupler, :EnergyLand, deepcopy(simLand.state.T_sfc[simLand.boundary]), simLand.grid, DateTime(0), u"J") # value on top of domainA for calculating upward flux into domainB
register_cpl_field!(coupler, :EnergyFluxAtmos, deepcopy(simAtmos.state.F_ρθ_accum[simAtmos.boundary]), simAtmos.grid, DateTime(0), u"J") # downward flux

compLand = (pre_step = preLand, component_model = simLand, post_step = postLand)
compAtmos = (pre_step = preAtmos, component_model = simAtmos, post_step = postAtmos)
component_list = (domainLand = compLand, domainAtmos = compAtmos)


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