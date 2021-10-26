#!/usr/bin/env julia --project
include("utilities/boilerplate.jl")
using CouplerMachine

########
# Set up parameters and initial conditions
########
include("parameters_initialconditions.jl")

FT = Float64

# Collects tracer values/fluxes for conservation checks
mutable struct LogThatFlux{F}
    A::F
    B::F
end

########
# Set up domain
########
# A: low atmos level (troposphere)
domainA = SphericalShell(radius = FT(1), height = FT(0.2))
# B: high atmos level (stratosphere)
domainB = SphericalShell(radius = domainA.radius + domainA.height, height = domainA.height)

gridA = DiscretizedDomain(
    domainA;
    elements = (vertical = 1, horizontal = 8),
    polynomial_order = (vertical = 0, horizontal = 3),
    overintegration_order = (vertical = 1, horizontal = 1),
)

gridB = DiscretizedDomain(
    domainB;
    elements = (vertical = 1, horizontal = 8),
    polynomial_order = (vertical = 0, horizontal = 3),
    overintegration_order = (vertical = 1, horizontal = 1),
)

########
# Set up model physics
######## 
physics = Physics(
    orientation = SphericalOrientation(),
    advection = NonLinearAdvection(),
    coriolis = DeepShellCoriolis{Float64}(Ω = parameters.Ω),
    gravity = Buoyancy{Float64}(α = parameters.α, g = 0.0),
    eos = BarotropicFluid{Float64}(ρₒ = parameters.ρₒ, cₛ = parameters.cₛ),
)

########
# Set up boundary conditions
########
bcsA = (
    bottom = (ρu = Impenetrable(FreeSlip()), ρθ = Insulating()),
    # top = (ρu = Impenetrable(FreeSlip()), ρθ = Insulating()), # for non-coupled testing
    top = (ρu = Impenetrable(FreeSlip()), ρθ = CoupledSecondaryBoundary()),
)
bcsB = (
    bottom = (ρu = Impenetrable(FreeSlip()), ρθ = CoupledPrimaryBoundary()),
    # bottom = (ρu = Impenetrable(FreeSlip()), ρθ = Insulating()), #for non-coupled testing
    top = (ρu = Impenetrable(FreeSlip()), ρθ = Insulating()),
)

########
# Set up model
########
modelA = ModelSetup(
    physics = physics,
    boundary_conditions = bcsA,
    initial_conditions = (ρ = ρ₀ᶜᵃʳᵗ, ρu = ρu⃗₀ᶜᵃʳᵗ, ρθ = ρθ₀ᶜᵃʳᵗA),
    numerics = (flux = RoeNumericalFlux(), flux_second_order = PenaltyNumFluxDiffusive()),
    parameters = parameters,
)

modelB = ModelSetup(
    physics = physics,
    boundary_conditions = bcsB,
    initial_conditions = (ρ = ρ₀ᶜᵃʳᵗ, ρu = ρu⃗₀ᶜᵃʳᵗ, ρθ = ρθ₀ᶜᵃʳᵗB),
    numerics = (flux = RoeNumericalFlux(), flux_second_order = PenaltyNumFluxDiffusive()),
    parameters = parameters,
)

########
# Set up time steppers (could be done automatically in simulation)
########
Δt = min_node_distance(gridA.numerical) / parameters.cₛ * 0.25
total_steps = 10000
start_time = 0
end_time = Δt * total_steps#30 * 24 * 3600
method = SSPRK22Heuns
callbacksA = (
    Info(),
    CFL(),
    #   VTKState(
    #     iteration = 250,#Int(floor(6*3600/Δt)), 
    #     filepath = "./out/A/"),
)
callbacksB = (
    Info(),
    CFL(),
    #   VTKState(
    #     iteration = 250,#Int(floor(6*3600/Δt)), 
    #     filepath = "./out/B/"),
)

########
# Set up simulation
########
nstepsA = 1
nstepsB = 1

epss = sqrt(eps(FT))
boundary_mask(param_set, xc, yc, zc) = @. abs((xc^2 + yc^2 + zc^2)^0.5 - domainB.radius) < epss

simA = CplSimulation(
    modelA;
    grid = gridA,
    timestepper = (method = method, timestep = Δt / nstepsA),
    time = (start = start_time, finish = end_time),
    nsteps = nstepsA,
    boundary_z = boundary_mask,
    callbacks = callbacksA,
)
simB = CplSimulation(
    modelB;
    grid = gridB,
    timestepper = (method = method, timestep = Δt / nstepsB),
    time = (start = start_time, finish = end_time),
    nsteps = nstepsB,
    boundary_z = boundary_mask,
    callbacks = callbacksB,
)

## Create a Coupler State object for holding import/export fields.
coupler = CplState()
coupler_register!(coupler, :EnergyA, deepcopy(simA.state.ρθ[simA.boundary]), simA.grid, DateTime(0), u"J") # value on top of domainA for calculating upward flux into domainB
coupler_register!(coupler, :EnergyFluxB, deepcopy(simB.state.F_ρθ_accum[simB.boundary]), simB.grid, DateTime(0), u"J") # downward flux

compA = (pre_step = preA, component_model = simA, post_step = postA)
compB = (pre_step = preB, component_model = simB, post_step = postB)
component_list = (domainA = compA, domainB = compB)
cpl_solver = CplSolver(
    component_list = component_list,
    coupler = coupler,
    coupling_dt = Δt,
    t0 = FT(start_time),
    fluxlog = LogThatFlux(zeros(total_steps), zeros(total_steps)),
)

########
# Run the simulation
########
evolve!(cpl_solver, total_steps)

########
# Check conservation
########
using Plots
fluxA = cpl_solver.fluxlog.A
fluxB = cpl_solver.fluxlog.B
fluxT = fluxA .+ fluxB
time = collect(1:1:total_steps)
rel_error = [((fluxT .- fluxT[1]) / fluxT[1])]
plot(time .* cpl_solver.dt, rel_error)
#plot(time .* cpl_solver.dt, [fluxA .- fluxA[1] fluxB .- fluxB[1]])
