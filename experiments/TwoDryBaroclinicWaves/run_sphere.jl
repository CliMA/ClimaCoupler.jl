#!/usr/bin/env julia --project
include("utilities/boilerplate.jl")
using CouplerMachine

########
# Set up parameters and initial conditions
########
include("parameters_initialconditions.jl")

FT = Float64

########
# Set up domain
########
# A: low atmos level (troposphere)
domainA = SphericalShell(
    radius = FT(planet_radius(param_set)), 
    height = FT(30e3),
)
# B: high atmos level (stratosphere)
domainB = SphericalShell(
    radius = planet_radius(param_set) + FT(30e3),
    height = (30e3),
)

gridA = DiscretizedDomain(
    domainA;
    elements = (vertical = 5, horizontal = 6),
    polynomial_order = (vertical = 3, horizontal = 6),
    overintegration_order = (vertical = 0, horizontal = 0),
)

gridB = DiscretizedDomain(
    domainB;
    elements = (vertical = 5, horizontal = 6),
    polynomial_order = (vertical = 3, horizontal = 6),
    overintegration_order = (vertical = 0, horizontal = 0),
)

########
# Set up model physics
########

ref_state = DryReferenceState(DecayingTemperatureProfile{FT}(param_set, FT(290), FT(220), FT(8e3)))

# total energy
eos     = TotalEnergy(γ = 1 / (1 - parameters.κ))   

physics = Physics(
    orientation = SphericalOrientation(),
    ref_state   = ref_state,
    eos         = eos,
    lhs         = (
        ESDGNonLinearAdvection(eos = eos),
        PressureDivergence(eos = eos),
    ),
    sources     = (
        DeepShellCoriolis{FT}(Ω = parameters.Ω),
    ),
)

linear_eos = linearize(physics.eos)
linear_physics = Physics(
    orientation = physics.orientation,
    ref_state   = physics.ref_state,
    eos         = linear_eos,
    lhs         = (
        ESDGLinearAdvection(),
        PressureDivergence(eos = linear_eos),
    ),
    sources     = (
        ThinShellGravityFromPotential(),
    ),
)

########
# Set up model
########

modelA = DryAtmosModel(
    physics = physics,
    boundary_conditions = (ExteriorBoundary(), CoupledSecondaryBoundary()),
    initial_conditions = (ρ = ρ₀ᶜᵃʳᵗ, ρu = ρu⃗₀ᶜᵃʳᵗ, ρe = ρeᶜᵃʳᵗ),
    numerics = (
        flux = RusanovNumericalFlux(),
    ),
    parameters = parameters,
)

linear_modelA = DryAtmosLinearModel(
    physics = linear_physics,
    boundary_conditions = modelA.boundary_conditions,
    initial_conditions = nothing,
    numerics = (
        flux = modelA.numerics.flux,
        direction = VerticalDirection()
    ),
    parameters = modelA.parameters,
)

modelB = DryAtmosModel(
    physics = physics,
    boundary_conditions = (CoupledPrimaryBoundary(), ExteriorBoundary()),
    initial_conditions = (ρ = ρ₀ᶜᵃʳᵗ, ρu = ρu⃗₀ᶜᵃʳᵗ, ρe = ρeᶜᵃʳᵗ),
    numerics = (
        flux = RusanovNumericalFlux(),
    ),
    parameters = parameters,
)

linear_modelB = DryAtmosLinearModel(
    physics = linear_physics,
    boundary_conditions = modelB.boundary_conditions,
    initial_conditions = nothing,
    numerics = (
        flux = modelB.numerics.flux,
        direction = VerticalDirection()
    ),
    parameters = modelB.parameters,
)
########
# Set up time steppers (could be done automatically in simulation)
########
# determine the time step construction
# element_size = (domain_height / numelem_vert)
# acoustic_speed = soundspeed_air(param_set, FT(330))
dx = minimum( [min_node_distance(gridA.numerical) , min_node_distance(gridB.numerical) ])
cfl = 3
Δt = cfl * dx / 330.0
start_time = 0
end_time = Δt*10#30 * 24 * 3600
method = ARK2GiraldoKellyConstantinescu
callbacks = (
  Info(),
  CFL(),
  VTKState(
    iteration = Int(floor(6*3600/Δt)), 
    filepath = "./out/"),
)

########
# Set up simulation
########
nstepsA = 1
nstepsB = 1

epss = sqrt(eps(FT))
boundary_mask( param_set, xc, yc, zc ) = @. abs(( xc^2 + yc^2 + zc^2 )^0.5 - planet_radius(param_set) - FT(30e3)) < epss
    
simA = CplSimulation(
    (modelA, linear_modelA,);
    grid = gridA,
    timestepper = (method = method, timestep = Δt / nstepsA),
    time        = (start = start_time, finish = end_time),
    nsteps      = nstepsA,
    boundary_z = boundary_mask,
    callbacks   = callbacks,
)
simB = CplSimulation(
    (modelB, linear_modelB,);
    grid = gridB,
    timestepper = (method = method, timestep = Δt / nstepsB),
    time        = (start = start_time, finish = end_time),
    nsteps      = nstepsB,
    boundary_z = boundary_mask,
    callbacks   = callbacks,
)

## Create a Coupler State object for holding imort/export fields.
coupler = CplState()
coupler_register!(coupler, :EnergyA, deepcopy(simA.state.ρe[simA.boundary]), simA.grid, DateTime(0), u"J") # value on top of domainA for calculating upward flux into domainB
coupler_register!(coupler, :EnergyFluxB, deepcopy(simB.state.F_ρe_accum[simB.boundary]), simB.grid, DateTime(0), u"J") # downward flux

compA = (pre_step = preA, component_model = simA, post_step = postA)
compB = (pre_step = preB, component_model = simB, post_step = postB)
component_list = (domainA = compA, domainB = compB)


cpl_solver = CplSolver(
    component_list = component_list,
    coupler = coupler,
    coupling_dt = Δt,
    t0 = FT(start_time),
)

########
# Run the simulation
########
numberofsteps = Int( round((end_time - start_time) / Δt))
evolve!(cpl_solver, numberofsteps)

nothing