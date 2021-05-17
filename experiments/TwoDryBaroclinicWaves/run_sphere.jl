#!/usr/bin/env julia --project
include("../interface/utilities/boilerplate.jl")

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

model = DryAtmosModel(
    physics = physics,
    boundary_conditions = (5, 6),
    initial_conditions = (ρ = ρ₀ᶜᵃʳᵗ, ρu = ρu⃗₀ᶜᵃʳᵗ, ρe = ρeᶜᵃʳᵗ),
    numerics = (
        flux = RusanovNumericalFlux(),
    ),
    parameters = parameters,
)

linear_model = DryAtmosLinearModel(
    physics = linear_physics,
    boundary_conditions = model.boundary_conditions,
    initial_conditions = nothing,
    numerics = (
        flux = model.numerics.flux,
        direction = VerticalDirection()
    ),
    parameters = model.parameters,
)

########
# Set up time steppers (could be done automatically in simulation)
########
# determine the time step construction
# element_size = (domain_height / numelem_vert)
# acoustic_speed = soundspeed_air(param_set, FT(330))
dx = min_node_distance(grid.numerical)
cfl = 3
Δt = cfl * dx / 330.0
start_time = 0
end_time = 30 * 24 * 3600
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
CplSimulation(model::Tuple;
    grid,
    timestepper,
    time,
    boundary_z = nothing,
    nsteps,
    callbacks)
simA = CplSimulation(
    (model, linear_model,);
    grid = gridA,
    timestepper = (method = method, timestep = Δt),
    time        = (start = start_time, finish = end_time),
    callbacks   = callbacks,
)
simB = CplSimulation(
    (model, linear_model,);
    grid = gridB,
    timestepper = (method = method, timestep = Δt),
    time        = (start = start_time, finish = end_time),
    callbacks   = callbacks,
)

## Create a Coupler State object for holding imort/export fields.
coupler = CplState()
register_cpl_field!(coupler, :EnergyFluxA, deepcopy(simA.state.ρe_accum[simA.boundary]), simA.grid, DateTime(0), u"J") # upward flux
register_cpl_field!(coupler, :EnergyFluxB, deepcopy(simB.state.ρe_accum[simB.boundary]), simB.grid, DateTime(0), u"J") # downward flux

compA = (pre_step = preA, component_model = simA, post_step = postA)
compB = (pre_step = preB, component_model = simB, post_step = postB)
component_list = (domainA = compA, domainB = compB)


cpl_solver = CplSolver(
    component_list = component_list,
    coupler = coupler,
    coupling_dt = couple_dt,
    t0 = 0.0,
)

########
# Run the simulation
########
evolve!(simulation)

nothing