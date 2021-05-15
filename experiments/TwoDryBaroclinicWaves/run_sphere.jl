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
    domain;
    elements = (vertical = 5, horizontal = 6),
    polynomial_order = (vertical = 3, horizontal = 6),
    overintegration_order = (vertical = 0, horizontal = 0),
)

gridB = DiscretizedDomain(
    domain;
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
method = 
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
simulation = Simulation(
    (model, linear_model,);
    grid = grid,
    timestepper = (method = method, timestep = Δt),
    time        = (start = start_time, finish = end_time),
    callbacks   = callbacks,
)

########
# Run the simulation
########
initialize!(simulation)
evolve!(simulation)

nothing