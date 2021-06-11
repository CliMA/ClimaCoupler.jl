#!/usr/bin/env julia --project
using CouplerMachine
using Unitful

include("utilities/boilerplate.jl")
include("numerics/timestepper_abstractions.jl")
include("balance_law_interface_slabocean.jl")

########
# Set up parameters and initial conditions
########
include("parameters_initialconditions_slabocean.jl")

FT = Float64

# Collects tracer valus/fluxes for conservation checks
mutable struct LogThatFlux{F}
    A::F
    B::F
end

########
# Set up domain
########

# A: low atmos level (ocean)
domainOcean = SphericalShell(
    radius = parameters.a,
    height = parameters.H,
)

# B: high atmos level (atmosphere)
domainAtmos = SphericalShell(
    radius = parameters.a,
    height = parameters.H,
)

gridOcean = DiscretizedDomain(
    domainOcean;
    elements = (vertical = 10, horizontal = 32),
    polynomial_order = (vertical = 2, horizontal = 2),
    overintegration_order = (vertical = 0, horizontal = 0),
)

gridAtmos = DiscretizedDomain(
    domainAtmos;
    elements = (vertical = 10, horizontal = 32),
    polynomial_order = (vertical = 2, horizontal = 2),
    overintegration_order = (vertical = 0, horizontal = 0),
)

########
# Set up model physics
######## 
FT = Float64

ref_state = DryReferenceState(
  DecayingTemperatureProfile{FT}(parameters, FT(290), FT(220), FT(8e3))
)

physicsOcean = Physics(
    orientation = SphericalOrientation(),
    parameters = parameters,
)

physicsAtmos = Physics(
    orientation = SphericalOrientation(),
    ref_state   = ref_state,
    eos         = MoistIdealGas(),
    lhs         = (
        NonlinearAdvection{(:ρ, :ρu, :ρe)}(),
        PressureDivergence(),
    ),
    sources     = (
        DeepShellCoriolis(),
        FluctuationGravity(),
        ZeroMomentMicrophysics(),
        HeldSuarezForcing(held_suarez_parameters),
        FluxAccumulator(),
    ),
    parameters = parameters,
)

linear_physicsAtmos = Physics(
    orientation = physicsAtmos.orientation,
    ref_state   = physicsAtmos.ref_state,
    eos         = physicsAtmos.eos,
    lhs         = (
        LinearAdvection{(:ρ, :ρu, :ρe)}(),
        LinearPressureDivergence(),
    ),
    sources     = (
        FluctuationGravity(),
    ),
    parameters = parameters,
)

# ########
# # Set up boundary conditions
# ########

oceanBC = ( Insulating(), CoupledSecondaryAtmosModelBC() )

atmosBC = (
        DefaultBC() ,  # (Insulating, FreeSlip) 
        CoupledPrimarySlabOceanBC( parameters.Cₑ, parameters.Cₗ),
        )

########
# Set up model
########

modelOcean = SlabOceanModelSetup(
    physics = physicsOcean,
    boundary_conditions = oceanBC,
    initial_conditions = (T_sfc = T_sfc₀,),
    numerics = (flux = CentralNumericalFluxFirstOrder(), flux_second_order = CentralNumericalFluxSecondOrder(), direction = EveryDirection()),
    #parameters = parameters,
)

modelAtmos = DryAtmosModel(
    physics = physicsAtmos,
    boundary_conditions = atmosBC, 
    initial_conditions = (ρ = ρ₀ᶜᵃʳᵗ, ρu = ρu⃗₀ᶜᵃʳᵗ, ρe = ρeᶜᵃʳᵗ, ρq = ρqᶜᵃʳᵗ),
    numerics = (flux = LMARSNumericalFlux(),),
)

linear_modelAtmos = DryAtmosModel(
    physics = linear_physicsAtmos,
    boundary_conditions = (DefaultBC(), DefaultBC()), # (Insulating, FreeSlip), (Insulating, FreeSlip) 
    initial_conditions = modelAtmos.initial_conditions,
    numerics = (flux = RefanovFlux(),),
)

########
# Set up time steppers (could be done automatically in simulation)
########
dx = min_node_distance(gridAtmos.numerical)
cfl = 5 # 13 for 10 days, 7.5 for 200+ days
Δt = cfl * dx / 330.0

total_steps = 1000
start_time = 0
end_time = Δt * total_steps#30 * 24 * 3600

methodOcean = SSPRK22Heuns
methodAtmos = IMEX() 

callbacks = (
  Info(),
  CFL(),
  VTKState(
    iteration = Int(floor(6*3600/Δt)), 
    filepath = "./output/SphereSlabOcean"),
)

########
# Set up simulation
########
nstepsOcean = 1
nstepsAtmos = 1

epss = sqrt(eps(FT))
boundary_mask( param_set, xc, yc, zc ) = @. abs(( xc^2 + yc^2 + zc^2 )^0.5 - param_set.a - param_set.H) < epss
    
simOcean = CplSimulation(
    modelOcean;
    grid = gridOcean,
    timestepper = (method = methodOcean, timestep = Δt / nstepsOcean),
    time        = (start = start_time, finish = end_time),
    nsteps      = nstepsOcean,
    boundary_z = boundary_mask,
    callbacks   = callbacks,
)
simAtmos = CplSimulation(
    (Explicit(modelAtmos), Implicit(linear_modelAtmos),);
    grid = gridAtmos,
    timestepper = (method = methodAtmos, timestep = Δt / nstepsAtmos),
    time        = (start = start_time, finish = end_time),
    nsteps      = nstepsAtmos,
    boundary_z = boundary_mask,
    callbacks   = callbacks,
)

## Create a Coupler State object for holding import/export fields.
coupler = CplState()
register_cpl_field!(coupler, :SeaSurfaceTemerature, deepcopy(simOcean.state.T_sfc[simOcean.boundary]), simOcean.grid, DateTime(0), u"K") # value on top of domainA for calculating upward flux into domainB
register_cpl_field!(coupler, :EnergyFluxAtmos, deepcopy(simAtmos.state.F_ρe_accum[simAtmos.boundary]), simAtmos.grid, DateTime(0), u"J") # downward flux

compOcean = (pre_step = preOcean, component_model = simOcean, post_step = postOcean)
compAtmos = (pre_step = preAtmos, component_model = simAtmos, post_step = postAtmos)
component_list = (domainOcean = compOcean, domainAtmos = compAtmos)

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
evolve!(cpl_solver, total_steps)

########
# Check conservation
########
using Plots
fluxA = cpl_solver.fluxlog.A
fluxB = cpl_solver.fluxlog.B

fluxT = fluxA .+ fluxB 
time = collect(1:1:total_steps)
rel_error = [ ((fluxT .- fluxT[1]) / fluxT[1]) ]
plot(time .* cpl_solver.dt, rel_error, ylabel = "rel. error = (fluxT - fluxT[1]) / fluxT[1]", xlabel = "time (s)")

using Statistics
plot(time .* cpl_solver.dt, [(fluxA .- fluxA[1]) (fluxB .- fluxB[1])],  label = ["Ocean Energy" "Atmos Energy"], xlabel = "time (s)", ylabel = "J / m2")
