# Import necessary packages
import ClimaComms
import ClimaCore as CC
import ClimaAtmos as CA
import ClimaLand as CL

# Choose device to run on (CPU or GPU)
context = ClimaComms.context()
device = ClimaComms.device(context)

# Choose float type
FT = Float64

# Parse optional config file - overwrites default values
!isnothing(config_file) && parse_inputs(config_file)

# Set up exchange space - requires component models to regrid to/from this space
space = create_sphere(FT, context) # including topography and land/sea mask

# Set up timestepping information
tspan = (; t_start = Float64(0), t_end = Float64(60 * 60 * 24)) # simulation length 1 day
dt_cpl = (400)

# Set up output directory
output_dir = "output"

# Set up parameters - need to be unified for calibration
# Default parameters, can be overwritten by TOML files or user inputs
parameters = init_params()

# Initialize component model simulations - maybe from restart
# Simulation instead of model allows each component to step/solve independently
atmos = CA.atmos_init(FT, space, output_dir)
land = CL.LandSimulation(FT, space, output_dir) # or stub_init()
ocean = ocean_init(FT, space, output_dir) # or stub_init()
ice = ice_init(FT, space, output_dir) # or stub_init()

coupled_model = CoupledModel(atmos, ocean, land, ice)

# Initialize coupled simulation
coupled_simulation = CoupledSimulation(atmos, ocean, land, seaice, tspan, dt_cpl)

# Run very short/coarse simulation to display which input files are required
check_artifacts(coupled_simulation)

# Run simulation
solve_coupler!(coupled_simulation)

# Postprocessing (optional - will bring in many deps) - maybe in extension
# Dispatch on config ID
postprocess(coupled_simulation)

# Q: How can we communicate to user which input files are needed?
# One option: run short, coarse simulation first to display which files are requires

"""
CoupledSimulation
- models (atmos, ocean, land, seaice)
- clock (tspan, dt_cpl, time/iteration)
- space (exchange space)
- land/sea mask (think about resolution/space)
- callback
- diagnostics/output writers
- turbulent flux partition
- coupler exchange fields
- logging
- config ID
- job ID
- initialized/running flags
- preprocessing
"""
# Similar simulation for components

# Optional user inputs:
# - callbacks
# - logging (verbose mode?) - option for MPI pids to write outputs to files
# - config ID, job ID - useful for running ensembles and distinguishing runs
# - initialized/running flags - e.g. to keep track of consistency between state/cache when starting and stopping
# - preprocessing (gapfilling, etc.) - could include restarts
