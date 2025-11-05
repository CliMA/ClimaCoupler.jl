# This test runs a small AMIP simulation twice times.
#
# - The first time the simulation is run for two steps
# - The second time the simulation is run for two steps, but restarting from the
#   first simulation
#
# Since the caches are not read from the restart file, the results will not be
# bit-wise identical to running the simulation without restarting. This test
# only checks that a simulation can be restarted from the state only without erroring.
#
# The content of the simulation is not the most important, but it helps if it
# has all of the complexity possible.

import ClimaComms
ClimaComms.@import_required_backends
import ClimaUtilities.OutputPathGenerator: maybe_wait_filesystem
import YAML
import Logging
using Test

# Uncomment the following for cleaner output (but more difficult debugging)
# Logging.disable_logging(Logging.Warn)

include("compare.jl")
include("../setup_run.jl")

comms_ctx = ClimaComms.context()
@info "Context: $(comms_ctx)"
ClimaComms.init(comms_ctx)

# Make sure that all MPI processes agree on the output_loc
tmpdir = ClimaComms.iamroot(comms_ctx) ? mktempdir(pwd()) : ""
tmpdir = ClimaComms.bcast(comms_ctx, tmpdir)
# Sometimes the shared filesystem doesn't work properly and the folder is not
# synced across MPI processes. Let's add an additional check here.
maybe_wait_filesystem(ClimaComms.context(), tmpdir)

# Parse the input config file as a dictionary
config_file = joinpath(@__DIR__, "amip_test.yml")
config_dict = get_coupler_config_dict(config_file)

# Four steps
two_steps = deepcopy(config_dict)

two_steps["dt"] = "180secs"
two_steps["dt_cpl"] = "180secs"
two_steps["t_end"] = "360secs"
two_steps["dt_rad"] = "180secs"
two_steps["checkpoint_dt"] = "360secs"
two_steps["coupler_output_dir"] = tmpdir
two_steps["job_id"] = "two_steps"

println("Simulating two steps")
cs_two_steps = setup_and_run(two_steps)

# Check that we can pick up a simulation by providing t_restart and restart_dir
println("Simulating two steps, options from command line")
two_steps_reading = deepcopy(two_steps)

two_steps_reading["t_end"] = "540secs"
two_steps_reading["detect_restart_files"] = true
two_steps_reading["restart_dir"] = cs_two_steps.dir_paths.checkpoints_dir
two_steps_reading["restart_t"] = 360
two_steps_reading["restart_cache"] = false
two_steps_reading["job_id"] = "two_steps_reading"

cs_two_steps_reading = setup_and_run(two_steps_reading)
@testset "Restarts from command line arguments" begin
    @test cs_two_steps_reading.tspan[1] == cs_two_steps.tspan[2]
end
