# This test runs a small CMIP simulation twice.
#
# - The first time the simulation is run for two coupler steps
# - The second time the simulation is run for two coupler steps, but restarting
#   from the first simulation
#
# Since the caches are not read from the restart file, the results will not be
# bit-wise identical to running the simulation without restarting. This test
# only checks that a simulation can be restarted from the state only without
# erroring.
#
# This test exercises the Oceananigans checkpointing path: checkpoint_model_state
# writes a JLD2 file via OC.checkpoint, and restart_model_state! restores it
# via OC.set!.

import ClimaComms
ClimaComms.@import_required_backends
import ClimaUtilities.OutputPathGenerator: maybe_wait_filesystem
import YAML
import Logging
using Test

# Uncomment the following for cleaner output (but more difficult debugging)
# Logging.disable_logging(Logging.Warn)

include("compare.jl")
include("../code_loading.jl")

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
config_file = joinpath(@__DIR__, "restart_state_only_cmip.yml")
config_dict = Input.get_coupler_config_dict(config_file)

# Two coupler steps
two_steps = deepcopy(config_dict)

two_steps["dt_atmos"] = "120secs"
two_steps["dt_cpl"] = "360secs"
two_steps["dt_land"] = "120secs"
two_steps["dt_ocean"] = "360secs"
two_steps["dt_seaice"] = "360secs"
two_steps["t_end"] = "720secs"
two_steps["checkpoint_dt"] = "720secs"
two_steps["coupler_output_dir"] = tmpdir
two_steps["job_id"] = "two_steps_cmip"
two_steps["save_cache"] = false

println("Simulating two steps")
cs_two_steps = setup_and_run(two_steps)

@testset "Cache files should not exist" begin
    checkpoints_dir = joinpath(tmpdir, two_steps["job_id"], "checkpoints")
    checkpoint_files = filter!(isfile, readdir(checkpoints_dir, join = true))
    are_not_cache_files = .!occursin.("cache", basename.(checkpoint_files))
    @test all(are_not_cache_files)
end

# Check that we can pick up a simulation by providing t_restart and restart_dir
println("Simulating two steps, options from command line")
two_steps_reading = deepcopy(two_steps)

two_steps_reading["t_end"] = "1080secs"
two_steps_reading["detect_restart_files"] = true
two_steps_reading["restart_dir"] = cs_two_steps.dir_paths.checkpoints_dir
two_steps_reading["restart_t"] = 720
two_steps_reading["restart_cache"] = false
two_steps_reading["job_id"] = "two_steps_cmip_reading"
two_steps_reading["save_cache"] = false
# calling setup_and_run bypasses get_coupler_config_dict which modifies t_start
Input.update_t_start_for_restarts!(two_steps_reading)

cs_two_steps_reading = setup_and_run(two_steps_reading)
@testset "CMIP restarts (state only)" begin
    @test cs_two_steps_reading.tspan[1] == cs_two_steps.tspan[2]
end
