# This test runs a small CMIP simulation four times.
#
# - The first time the simulation is run for four steps
# - The second time the simulation is run for two steps
# - The third time the simulation is run for two steps, but restarting from the
#   second simulation
#
# After all these simulations are run, we compare the first and last runs. They
# should be bit-wise identical.
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

include("compare_cmip.jl")
include("../CMIP/code_loading.jl")

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
config_file = joinpath(@__DIR__, "restart_cmip.yml")
config_dict = Input.get_coupler_config_dict(config_file)

# Four steps
four_steps = deepcopy(config_dict)

four_steps["dt_cpl"] = "360secs"
four_steps["dt_ocean"] = "360secs"
four_steps["dt_seaice"] = "360secs"
four_steps["t_end"] = "1440secs"
four_steps["coupler_output_dir"] = tmpdir
four_steps["checkpoint_dt"] = "1440secs"
four_steps["job_id"] = "four_steps"

println("Simulating four steps")
cs_four_steps = setup_and_run(four_steps)

# Check that we can pick up a simulation by providing t_restart and restart_dir
println("Simulating four steps, options from command line")
four_steps_reading = deepcopy(four_steps)

four_steps_reading["t_end"] = "1800secs"
four_steps_reading["detect_restart_files"] = true
four_steps_reading["restart_dir"] = cs_four_steps.dir_paths.checkpoints_dir
four_steps_reading["restart_t"] = 1440
four_steps_reading["job_id"] = "four_steps_reading"
Input.update_t_start_for_restarts!(four_steps_reading)

cs_four_steps_reading = setup_and_run(four_steps_reading)
@testset "CMIP restarts (state and cache)" begin
    @test cs_four_steps_reading.tspan[1] == cs_four_steps.tspan[2]
end

# Two steps + two steps (2 × 360 s = 720 s each half)
two_steps = deepcopy(config_dict)

two_steps["dt_cpl"] = "360secs"
two_steps["dt_ocean"] = "360secs"
two_steps["dt_seaice"] = "360secs"
two_steps["t_end"] = "720secs"
two_steps["coupler_output_dir"] = tmpdir
# restart_cmip.yml sets dt_nogw/dt_ogw to 360 s so checkpoint_dt (360 s here, 1440 s in
# four_steps) is an integer multiple of the GW callback periods (ClimaAtmos checks
# checkpoint_dt / dt_nogw and checkpoint_dt / dt_ogw).
two_steps["checkpoint_dt"] = "360secs"
two_steps["job_id"] = "two_steps"

# Copying since setup_and_run changes its content
println("Simulating two steps")
cs_two_steps1 = setup_and_run(two_steps)

println("Restarting from checkpoint, initialization only")
# Construct a restarted CoupledSimulation at t = 720s, but do not advance it.
restart_init = deepcopy(two_steps)
restart_init["t_end"] = "720secs" # equal to t_start after update_t_start_for_restarts!
restart_init["detect_restart_files"] = true
restart_init["restart_dir"] = cs_two_steps1.dir_paths.checkpoints_dir
restart_init["restart_t"] = 720
restart_init["restart_cache"] = true
restart_init["job_id"] = "two_steps_restart_init_only"
Input.update_t_start_for_restarts!(restart_init)
cs_two_steps_restart_init = Interfacer.CoupledSimulation(restart_init)

@testset "Restart initialization matches checkpointed state" begin
    @test cs_two_steps_restart_init.tspan[1] == cs_two_steps1.tspan[2]

    # Compare prognostic states after restart initialization (including coupler flux init)
    # to the end-of-run state from the pre-restart segment (the checkpointed step).
    @test compare(
        cs_two_steps1.model_sims.atmos_sim.integrator.u,
        cs_two_steps_restart_init.model_sims.atmos_sim.integrator.u,
    )
    @test compare(
        cs_two_steps1.model_sims.land_sim.integrator.u,
        cs_two_steps_restart_init.model_sims.land_sim.integrator.u,
    )
    # Ice/ocean: compare prognostic Fields only. Walking the full model hits
    # reconstructed grid/cache internals (not restored from JLD2) and is not a
    # meaningful restart check. `compare_cmip.jl` compares Field data via `parent`.
    @test compare(
        OC.prognostic_fields(cs_two_steps1.model_sims.ice_sim.ice.model),
        OC.prognostic_fields(cs_two_steps_restart_init.model_sims.ice_sim.ice.model),
    )
    @test compare(
        OC.prognostic_fields(cs_two_steps1.model_sims.ocean_sim.ocean.model),
        OC.prognostic_fields(cs_two_steps_restart_init.model_sims.ocean_sim.ocean.model),
    )
end
