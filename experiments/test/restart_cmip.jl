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

four_steps["t_end"] = "960secs"
four_steps["coupler_output_dir"] = tmpdir
four_steps["checkpoint_dt"] = "960secs"
four_steps["job_id"] = "four_steps"

println("Simulating four steps")
cs_four_steps = setup_and_run(four_steps)

# Check that we can pick up a simulation by providing t_restart and restart_dir
println("Simulating four steps, options from command line")
four_steps_reading = deepcopy(four_steps)

four_steps_reading["t_end"] = "1200secs"
four_steps_reading["detect_restart_files"] = true
four_steps_reading["restart_dir"] = cs_four_steps.dir_paths.checkpoints_dir
four_steps_reading["restart_t"] = 960
four_steps_reading["job_id"] = "four_steps_reading"
Input.update_t_start_for_restarts!(four_steps_reading)

cs_four_steps_reading = setup_and_run(four_steps_reading)
@testset "CMIP restarts (state and cache)" begin
    @test cs_four_steps_reading.tspan[1] == cs_four_steps.tspan[2]
end

# Now, two steps plus one
two_steps = deepcopy(config_dict)

two_steps["t_end"] = "480secs"
two_steps["coupler_output_dir"] = tmpdir
two_steps["checkpoint_dt"] = "480secs"
two_steps["job_id"] = "two_steps"

# Copying since setup_and_run changes its content
println("Simulating two steps")
cs_two_steps1 = setup_and_run(two_steps)

println("Reading and simulating last two steps")
# Two additional steps
two_steps["t_end"] = "960secs"
two_steps["detect_restart_files"] = true
Input.update_t_start_for_restarts!(two_steps)
cs_two_steps2 = setup_and_run(two_steps)

@testset "Restarts" begin
    # We put cs_four_steps.fields in a NamedTuple so that we can start the recursion in compare
    @test compare(
        (; coupler_fields = cs_four_steps.fields),
        (; coupler_fields = cs_two_steps2.fields),
    )

    @test compare(
        cs_four_steps.model_sims.atmos_sim.integrator.u,
        cs_two_steps2.model_sims.atmos_sim.integrator.u,
    )

    @test compare(
        cs_four_steps.model_sims.atmos_sim.integrator.p,
        cs_two_steps2.model_sims.atmos_sim.integrator.p,
        ignore = [
            :walltime_estimate,                 # Stateful
            :output_dir,                        # Changes across runs
            :scratch,                           # Irrelevant
            :ghost_buffer,                      # Irrelevant
            :hyperdiffusion_ghost_buffer,       # Irrelevant
            :data_handler,                      # Stateful
            :face_clear_sw_direct_flux_dn,      # Not filled by RRTGMP
            :face_sw_direct_flux_dn,            # Not filled by RRTGMP
            :rc,                                # CUDA internal object
        ],
    )

    @test compare(
        cs_four_steps.model_sims.land_sim.integrator.u,
        cs_two_steps2.model_sims.land_sim.integrator.u,
    )
    @test compare(
        cs_four_steps.model_sims.land_sim.integrator.p,
        cs_two_steps2.model_sims.land_sim.integrator.p,
        ignore = [:dss_buffer_3d, :dss_buffer_2d, :rc],
    )

    @test compare(
        cs_four_steps.model_sims.ice_sim.ice.model,
        cs_two_steps2.model_sims.ice_sim.ice.model,
        ignore = [:clock, :parent, :ptr],
    )

    @test compare(
        cs_four_steps.model_sims.ocean_sim.ocean.model,
        cs_two_steps2.model_sims.ocean_sim.ocean.model,
        ignore = [:clock, :parent, :ptr],
    )
end
