# This test runs a small AMIP simulation three times.
#
# - The first time the simulation is run for three steps
# - The second time the simulation is run for two steps
# - The third time the simulation is run for one step, but restarting from the
#   second time
#
# After all these simulations are run, we compare the first and last runs. They
# should be bit-wise identical.
#
# The content of the simulation is not the most important, but it helps if has a
# all the complexity possible.

import ClimaComms
ClimaComms.@import_required_backends
import ClimaUtilities.OutputPathGenerator: maybe_wait_filesystem
import YAML
import Logging
using Test

# Uncomment the following for cleaner output (but more difficult debugging)
# Logging.disable_logging(Logging.Warn)

include("compare.jl")
comms_ctx = ClimaComms.context()
@info "Context: $(comms_ctx)"
ClimaComms.init(comms_ctx)

# Make sure that all MPI processes agree on the output_loc
tmpdir = ClimaComms.iamroot(comms_ctx) ? mktempdir(pwd()) : ""
tmpdir = ClimaComms.bcast(comms_ctx, tmpdir)
# Sometimes the shared filesystem doesn't work properly and the folder is not
# synced across MPI processes. Let's add an additional check here.
maybe_wait_filesystem(ClimaComms.context(), tmpdir)

base_config_file = joinpath(@__DIR__, "amip_test.yml")
base_config = YAML.load_file(base_config_file)

# Three steps
three_steps = deepcopy(base_config)

three_steps["dt"] = "180secs"
three_steps["dt_cpl"] = "180secs"
three_steps["t_end"] = "540secs"
three_steps["dt_rad"] = "180secs"
three_steps["coupler_output_dir"] = joinpath(tmpdir, "three_steps")
three_steps["checkpoint_dt"] = "180secs"

YAML.write_file(joinpath(tmpdir, "three_steps.yml"), three_steps)

push!(ARGS, "--config_file", joinpath(tmpdir, "three_steps.yml"))
push!(ARGS, "--job_id", "three_steps")
module TestThree
include("../run_amip.jl")
end

empty!(ARGS)
println("Simulating two steps")

# Now, two steps plus one
two_steps = deepcopy(base_config)

two_steps["dt"] = "180secs"
two_steps["dt_cpl"] = "180secs"
two_steps["t_end"] = "360secs"
two_steps["dt_rad"] = "180secs"
two_steps["coupler_output_dir"] = joinpath(tmpdir, "two_steps")
two_steps["checkpoint_dt"] = "180secs"

YAML.write_file(joinpath(tmpdir, "two_steps.yml"), two_steps)

push!(ARGS, "--config_file", joinpath(tmpdir, "two_steps.yml"))
push!(ARGS, "--job_id", "two_steps")
module TestTwo1
include("../run_amip.jl")
end

empty!(ARGS)
println("Reading and simulating last step")

# Restart (just re-run from the same folder)
restarted = YAML.load_file(joinpath(tmpdir, "two_steps.yml"))
restarted["t_end"] = "540secs"
YAML.write_file(joinpath(tmpdir, "restarted.yml"), restarted)
push!(ARGS, "--config_file", joinpath(tmpdir, "restarted.yml"))
push!(ARGS, "--job_id", "two_steps")
module TestTwo2
include("../run_amip.jl")
end

@test compare(TestThree.cs.model_sims.atmos_sim.integrator.u, TestTwo2.cs.model_sims.atmos_sim.integrator.u)

@test compare(
    TestThree.cs.model_sims.atmos_sim.integrator.p,
    TestTwo2.cs.model_sims.atmos_sim.integrator.p,
    ignore = [
        :walltime_estimate,
        :output_dir,
        :scratch,
        :ghost_buffer,
        :hyperdiffusion_ghost_buffer,
        :data_handler,
        :face_clear_sw_direct_flux_dn,
        :face_sw_direct_flux_dn,
        :rc,
    ],
)


@testset "Restarts" begin
    @test compare(TestThree.cs.model_sims.ice_sim.integrator.u, TestTwo2.cs.model_sims.ice_sim.integrator.u)

    @test compare(TestThree.cs.model_sims.land_sim.integrator.u, TestTwo2.cs.model_sims.land_sim.integrator.u)
    @test compare(
        TestThree.cs.model_sims.land_sim.integrator.p,
        TestTwo2.cs.model_sims.land_sim.integrator.p,
        ignore = [:dss_buffer_3d, :dss_buffer_2d, :rc],
    )

    @test compare(
        TestThree.cs.model_sims.land_sim.integrator.p,
        TestTwo2.cs.model_sims.land_sim.integrator.p;
        ignore = [:dss_buffer_3d, :dss_buffer_2d, :rc],
    )

    # Ignoring SST_timevaryinginput because it contains closures (which should be
    # reinitialized correctly). We have to remove it from the type, otherwise
    # compare will not work
    function delete(nt::NamedTuple, fieldnames...)
        return (; filter(p -> !(first(p) in fieldnames), collect(pairs(nt)))...)
    end

    ocean_cache_three = delete(TestThree.cs.model_sims.ocean_sim.cache, :SST_timevaryinginput)
    ocean_cache_two2 = delete(TestTwo2.cs.model_sims.ocean_sim.cache, :SST_timevaryinginput)

    @test compare(ocean_cache_three, ocean_cache_two2)
end
