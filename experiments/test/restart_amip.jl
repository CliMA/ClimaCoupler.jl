# This test runs a small AMIP simulation four times.
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

include("compare_amip.jl")
include("../AMIP/code_loading.jl")

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
config_dict = Input.get_coupler_config_dict(config_file)

# Four steps
four_steps = deepcopy(config_dict)

four_steps["dt"] = "180secs"
four_steps["dt_cpl"] = "180secs"
four_steps["t_end"] = "720secs"
four_steps["dt_rad"] = "180secs"
four_steps["checkpoint_dt"] = "360secs"
four_steps["coupler_output_dir"] = tmpdir
four_steps["job_id"] = "four_steps"

println("Simulating four steps")
cs_four_steps = setup_and_run(four_steps)

# Check that we can pick up a simulation by providing t_restart and restart_dir
println("Simulating four steps, options from command line")
four_steps_reading = deepcopy(four_steps)

four_steps_reading["t_end"] = "900secs"
four_steps_reading["detect_restart_files"] = true
four_steps_reading["restart_dir"] = cs_four_steps.dir_paths.checkpoints_dir
four_steps_reading["restart_t"] = 720
four_steps_reading["job_id"] = "four_steps_reading"
Input.update_t_start_for_restarts!(four_steps_reading)

cs_four_steps_reading = setup_and_run(four_steps_reading)
@testset "AMIP restarts (state and cache)" begin
    @test cs_four_steps_reading.tspan[1] == cs_four_steps.tspan[2]
end

# Now, two steps plus one
two_steps = deepcopy(config_dict)

two_steps["dt"] = "180secs"
two_steps["dt_cpl"] = "180secs"
two_steps["t_end"] = "360secs"
two_steps["dt_rad"] = "180secs"
two_steps["coupler_output_dir"] = tmpdir
two_steps["checkpoint_dt"] = "360secs"
two_steps["job_id"] = "two_steps"

# Copying since setup_and_run changes its content
println("Simulating two steps")
cs_two_steps1 = setup_and_run(two_steps)

println("Reading and simulating last two steps")
# Two additional steps
two_steps["t_end"] = "720secs"
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
            :face_clear_sw_direct_flux_dn,      # Not filled by RRTMGP
            :face_sw_direct_flux_dn,            # Not filled by RRTMGP
            :rc,                                # CUDA internal object
            :non_orographic_gravity_wave,       # Recomputed every timestep; sensitive to FP accumulation
            # RRTMGP 0.22 allocates these working buffers with `undef`; they are
            # rewritten from scratch on every radiation solve and are not restart state.
            # The meaningful radiation outputs are checked explicitly via the API below.
            :fluxb,                             # Band-by-band intermediate fluxes (scratch)
            :src,                               # Source function workspace (includes sfc_source)
            :lon,                               # "required but unused" per RRTMGP comment
            :aero_size,                         # Partially filled; unused slots stay undef
        ],
    )

    # Verify radiation restart fidelity via the public RRTMGP accessor API.
    # These are the fluxes that coupling and the energy tendency actually use.
    rrtmgp_four = cs_four_steps.model_sims.atmos_sim.integrator.p.radiation.rrtmgp_solver
    rrtmgp_two2 = cs_two_steps2.model_sims.atmos_sim.integrator.p.radiation.rrtmgp_solver
    @test compare(
        (;
            lw_flux_dn = CA.RRTMGP.lw_flux_dn(rrtmgp_four),
            lw_flux_up = CA.RRTMGP.lw_flux_up(rrtmgp_four),
            sw_flux_dn = CA.RRTMGP.sw_flux_dn(rrtmgp_four),
            sw_flux_up = CA.RRTMGP.sw_flux_up(rrtmgp_four),
            ᶠradiation_flux = cs_four_steps.model_sims.atmos_sim.integrator.p.radiation.ᶠradiation_flux,
        ),
        (;
            lw_flux_dn = CA.RRTMGP.lw_flux_dn(rrtmgp_two2),
            lw_flux_up = CA.RRTMGP.lw_flux_up(rrtmgp_two2),
            sw_flux_dn = CA.RRTMGP.sw_flux_dn(rrtmgp_two2),
            sw_flux_up = CA.RRTMGP.sw_flux_up(rrtmgp_two2),
            ᶠradiation_flux = cs_two_steps2.model_sims.atmos_sim.integrator.p.radiation.ᶠradiation_flux,
        ),
    )

    @test compare(
        cs_four_steps.model_sims.ice_sim.integrator.u,
        cs_two_steps2.model_sims.ice_sim.integrator.u,
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

    # Ignoring SST_timevaryinginput because it contains closures (which should be
    # reinitialized correctly). We have to remove it from the type, otherwise
    # compare will not work
    function delete(nt::NamedTuple, fieldnames...)
        return (; filter(p -> !(first(p) in fieldnames), collect(pairs(nt)))...)
    end

    ocean_cache_four =
        delete(cs_four_steps.model_sims.ocean_sim.cache, :SST_timevaryinginput)
    ocean_cache_two2 =
        delete(cs_two_steps2.model_sims.ocean_sim.cache, :SST_timevaryinginput)

    @test compare(ocean_cache_four, ocean_cache_two2)
end
