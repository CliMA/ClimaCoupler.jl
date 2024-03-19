using ClimaComms
using Logging

redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))
import ClimaAtmos as CA
import Random
using ClimaCoupler
Random.seed!(1234)

pkg_dir = pkgdir(ClimaCoupler)
atmos_config_file = joinpath(pkg_dir, "test/component_model_tests/climaatmos_standalone/longrun_aquaplanet_dyamond.yml")
comms_ctx = ClimaComms.context(ClimaComms.device())
@show typeof(comms_ctx)

config = CA.AtmosConfig(atmos_config_file; comms_ctx = comms_ctx)

OUTPUT_DIR = joinpath(pkg_dir, "test/component_model_tests/climaatmos_standalone/output/longrun_aquaplanet_dyamond_artifacts/")
mkpath(OUTPUT_DIR)

if !(@isdefined config)
    config = CA.AtmosConfig()
end
simulation = CA.get_simulation(config)
(; integrator) = simulation
sol_res = CA.solve_atmos!(simulation)
