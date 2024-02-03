using ClimaComms
using Logging
using ClimaAtmos

redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))
import ClimaAtmos as CA
import Random
Random.seed!(1234)

if !(@isdefined config)
    config = CA.AtmosConfig()
end
simulation = CA.get_simulation(config)
(; integrator) = simulation
sol_res = CA.solve_atmos!(simulation)
