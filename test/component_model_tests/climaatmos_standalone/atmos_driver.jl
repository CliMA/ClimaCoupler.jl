import ClimaAtmos as CA
import Random

redirect_stderr(IOContext(stderr, :stacktrace_types_limited => Ref(false)))
Random.seed!(1234)

if !(@isdefined config)
    config = CA.AtmosConfig()
end
simulation = CA.get_simulation(config)
(; integrator) = simulation
sol_res = CA.solve_atmos!(simulation)
