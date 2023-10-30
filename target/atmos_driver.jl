using ClimaComms
using Logging
using ClimaAtmos

const is_distributed = get(ENV, "CLIMACORE_DISTRIBUTED", "") == "MPI"
if is_distributed
    const comms_ctx = ClimaComms.context(ClimaComms.CPUSingleThreaded())
    const pid, nprocs = ClimaComms.init(comms_ctx)

    if ClimaComms.iamroot(comms_ctx)
        Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Info))
    else
        Logging.global_logger(Logging.NullLogger())
    end
else
    using TerminalLoggers: TerminalLogger
    const comms_ctx = ClimaComms.SingletonCommsContext()
    const pid, nprocs = ClimaComms.init(comms_ctx)
end


ClimaComms.barrier(comms_ctx)

# this tests the standalone runs, using the hybrid/driver
import Random
Random.seed!(1234)

if !(@isdefined config)
    config = ClimaAtmos.AtmosConfig(comms_ctx = comms_ctx)
end

ClimaComms.iamroot(comms_ctx) ? @info(config) : nothing
integrator = ClimaAtmos.get_integrator(config)
ClimaComms.barrier(comms_ctx)
sol_res = ClimaAtmos.solve_atmos!(integrator)
