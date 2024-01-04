# TODO: const comms_ctx cannot be overwritten in ClimaAtmos - add if statement to clima atmos
using ClimaComms
using Logging

const is_distributed = get(ENV, "CLIMACORE_DISTRIBUTED", "") == "MPI"
if is_distributed
    const comms_ctx = ClimaComms.context(ClimaComms.CPUSingleThreaded())
    const pid, nprocs = ClimaComms.init(comms_ctx)

    #=
    if ClimaComms.iamroot(comms_ctx)
        Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Info))
    else
        Logging.global_logger(Logging.NullLogger())
    end
    =#
else
    using TerminalLoggers: TerminalLogger
    const comms_ctx = ClimaComms.SingletonCommsContext()
    const pid, nprocs = ClimaComms.init(comms_ctx)
end
