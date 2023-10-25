using Pkg
Pkg.activate("../experiments/AMIP/modular/")

import Pkg; Pkg.add("ClimaComms")
import Pkg; Pkg.add("Logging")
import Pkg; Pkg.add("ClimaAtmos")
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


Pkg.activate(joinpath(pkgdir(ClimaAtmos), "examples"))
ClimaComms.barrier(comms_ctx)

include(joinpath(pkgdir(ClimaAtmos), "examples/hybrid/driver.jl"))
