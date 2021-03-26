module CouplerMachine

using ClimateMachine

include(joinpath("Coupling", "Coupling.jl"))
include(joinpath("Coupling", "CplSolver.jl"))
end
