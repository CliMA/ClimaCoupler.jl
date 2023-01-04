using SafeTestsets

@safetestset "Clock tests" begin
    include("CoupledSimulations/clock.jl")
end
@safetestset "Regridder tests" begin
    include("regridder_tests.jl")
end
@safetestset "ConservationChecker tests" begin
    include("conservation_checker_tests.jl")
end
@safetestset "BCReader tests" begin
    include("bcreader_tests.jl")
end
@safetestset "Utilities tests" begin
    include("utilities_tests.jl")
end
@safetestset "TimeManager tests" begin
    include("time_manager_tests.jl")
end
@safetestset "Diagnostics tests" begin
    include("diagnostics_tests.jl")
end
@safetestset "CouplerState tests" begin
    include("CouplerState/cplstate_interface.jl")
end
# include("CoupledSimulations/cplsolver.jl")
