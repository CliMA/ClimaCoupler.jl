using SafeTestsets

@safetestset "Aqua tests" begin
    include("aqua.jl")
end
@safetestset "Clock tests" begin
    include("CoupledSimulations/clock.jl")
end
@safetestset "Interfacer tests" begin
    include("interfacer_tests.jl")
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
@safetestset "FieldExchanger tests" begin
    include("field_exchanger_tests.jl")
end
@safetestset "FluxCalculator tests" begin
    include("flux_calculator_tests.jl")
end
@safetestset "Diagnostics tests" begin
    include("diagnostics_tests.jl")
end
@safetestset "PostProcessor tests" begin
    include("postprocessor_tests.jl")
end
@safetestset "Checkpointer tests" begin
    include("checkpointer_tests.jl")
end
@safetestset "CouplerState tests" begin
    include("CouplerState/cplstate_interface.jl")
end
# include("CoupledSimulations/cplsolver.jl")
