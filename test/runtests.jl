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
@safetestset "TimeManager tests" begin
    include("time_manager_tests.jl")
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
@safetestset "component test: bucket" begin
    include("component_model_tests/bucket_tests.jl")
end
@safetestset "component model test: ClimaAtmos" begin
    include("component_model_tests/climaatmos_tests.jl")
end
@safetestset "component model test: prescr. sea ice" begin
    include("component_model_tests/prescr_seaice_tests.jl")
end
@safetestset "component model test: eisenman sea ice" begin
    include("component_model_tests/eisenman_seaice_tests.jl")
end
@safetestset "component model test: slab ocean" begin
    include("component_model_tests/slab_ocean_tests.jl")
end


# include("CoupledSimulations/cplsolver.jl")
