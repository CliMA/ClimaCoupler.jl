using SafeTestsets
import ClimaComms
ClimaComms.@import_required_backends

# load extensions for AMIP and plotting
include(joinpath("..", "AMIP", "code_loading.jl"))

@safetestset "component model test: ClimaAtmos" begin
    include("component_model_tests/climaatmos_tests.jl")
end
@safetestset "component model test: ClimaLand integrated model" begin
    include("component_model_tests/climaland_tests.jl")
end
@safetestset "surface radiative flux consistency tests" begin
    include("fluxes_test.jl")
end
@safetestset "debug diagnostics: debug plots" begin
    include("debug_plots_tests.jl")
end
@safetestset "AMIP test" begin
    include("amip_test.jl")
end
@safetestset "time incrementing test" begin
    include("time_increment_test.jl")
end
