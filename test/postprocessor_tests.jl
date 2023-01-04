#= 
    Unit tests for ClimaCoupler Diagnostics module
=#
using Test
using Dates
using ClimaComms
using ClimaCoupler: Utilities, TestHelper

FT = Float64
get_var(cs::CoupledSimulation, ::Val{:x}) = FT(1)

@testset "postprocess" begin

    
end
