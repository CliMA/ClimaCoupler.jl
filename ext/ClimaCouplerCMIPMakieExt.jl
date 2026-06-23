"""
    ClimaCouplerCMIPMakieExt

This module contains code for plotting output of coupled simulations
including Oceananigans.

Currently, it includes:
- diagnostics plots for Oceananigans
- movies of Oceananigans JLD2 diagnostics
- interactive globe viewing of ocean diagnostics
"""
module ClimaCouplerCMIPMakieExt

using ClimaCoupler
import ClimaCoupler: Plotting

using CairoMakie
using ClimaCoreMakie
using GeoMakie
using JLD2
using Makie
using Printf
import Oceananigans as OC

include(joinpath(@__DIR__, "ClimaCouplerCMIPMakieExt", "ocean_diagnostics_movies.jl"))

"""
    print_extrema(field)

Return a tuple `(string, has_nan::Bool)` where:
- `string` is the minimum and maximum values of an Oceananigans field formatted as a string
- `has_nan` is true if the extrema contain NaN values
"""
function Plotting.print_extrema(field::OC.Field)
    min_val = minimum(field)
    max_val = maximum(field)

    # Check for NaN values
    has_nan = isnan(min_val) || isnan(max_val)

    min = Printf.@sprintf("%.2E", min_val)
    max = Printf.@sprintf("%.2E", max_val)
    return (" [$min, $max]", has_nan)
end
function Plotting.print_extrema(operation::OC.AbstractOperations.AbstractOperation)
    evaluated_field = OC.Field(operation)
    OC.compute!(evaluated_field)
    return Plotting.print_extrema(evaluated_field)
end

end
