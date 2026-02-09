"""
    ClimaCouplerOceananigansMakieExt

This module contains code for plotting output of coupled simulations
including Oceananigans.

Currently, it includes:
- diagnostics plots for Oceananigans
"""
module ClimaCouplerOceananigansMakieExt

using ClimaCoupler
import ClimaCoupler: Plotting

using CairoMakie
using ClimaCoreMakie
using GeoMakie
using Makie
using Printf
import Oceananigans as OC

"""
    Plotting.debug_plot!(ax, fig, field, i, j)

Plot a heatmap of the provided Oceananigans field or operation.
This is intended to be used as part of the debug plotting system.
"""
function Plotting.debug_plot!(ax, fig, field::OC.Field, i, j)
    grid = field.grid
    hm = CairoMakie.heatmap!(ax, view(field, :, :, grid.Nz))
    Makie.Colorbar(fig[i, j * 2], hm)
    return nothing
end
function Plotting.debug_plot!(ax, fig, field::OC.AbstractOperations.AbstractOperation, i, j)
    # Evaluate the operation to get a field
    evaluated_field = OC.Field(field)
    OC.compute!(evaluated_field)

    # Plot the evaluated field
    Plotting.debug_plot!(ax, fig, evaluated_field, i, j)
    return nothing
end

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
