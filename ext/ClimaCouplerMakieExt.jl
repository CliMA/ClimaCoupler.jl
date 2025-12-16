"""
    ClimaCouplerMakieExt

This module contains code for plotting simulation output.

Currently, it includes:
- diagnostics plots
- plots of model and coupler fields for debugging
- leaderboard comparing to observations
- parameter calibration plots
"""
module ClimaCouplerMakieExt

using ClimaCoupler
import ClimaCoupler: Plotting

using CairoMakie
using ClimaCoreMakie
using GeoMakie
using Makie

# Diagnostics plots
include(joinpath("ClimaCouplerMakieExt", "diagnostics_plots.jl"))

# Debug plots
include(joinpath("ClimaCouplerMakieExt", "debug_plots.jl"))

# Leaderboard
include(joinpath("ClimaCouplerMakieExt", "leaderboard", "data_sources.jl"))
include(joinpath("ClimaCouplerMakieExt", "leaderboard", "leaderboard.jl"))

end
