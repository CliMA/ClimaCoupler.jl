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

# Instantaneous snapshot plots
include(joinpath("ClimaCouplerMakieExt", "snapshot_plots.jl"))

# Conservation plots
include(joinpath("ClimaCouplerMakieExt", "conservation_plots.jl"))

# Leaderboard
include(joinpath("ClimaCouplerMakieExt", "leaderboard", "data_sources.jl"))
include(joinpath("ClimaCouplerMakieExt", "leaderboard", "leaderboard.jl"))

end
