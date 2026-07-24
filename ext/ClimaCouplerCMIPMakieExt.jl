"""
    ClimaCouplerCMIPMakieExt

This module contains code for plotting the diagnostics of the Oceananigans-based
ocean and sea ice components of a coupled simulation.

The ocean and sea ice write their surface diagnostics as Oceananigans JLD2
`FieldTimeSeries` output (`ocean_surface`, `seaice_surface`), on their native
(possibly curvilinear, e.g. tripolar) grids. This extension reads those series
and plots each field on a global map using the shared plot styling (Robinson
projection, coastlines, and per-variable colormaps from the `Plotting` module),
so that ocean/sea ice diagnostics match the atmosphere/land diagnostics and the
instantaneous snapshots.
"""
module ClimaCouplerCMIPMakieExt

using ClimaCoupler
import ClimaCoupler: Plotting

using CairoMakie
using GeoMakie
using Makie
using Printf
import Oceananigans as OC

include(joinpath("ClimaCouplerCMIPMakieExt", "ocean_seaice_diagnostics_plots.jl"))

end
