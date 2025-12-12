"""
    ClimaCouplerPlotsExt

Extension module for ClimaCoupler that provides plotting functionality when Makie is available.
This extension loads the plotting submodules into the Postprocessor module.
"""
module ClimaCouplerPlotsExt

using ..ClimaCoupler
using Makie
using CairoMakie
using ClimaCoreMakie

include("ClimaCouplerPlotsExt/diagnostics_plots.jl")
include("ClimaCouplerPlotsExt/debug_plots.jl")

end # module ClimaCouplerPlotsExt
