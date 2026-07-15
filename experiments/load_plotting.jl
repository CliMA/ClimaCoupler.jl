#=
Trigger ClimaCouplerMakieExt (and ClimaCouplerCMIPMakieExt when Oceananigans is loaded).

Include this file whenever plotting is needed — either before `run!` (mid-run /
callback plots) or after `run!` (postprocess-only). Keeping this out of
`code_loading.jl` avoids paying Makie compilation cost on every simulation.

Usage:
    include(joinpath(@__DIR__, "load_plotting.jl"))            # from experiments/
    include(joinpath(@__DIR__, "..", "load_plotting.jl"))       # from AMIP/ or CMIP/
=#

using CairoMakie, ClimaCoreMakie, GeoMakie, Makie, Poppler_jll, Printf
