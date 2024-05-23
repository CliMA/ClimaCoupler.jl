# Hierarchy Experiments for Global Climate Models

This directory contains a series of experiments that demonstrate the use of hierarchical modeling using ClimaEarth, with a focus on the atmospheric component. The hierarchy spans from a simple dry atmosphere to a more complex moist atmosphere with clouds and an Earth-like surface. The experiments are designed to be run in sequence, with each experiment building on the previous one. Each experiment is a self-contained run script that contains all necessary configurations.

## Experiments
- Dry Held-Suarez
- Moist Held-Suarez
- Cloudless Aquaplanet
- Cloudy Aquaplanet
- Cloudy Slabplanet

## Postprocessing
We provide a simple postprocessing script `climate_plots.jl` that can be used to visualize the results of each experiment, with some helper functions in `plot_helper.jl`.

## Associated publication
- (in preparation)
