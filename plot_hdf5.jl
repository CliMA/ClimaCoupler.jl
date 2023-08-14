using Plots
using ClimaCorePlots
using ClimaCore.InputOutput # https://clima.github.io/ClimaCore.jl/dev/api/#InputOutput
datafile = ClimaCore.InputOutput.HDF5Reader("<filename>")
variable = ClimaCore.InputOutput.read_field(datafile, "Y")
Y |> propertynames # Explore variables available in HDF5 Dataset
Plots.plot(variable)  # Plots on an unwrapped cubed-sphere 

## TODO: reconfigure plots from amip pipeline run