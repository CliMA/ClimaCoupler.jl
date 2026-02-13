"""
    ClimaCouplerCMIPExt

Module containing CMIP component models (Oceananigans and ClimaSeaIce models)
which extend the ClimaCoupler.jl simulation interface.

This extension is loaded when Oceananigans, ClimaOcean, ClimaSeaIce,
and KernelAbstractions are loaded with either `import` or `using`.
"""
module ClimaCouplerCMIPExt

import ClimaCoupler
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import Oceananigans as OC
import ClimaOcean as CO
import ClimaSeaIce as CSI
import ClimaCore as CC
import ClimaParams as CP
using KernelAbstractions: @kernel, @index, @inbounds

# Include helper functions first (used by both oceananigans.jl and clima_seaice.jl)
include("ClimaCouplerCMIPExt/climaocean_helpers.jl")

# Include skin temperature utilities
include("ClimaCouplerCMIPExt/skin_temperature.jl")

# Include the model files
include("ClimaCouplerCMIPExt/oceananigans.jl")
include("ClimaCouplerCMIPExt/clima_seaice.jl")

end # module ClimaCouplerCMIPExt
