"""
    ClimaCouplerCMIPExt

Module containing CMIP component models (Oceananigans and ClimaSeaIce models)
which extend the ClimaCoupler.jl simulation interface.

This extension is loaded when Oceananigans, ClimaOcean, ClimaSeaIce,
and KernelAbstractions are loaded with either `import` or `using`.

For more information about the CMIP component models, please see the
"Available component models" section of the ClimaCoupler documentation,
or see the [Oceananigans documentation](https://clima.github.io/OceananigansDocumentation/stable/)
and [ClimaSeaIce documentation](https://clima.github.io/ClimaSeaIceDocumentation/dev/).
"""
module ClimaCouplerCMIPExt

import ClimaCoupler
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities, Plotting
import Oceananigans as OC
import ClimaOcean as CO
import ClimaSeaIce as CSI
import ClimaCore as CC
import ClimaComms
import ClimaParams as CP
import ConservativeRegridding as CR
import Adapt
import LinearAlgebra
import Logging
using KernelAbstractions: @kernel, @index, @inbounds

# Include helper functions first (used by both oceananigans.jl and clima_seaice.jl)
include("ClimaCouplerCMIPExt/climaocean_helpers.jl")

# Include skin temperature utilities
include("ClimaCouplerCMIPExt/skin_temperature.jl")

# Include the model files
include("ClimaCouplerCMIPExt/oceananigans.jl")
include("ClimaCouplerCMIPExt/clima_seaice.jl")

# Conservative regridding + `Interfacer.remap!`/`remap` extensions (after simulation types exist)
include("ClimaCouplerCMIPExt/cmip_coupler_regrid.jl")

end # module ClimaCouplerCMIPExt
