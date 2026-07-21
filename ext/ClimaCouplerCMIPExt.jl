"""
    ClimaCouplerCMIPExt

Module containing CMIP component models (Oceananigans and ClimaSeaIce models)
which extend the ClimaCoupler.jl simulation interface.

This extension is loaded when Oceananigans, ClimaOcean, ClimaSeaIce,
ClimaAtmos, KernelAbstractions, ConservativeRegridding, and Adapt are loaded
with either `import` or `using`.


For more information about the CMIP component models, please see the
"Available component models" section of the ClimaCoupler documentation,
or see the [Oceananigans documentation](https://clima.github.io/OceananigansDocumentation/stable/)
and [ClimaSeaIce documentation](https://clima.github.io/ClimaSeaIceDocumentation/dev/).
"""
module ClimaCouplerCMIPExt

import ClimaCoupler
import ClimaCoupler:
    Checkpointer,
    FieldExchanger,
    FluxCalculator,
    Interfacer,
    TimeManager,
    Utilities,
    Plotting
import Oceananigans as OC
import ClimaOcean as CO
import ClimaSeaIce as CSI
import ClimaAtmos as CA # for basis conversions (projected_vector_data)
import ClimaCore as CC
import ClimaParams as CP
using KernelAbstractions: @kernel, @index, @inbounds

import Adapt # for ConservativeRegridding
import ClimaCore as CC # for ConservativeRegriddingClimaCoreExt
import ConservativeRegridding as CR
import SparseArrays # for converting Regridder element types

get_ConservativeRegriddingCCExt() =
    Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)

get_ConservativeRegriddingOCExt() =
    Base.get_extension(CR, :ConservativeRegriddingOceananigansExt)

# Exchange-grid geometry/weights and per-polygon flux machinery, used by the
# ocean and sea-ice models below
include("ClimaCouplerCMIPExt/exchange_grid.jl")
include("ClimaCouplerCMIPExt/exchange_fluxes.jl")

# Include the model files first so their types are available to climaocean_helpers.jl
include("ClimaCouplerCMIPExt/oceananigans.jl")
include("ClimaCouplerCMIPExt/clima_seaice.jl")
include("ClimaCouplerCMIPExt/climaocean_helpers.jl")

# Include skin temperature utilities
include("ClimaCouplerCMIPExt/skin_temperature.jl")

include("ClimaCouplerCMIPExt/ocean_diagnostics.jl")
include("ClimaCouplerCMIPExt/seaice_diagnostics.jl")

end # module ClimaCouplerCMIPExt
