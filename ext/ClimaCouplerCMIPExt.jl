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
import ConservativeRegridding as CR
import Adapt # for ConservativeRegridding

"""
    OceananigansSimulation{SIM, A, OPROP, REMAP, SIC}

The ClimaCoupler simulation object used to run with Oceananigans.
This type is used by the coupler to indicate that this simulation
is a surface/ocean simulation for dispatch.

It contains the following objects:
- `ocean::SIM`: The Oceananigans simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
- `ocean_properties::OPROP`: A NamedTuple of ocean properties and parameters (including COARE3 roughness params).
- `remapping::REMAP`: Objects needed to remap from the exchange (spectral) grid to Oceananigans spaces.
- `ice_concentration::SIC`: An Oceananigans Field representing the sea ice concentration on the ocean/sea ice grid.
"""
struct OceananigansSimulation{SIM, A, OPROP, REMAP, SIC} <:
       Interfacer.AbstractOceanSimulation
    ocean::SIM
    area_fraction::A
    ocean_properties::OPROP
    remapping::REMAP
    ice_concentration::SIC
end

"""
    ClimaSeaIceSimulation{SIM, A, REMAP, NT, IP}

The ClimaCoupler simulation object used to run with ClimaSeaIce.
This type is used by the coupler to indicate that this simulation
is a surface/sea ice simulation for dispatch.

It contains the following objects:
- `ice::SIM`: The ClimaSeaIce simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
- `remapping::REMAP`: Objects needed to remap from the exchange (spectral) grid to Oceananigans spaces.
- `ocean_ice_interface::NT`: A NamedTuple containing fluxes between the ocean and sea ice, computed at each coupling step,
                             the interfacial temperature and salinity, and the flux formulation used to compute the fluxes.
- `ice_properties::IP`: A NamedTuple of sea ice properties, including melting speed, Stefan-Boltzmann constant,
    and the Celsius to Kelvin conversion constant.
"""
struct ClimaSeaIceSimulation{SIM, A, REMAP, NT, IP} <: Interfacer.AbstractSeaIceSimulation
    ice::SIM
    area_fraction::A
    remapping::REMAP
    ocean_ice_interface::NT
    ice_properties::IP
end

# Include helper functions first (used by both oceananigans.jl and clima_seaice.jl)
include("ClimaCouplerCMIPExt/climaocean_helpers.jl")

# Include the model files
include("ClimaCouplerCMIPExt/oceananigans.jl")
include("ClimaCouplerCMIPExt/clima_seaice.jl")

end # module ClimaCouplerCMIPExt
