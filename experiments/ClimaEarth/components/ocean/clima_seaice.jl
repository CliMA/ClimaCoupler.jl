import Oceananigans as OC
import ClimaOcean as CO
import ClimaCoupler: Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities
import ClimaComms
import ClimaCore as CC
import Thermodynamics as TD
import ClimaOcean.EN4: download_dataset
using KernelAbstractions: @kernel, @index, @inbounds

"""
    ClimaSeaIceSimulation{SIM, A, OPROP, REMAP}

The ClimaCoupler simulation object used to run with Oceananigans.
This type is used by the coupler to indicate that this simulation
is an surface/ocean simulation for dispatch.

It contains the following objects:
- `sea_ice::SIM`: The ClimaSeaIce simulation object.
- `area_fraction::A`: A ClimaCore Field representing the surface area fraction of this component model on the exchange grid.
- `melting_speed::MS`: An constant characteristic speed for melting/freezing.
- `remapping::REMAP`: Objects needed to remap from the exchange (spectral) grid to Oceananigans spaces.
"""
struct ClimaSeaIceSimulation{SIM,A,MS,REMAP} <: Interfacer.SeaIceModelSimulation
    sea_ice::SIM
    area_fraction::A
    melting_speed::MS
    remapping::REMAP
end

"""
    ClimaSeaIceSimulation()

Creates an OceananigansSimulation object containing a model, an integrator, and
a surface area fraction field.
This type is used to indicate that this simulation is an ocean simulation for
dispatch in coupling.

Specific details about the default model configuration
can be found in the documentation for `ClimaOcean.ocean_simulation`.
"""
function ClimaSeaIceSimulation(area_fraction, ocean; output_dir)
    # Initialize the sea ice with the same grid as the ocean
    # Initially no sea ice is present
    grid = ocean.ocean.model.grid
    advection = ocean.ocean.model.advection.T
    sea_ice = CO.sea_ice_simulation(grid, ocean.ocean; advection)

    melting_speed = 1e-4
    remapping = ocean.remapping

    # Before version 0.96.22, the NetCDFWriter was broken on GPU
    if arch isa OC.CPU || pkgversion(OC) >= v"0.96.22"
        # TODO maybe this is broken
        # Save all tracers and velocities to a NetCDF file at daily frequency
        outputs = prognostic_fields(sea_ice.model)
        netcdf_writer = OC.NetCDFWriter(
            sea_ice.model,
            outputs;
            schedule = OC.TimeInterval(86400), # Daily output
            filename = joinpath(output_dir, "seaice_diagnostics.nc"),
            indices = (:, :, grid.Nz),
            overwrite_existing = true,
            array_type = Array{Float32},
        )
        sea_ice.output_writers[:diagnostics] = netcdf_writer
    end

    sim = ClimaSeaIceSimulation(sea_ice, area_fraction, melting_speed, remapping)
    return sim
end
