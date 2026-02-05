"""
    ClimaCouplerClimaLandExt

This module contains code for extending the ClimaCoupler interface for
ClimaLand.
"""
module ClimaCouplerClimaLandExt

import ClimaCore as CC
import ClimaLand as CL
import ClimaParams as CP

import ClimaLand.Parameters as LP
import Dates
import ClimaUtilities.TimeVaryingInputs:
    LinearInterpolation, PeriodicCalendar, TimeVaryingInput
import ClimaUtilities.SpaceVaryingInputs: SpaceVaryingInput
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities, Plotting
import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaDiagnostics as CD
import ClimaUtilities.TimeManager: ITime
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD

import Statistics
import ClimaComms
import ClimaUtilities.TimeManager: ITime
using NCDatasets
import StaticArrays

include("climaland/climaland_helpers.jl")
include("climaland/climaland_bucket.jl")
include("climaland/climaland_integrated.jl")

end # module
