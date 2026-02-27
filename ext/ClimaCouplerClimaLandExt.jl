"""
    ClimaCouplerClimaLandExt

This module contains code for extending the ClimaCoupler interface for
ClimaLand. This includes two land model options: a simple bucket or the more
complex integrated land model. For more information about each model, please see
`experiments/ClimaEarth/README.md`
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
import ClimaUtilities.ClimaArtifacts: @clima_artifact
import ClimaCoupler:
    Checkpointer, FieldExchanger, FluxCalculator, Interfacer, Utilities, Plotting
import SciMLBase
import ClimaTimeSteppers as CTS
import ClimaDiagnostics as CD
import ClimaUtilities.TimeManager: ITime
import SurfaceFluxes as SF
import SurfaceFluxes.Parameters as SFP
import Thermodynamics as TD
import ClimaComms
import ClimaUtilities.TimeManager: ITime
using NCDatasets
import StaticArrays
import Interpolations

include("ClimaCouplerClimaLandExt/climaland_helpers.jl")
include("ClimaCouplerClimaLandExt/climaland_bucket.jl")
include("ClimaCouplerClimaLandExt/climaland_integrated.jl")

end # module
