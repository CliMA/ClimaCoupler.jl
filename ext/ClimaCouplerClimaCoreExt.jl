"""
    ClimaCouplerClimaCoreExt

Extension providing spectrum computation for ClimaAnalysis OutputVars using ClimaCoreSpectra.

When both ClimaCore and ClimaCoreSpectra are loaded, this extension provides `compute_spectrum`,
which computes the spherical power spectrum of a 3D (lon, lat, z) or 4D (time, lon, lat, z)
ClimaAnalysis.OutputVar, as used in ClimaAtmos ci_plots and diagnostics.
"""
module ClimaCouplerClimaCoreExt

using ClimaCoupler
import ClimaAnalysis
import ClimaCore
import ClimaCoreSpectra: power_spectrum_2d

include(joinpath("ClimaCouplerClimaCoreExt", "spectra_utilities.jl"))

end
