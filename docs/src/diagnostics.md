# Diagnostics

ClimaCoupler.jl utilizes other packages in the CliMA ecosystem to generate
and visualize diagnostics, namely ClimaDiagnostics.jl and ClimaAnalysis.jl.

## Default AMIP diagnostics
We output a set of default diagnostics for all AMIP runs.
These currently include the following fields:

Atmospheric quantiies:
- air temperature at the bottom of the atmosphere (3D)
- eastward near-surface wind (3D)
- specific humidity (3D)
- mass fraction of cloud liquid water (3D)
- net top-of-atmosphere fluxes (3D)
- precipitation (2D)
- surface temperature (2D)

Coupler quantities
- turbulent energy fluxes (2D)

These diagnostics are all averaged over a period of time that depends
on the length of the overall simulation according to the following rule:
- simulation length >= 90 days: 30-day mean
- simulation length >= 30 days and < 90 days: 10-day mean
- simulation length >= 1 day and < 30 days: 1-day mean
- simulation length < 1 day: 1-hour mean

## How to add a new diagnostic variable
### Adding a diagnostic for a ClimaCoupler quantity
For diagnostics that come from coupler fields or that are computed using input
from multiple component models, we set up the diagnostics by directly creating
ClimaDiagnostics.jl objects.

Specifically, we first create a `DiagnosticVariable` object containing the variable's name,
units, any comments, and the function to compute it. This is then used to create a
`ScheduledDiagnostic` variable, which includes saving and output time information.
Once we have created a `ScheduledDiagnostic` for each variable we're interested in,
we collect them in a vector and pass this to our `DiagnosticsHandler` object.

An example of this process for the combined turbulent energy flux, `F_turb_energy`, can be found in
`experiments/ClimaEarth/user_io/coupler_diagnostics.jl`.

For more information about this process, please see the
ClimaDiagnostics.jl [documentation](https://clima.github.io/ClimaDiagnostics.jl/dev/).

### Adding a diagnostic for a CliMA component model quantity
Many of our current diagnostics are values that we access directly from a component model.
To add a new diagnostic of this kind, you can add a new method `add_diagnostic_variable!`
extending this function from the component model package's Diagnostics module.

For more information about this function and the form it takes, please see the
ClimaDiagnostics.jl [documentation](https://clima.github.io/ClimaDiagnostics.jl/dev/).

#### ClimaAtmos.jl
To add a new diagnostic for the ClimaAtmos.jl atmosphere model, you can add this new method
for `ClimaAtmos.Diagnostics.add_diagnostic_variable!` in
`components/atmosphere/climaatmos_extra_diags.jl`.
The existing diagnostics in that file can be used as templates.

For more information about ClimaAtmos diagnostics, and to see the default atmospheric diagnostics,
please see that package's [documentation](https://clima.github.io/ClimaAtmos.jl/dev/diagnostics/).

#### ClimaLand.jl
To add a new diagnostic for the ClimaLand.jl bucket model, you can add this new method
for `ClimaLand.Diagnostics.add_diagnostic_variable!` in
`components/land/climaland_bucket_extra_diags.jl` (which doesn't exist at the time of writing).

For more information about ClimaLand diagnostics, and to see the default land diagnostics,
please see that package's [documentation](https://clima.github.io/ClimaLand.jl/dev/diagnostics/users_diagnostics/).

## Visualizing diagnostics
ClimaCoupler.jl uses ClimaAnalysis.jl to parse and visualize the outputs saved
using ClimaDiagnostics.jl.

For more information about ClimaAnalysis.jl, please see that package's [documentation](https://clima.github.io/ClimaAnalysis.jl/dev/).
