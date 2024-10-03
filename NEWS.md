ClimaCoupler.jl Release Notes
============================

`main`
-------

### ClimaEarth features

### Sea-surface temperature and sea ice concentration data can now be automatically downloaded

Sea-surface temperature and sea ice concentration require external files. Now, a
low-resolution version of such files is automatically downloaded when a
higher-resolution version is not available. Please, refer to
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts) for more information.

### Code cleanup
#### Remove ClimaCoupler.Diagnostics module - PR [#953](https://github.com/CliMA/ClimaCoupler.jl/pull/953)
The ClimaCoupler Diagnostics module had become redundant with
ClimaDiagnostics.jl, a package designed to provide robust
diagnostics across the CliMA ecosystem.
Here we remove ClimaCoupler.Diagnostics and instead use
ClimaDiagnostics.

We're able to retrieve most of the diagnostics
we want directly from ClimaAtmos and ClimaLand, but also want
some that come from coupler-computed quantities, such as
`F_turb_energy`. In this PR we add this coupler quantity
to our output diagnostics using the ClimaDiagnostics interface.

This PR also removes the AMIP paperplots function, but this
functionality is replaced by the generalized `make_plots` function.
