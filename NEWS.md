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

### A higher resolution land-sea mask is now used and automatically downloaded - PR [#1006](https://github.com/CliMA/ClimaCoupler.jl/pull/1006)

A 60 arcsecond land-sea mask constructed from topographic data is now used.
Topographic data is automatically downloaded and a land-sea mask is constructed
by identifying where elevation is greater than 0. Note, this can lead to
misidentification of ocean in some areas of the globe that are inland but below
sea level (Dead Sea, Death Valley, ...).



### Code cleanup
#### Output path update - PR [#1058](https://github.com/CliMA/ClimaCoupler.jl/pull/1058)
Previously, ClimaEarth simulation outputs were saved in a path
`experiments/ClimaEarth/output/$mode_name/$job_id/artifacts/`.
This PR removes `mode_name` has from this pattern, so output will now be in
`experiments/ClimaEarth/output/$job_id/artifacts/`.
Note that any external scripts that assume an output path will need to be updated.

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

#### Remove PostProcessor module - PR [#1022](https://github.com/CliMA/ClimaCoupler.jl/pull/1022)
After switching to use ClimaDiagnostics.jl and ClimaAnalysis.jl
for our diagnostics, the PostProcessor module is not needed.
This PR follows directly from the Diagnostics module removal.

### Maintenance
- Update to JuliaFormatter v2. PR [#1024](https://github.com/CliMA/ClimaCoupler.jl/pull/1024)
- Update CI to use Julia v1.11. Introduce Manifest files for Julia 1.11, in addition to the existing generic Manifests. PR [#1026](https://github.com/CliMA/ClimaCoupler.jl/pull/1026)

#### Various ClimaEarth cleanup - PR [#1070](https://github.com/CliMA/ClimaCoupler.jl/pull/1070)
This PR does a few cleanup tasks in the experiments/ClimaEarth/ directory:
- Update the ClimaEarth README.
- Delete the `viz_explorer.jl`, which was previously used to create animations of simulation fields, but is not currently being used. Note that these animations will no longer appear in buildkite output.
- Move functions in the `io_helpers.jl` file to the Utilities module and delete this file.
This should not change any behavior.
- Move `checkpoint_sims` function to the Checkpointer module and delete the `user_logging.jl` file.
This should not change any behavior.
