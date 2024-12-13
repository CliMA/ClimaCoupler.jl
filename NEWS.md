ClimaCoupler.jl Release Notes
============================

`main`
-------

### ClimaEarth features

### Read bucket initial conditions from NetCDF files

Added functionality to allow the bucket initial conditions to be overwritten by interpolated NetCDF datasets.
To use this feature from the YAML interface, just pass the path of the file to `land_initial_condition`. 
We expect the file to contain the following variables:
`W`, for subsurface water storage (2D),
`Ws`, for surface water content (2D),
`T`, for soil temperature (3D),
`S`, for snow water equivalent (2D).

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

### Leaderboard for variables over longitude, latitude, time, and pressure - PR [#1094](https://github.com/CliMA/ClimaCoupler.jl/pull/1094)

As a part of the post processing pipeline, bias plots for variables at the
pressure levels of 850.0, 500.0, 250.0 hPa and bias plots over latitude and
pressure levels are being created.

### Code cleanup

#### Output path updates - PRs [#1058](https://github.com/CliMA/ClimaCoupler.jl/pull/1058),
    [#1106](https://github.com/CliMA/ClimaCoupler.jl/pull/1106),
    [#1123](https://github.com/CliMA/ClimaCoupler.jl/pull/1123)


Previously, ClimaEarth simulation outputs were saved in a path
`experiments/ClimaEarth/output/$mode_name/$job_id/artifacts/`. Now, `ClimaEarth`
creates output folders with an increment (increasing the counter every time the
simulation is run). This is in preparation to restarts. The output now looks
like
```
coupler_output_dir_amip/
├── checkpoints
│       └── checkpoints for the various models
├── artifacts
│       └── plots produced by the postprocessing step
├── output_0000/
│   ├── clima_atmos/
│   │   └── output of the atmos model
│   └── clima_land/
│   │   └── output of the land model
│   └── coupler/
│   │   └── output of coupler quantities
│   └── ocean/
│       └── output of the ocean model
├── output_0001/
│   └── ... component model outputs in their folders ...
├── output_0002/
│   └── ... component model outputs in their folders ...
└── output_active -> output_0002/
``
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
