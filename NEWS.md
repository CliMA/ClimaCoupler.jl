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


#### Remove ClimaCoupler.Regridder module - PR [#1109](https://github.com/CliMA/ClimaCoupler.jl/pull/1109)
This PR removes the Regridder module. Most of the functions that were
inside the module are available in [ClimaUtilities](https://github.com/CliMA/ClimaUtilities.jl/). Some functions were moved
to other modules, and some unused functions were deleted.

The functions:
- `ClimaCoupler.Regridder.dummmy_remap!`
- `ClimaCoupler.Regridder.update_surface_fractions!`
- `ClimaCoupler.Regridder.combine_surfaces!`
- `ClimaCoupler.Regridder.binary_mask`

are now:
- `ClimaCoupler.FieldExchanger.dummmy_remap!`
- `ClimaCoupler.FieldExchanger.update_surface_fractions!`
- `ClimaCoupler.FieldExchanger.combine_surfaces!`
- `ClimaCoupler.Utilities.binary_mask`

The following functions were removed from ClimaCoupler and have an equivalent in [ClimaUtilities](https://github.com/CliMA/ClimaUtilities.jl/tree/main):

- `ClimaCoupler.Regridder.write_to_hdf5`
- `ClimaCoupler.Regridder.read_from_hdf5`
- `ClimaCoupler.Regridder.hdwrite_regridfile_rll_to_cgll`
- `ClimaCoupler.Regridder.reshape_cgll_sparse_to_field!`
- `ClimaCoupler.Regridder.get_time`

Note that the the `hdwrite_regridfile_rll_to_cgll` in [ClimaUtilities](https://github.com/CliMA/ClimaUtilities.jl/tree/main) does not support 3d fields, but the removed function from the Regridder module did.

The following functions were deleted and do not have an equivalent in [ClimaUtilities](https://github.com/CliMA/ClimaUtilities.jl/tree/main):

- `ClimaCoupler.Regridder.remap_field_cgll_to_rll`
- `ClimaCoupler.Regridder.land_fraction`
- `ClimaCoupler.Regridder.combine_surfaces_from_sol!`
- `ClimaCoupler.Regridder.write_datafile_cc`
- `ClimaCoupler.Regridder.read_remapped_field`
- `ClimaCoupler.Regridder.get_coords`

All the above functions can be found in commit
`9e5bf061f34659188485f066bc322c77bcc0f1fa`

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
