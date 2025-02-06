ClimaCoupler.jl Release Notes
============================

`main`
-------

### ClimaCoupler features

#### Postprocessing no longer uses `sim_mode` PR[#1153](https://github.com/CliMA/ClimaCoupler.jl/pull/1153)
Postprocessing now consists of a single function that's used for all simulation
modes. Note that now all available diagnostics will be plotted at the end of
the simulation, where before we specified a subset to plot based on the
simulation type.
This PR also removes the `plot_diagnostics` config option. `use_coupler_diagnostics`
should be used instead.

#### PrescribedOceanSimulation added, includes SST update PR[#1144](https://github.com/CliMA/ClimaCoupler.jl/pull/1144)
This PR is similar to #1141 below, but applies to SST rather than SIC.
SST is now updated within `PrescribedOceanSimulation` methods, rather
than in the driver. This new simulation type reads in prescribed SST data
and stores it directly in its cache without any processing.
Note that this PR results in changed order of operations and slightly different
behavior for AMIP simulations - see the PR description for more details.

#### Prescribed SIC updated in PrescribedIceSimulation PR[#1141](https://github.com/CliMA/ClimaCoupler.jl/pull/1141)
Previously, SIC was read in directly from the driver during init and in the
coupling loop. Now, it is read internally within the `PrescribedIceSimulation`
object initialization and `step!`. Note that this results in changed order of
operations and slightly different behavior for AMIP simulations - see the PR
description for more details.

#### CO2 is now updated by ClimaAtmos directly PR[#1143](https://github.com/CliMA/ClimaCoupler.jl/pull/1143)

Previously, CO2 was read in directly from the driver during init and in the
coupling loop. Now, it is read internally by ClimaAtmos. Note that this results
in changed order of operations and slightly different behavior for AMIP
simulations - see the PR description for more details.


#### Simplified callbacks PR [#1121](https://github.com/CliMA/ClimaCoupler.jl/pull/1121)

Callbacks were also reworked, and the previous system was removed. Here is an example of the
new way to create a callback that runs every month:
```julia
import ClimaCouper.TimeManager: Callback, maybe_trigger_callback
import ClimaDiagnostics.Schedules: EveryCalendarDtSchedule # This will be moved to ClimaUtilities
import Dates

start_date = Dates.DateTime(1912, 4, 15)
func(cs) = println(maxima(cs.fields.T_s))
schedule = EveryCalendarDtSchedule(Dates.Month(1); start_date)
callback = Callback(schedule, func)

# Then, call it with
maybe_trigger_callback(callback, cs, t)
```

#### Constructor renaming PR[#1135](https://github.com/CliMA/ClimaCoupler.jl/pull/1135)

Simulation constructor functions have been renamed to use the simulation name
itself, following general convention for constructor naming. For example,
`atmos_init` is now `ClimaAtmosSimulation`, and `bucket_init` is now `BucketSimulation`.

### Maintenance

#### Experiment-related tests moved into `experiments/`

Tests for code in the `experiments/` folder, including for component models and
debug plotting have been moved into the `experiments/ClimaEarth/` and
`experiments/ClimaCore` folders. This allows us to remove packages from the `test/`
environment and keep the `experiments/` environments more self-contained.
Github Actions have been added for the new `experiments/ClimaCore/` tests.

v0.1.2
-------
### ClimaEarth features

#### Read bucket initial conditions from NetCDF files

Added functionality to allow the bucket initial conditions to be overwritten by interpolated NetCDF datasets.
To use this feature from the YAML interface, just pass the path of the file to `land_initial_condition`.
We expect the file to contain the following variables:
`W`, for subsurface water storage (2D),
`Ws`, for surface water content (2D),
`T`, for soil temperature (3D),
`S`, for snow water equivalent (2D).

#### Sea-surface temperature and sea ice concentration data can now be automatically downloaded

Sea-surface temperature and sea ice concentration require external files. Now, a
low-resolution version of such files is automatically downloaded when a
higher-resolution version is not available. Please, refer to
[ClimaArtifacts](https://github.com/CliMA/ClimaArtifacts) for more information.

#### A higher resolution land-sea mask is now used and automatically downloaded - PR [#1006](https://github.com/CliMA/ClimaCoupler.jl/pull/1006)

A 60 arcsecond land-sea mask constructed from topographic data is now used.
Topographic data is automatically downloaded and a land-sea mask is constructed
by identifying where elevation is greater than 0. Note, this can lead to
misidentification of ocean in some areas of the globe that are inland but below
sea level (Dead Sea, Death Valley, ...).

#### Leaderboard for variables over longitude, latitude, time, and pressure - PR [#1094](https://github.com/CliMA/ClimaCoupler.jl/pull/1094)

As a part of the post processing pipeline, bias plots for variables at the
pressure levels of 850.0, 500.0, 250.0 hPa and bias plots over latitude and
pressure levels are being created.

### Breaking changes

#### `hourly_checkpoint_dt` is now just `checkpoint_dt` - PR[#1121](https://github.com/CliMA/ClimaCoupler.jl/pull/1121)

Previously, the checkpointing frequency had to be specified in hours and it was
not possible to checkpoint with smaller frequencies. Now, the argument
`hourly_checkpoint_dt` was renamed to `checkpoint_dt` and any frequency can be
specified (e.g., "2months", "200secs").

### Code cleanup

#### Abstract types representing simulation modes & postprocessing changes - PR [#1120](https://github.com/CliMA/ClimaCoupler.jl/pull/1120)

The available simulation modes are now represented by the following abstract types:

- `ClimaCoupler.Interfacer.AMIPMode`
- `ClimaCoupler.Interfacer.SlabplanetMode`
- `ClimaCoupler.Interfacer.SlabplanetAquaMode`
- `ClimaCoupler.Interfacer.SlabplanetTerraMode`
- `ClimaCoupler.Interfacer.SlabplanetEisenmanMode`

All of the above types are subtypes of the abstract
`ClimaCoupler.Interfacer.AbstractSlabplanetSimulationMode`, and all of them except
`ClimaCoupler.Interfacer.AMIPMode` are subtypes of `ClimaCoupler.Interfacer.AbstractSlabplanetSimulationMode`.

These types are used in `experiments/ClimaEarth/run_amip.jl` instead of representing the
simulation mode as a string.

The postprocessing in `experiments/ClimaEarth/run_amip.jl` is now moved into functions in
the new `experiments/ClimaEarth/user_io/postprocessing.jl`. When the simulation is complete,
`postprocess_sim` is called using the type representing the simulation mode for dispatch.
`postprocess_sim` has one method for all slabplanet simulation modes and another for the AMIP
simulation mode. All postprocessing common to all simulation modes is done in the `common_postprocessing`
function, which is called by both `postprocess_sim` methods.

#### Output path updates - PRs [#1058][1], [#1106][2], [#1123][3]

[1]: https://github.com/CliMA/ClimaCoupler.jl/pull/1058
[2]: https://github.com/CliMA/ClimaCoupler.jl/pull/1106
[3]: https://github.com/CliMA/ClimaCoupler.jl/pull/1123
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
```
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
