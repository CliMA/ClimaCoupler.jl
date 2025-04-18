ClimaCoupler.jl Release Notes
============================

`main`
-------

### ClimaCoupler features

#### Removed hierarchy experiments. PR[#1277](https://github.com/CliMA/ClimaCoupler.jl/pull/1277)

The hierarchy experiments have been removed. The last commit that contains them
is
[a6557a3](https://github.com/CliMA/ClimaCoupler.jl/commit/a6557a3bd5853e099429c6f3dda4644c3e28c0d0).

#### Switch to `PartitionedStateFluxes` by default. PR[#1117](https://github.com/CliMA/ClimaCoupler.jl/pull/1117)

Fixed `PartitionedStateFluxes` option. Now `PartitionedStateFluxes` is the
default: instead of combining the surface states and computing fluxes once, we
compute surface fluxes for each component and combine them. Results might be
different. The `CombinedStateFluxes` option will be removed very soon.

#### Split `setup_and_run` in multiple functions. PR[#1251](https://github.com/CliMA/ClimaCoupler.jl/pull/1251)

`setup_and_run` was split into three functions:
- `CoupledSimulation`, which takes a dictionary of a file path and constructs a coupled simulation
- `run!`, which evolves the model
- `postprocess`, which makes plots

In addition to this, `step!` was introduced to take a single step forward with
the coupled simulation.

The function `setup_and_run` is still available.

This change also renames `calendar_dt` to `diagnostics_dt` to make it clearer
that it refers to diagnostics.

#### Some misc. cleanup PR[#1244](https://github.com/CliMA/ClimaCoupler.jl/pull/1244)
Changes include
- Land simulation constructors no longer take in `domain_type`, which was unused.
- `SurfaceModelSimulation`s no longer have a `domain` field, which were unused.
- `CoupledSimulation` now stores the current time as `t`, and does not store the
current date. Its `dates` field is replaced with `date0`.
- `TimeManager` `strdate_to_datetime` and `datetime_to_strdate` are removed, which were unused.

#### Add support for relative parameter filepaths PR[#1228](https://github.com/CliMA/ClimaCoupler.jl/pull/1228)
Changed TOML parameter file handling to prepend the `pkgdir(ClimaCoupler)`
if no file is found at the relative filepath. Before this change, all files
were assumed to be within the `ClimaCoupler` or `ClimaAtmos` repositories.

#### Add support for parameter files in `BucketSimulation` PR[#1217](https://github.com/CliMA/ClimaCoupler.jl/pull/1217)
Add a keyword argument `parameter_files` to `BucketSimulation` to enable
calibration in a coupled simulation, passed via the `"coupler_toml"` argument.

#### Add `ClimaLandSimulation` object PR[#1199](https://github.com/CliMA/ClimaCoupler.jl/pull/1199)
Add methods to support running `ClimaLand.LandModel` in a coupled simulation.
Also add tests to verify the constructor setup and taking a step.
This type is not yet tested within a coupled simulation, but much
of the necessary software infrastructure is added in this PR.

#### Add default `get_field` methods for surface models PR[#1210](https://github.com/CliMA/ClimaCoupler.jl/pull/1210)
Add default methods for `get_field` methods that are commonly
not extended for surface models. These return reasonable default
values, and can be extended by surface models that won't use the
defaults (e.g. the full land model).

#### Add coupler fields based on simulation type PR[#1207](https://github.com/CliMA/ClimaCoupler.jl/pull/1207)
Previously, the coupler fields were hardcoded to be the same for all
simulations, independent of what components were included. Now, each
component model specifies the coupler fields it requires for coupling,
and these are used to construct the set of coupler fields.
TOA radiation and net precipitation are added only if conservation is enabled.
The coupler fields are also now stored as a ClimaCore Field of NamedTuples,
rather than as a NamedTuple of ClimaCore Fields.

#### Restart simulations with JLD2 files PR[#1179](https://github.com/CliMA/ClimaCoupler.jl/pull/1179)

`ClimaCoupler` can now use `JLD2` files to save state and cache for its model
component, allowing it to restart from saved checkpoints. Some restrictions
apply:

- The number of MPI processes has to remain the same across checkpoints
- Restart files are generally not portable across machines, julia versions, and package versions
- Adding/changing new component models will probably require adding/changing code

Please, refer to the
[documentation](https://clima.github.io/ClimaCoupler.jl/dev/checkpointer/) for
more information.

#### Remove extra `get_field` functions PR[#1203](https://github.com/CliMA/ClimaCoupler.jl/pull/1203)
Removes the `get_field` functions for `air_density` for all models, which
were unused except for the `BucketSimulation` method, which is replaced by a
function call to `extrapolate_ρ_to_sfc`. Also removes the `get_field` function
for `ClimaAtmosSimulation` `air_temperature`, which was unused.

#### Coupler fields surface variable names all use `_sfc` PR[#1195](https://github.com/CliMA/ClimaCoupler.jl/pull/1195)
`T_S`, `z0m_S`, and `z0b_S` are renamed to `T_sfc`, `z0m_sfc`, and `z0b_sfc`.
This makes them consistent with the other surface fields, e.g. `ρ_sfc` and `q_sfc`.

#### Driver function `setup_and_run` added PR[#1178](https://github.com/CliMA/ClimaCoupler.jl/pull/1178)
A new function `setup_and_run` is added, which takes in a path to a config file,
and contains all the code to initialize component models and run the simulation.
This function is needed for calibration to easily run multiple simulations
with different parameters. This function is moved to a different file
`experiments/ClimaEarth/setup_run.jl`, so it can be included by the `run_amip.jl`
driver, or by another driver used for calibration.

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
maybe_trigger_callback(callback, cs)
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
