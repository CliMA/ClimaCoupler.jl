ClimaCoupler.jl Release Notes
============================

`main`
-------

### ClimaCoupler features

#### Add a SimCoordinator module PR[#1676](https://github.com/CliMA/ClimaCoupler.jl/pull/1676)
Move the functions `run!(cs)` and `step!(cs)` into a new module, SimCoordinator.
These are exported at the top-level ClimaCoupler.jl, so they can continue to be used
as before. This PR also moves `postprocess(cs)` into the Plotting module.

#### Update interfaces to surface flux calculator PR[#1646](https://github.com/CliMA/ClimaCoupler.jl/pull/1646)
Uses new interface in SurfaceFluxes v0.15 to calculate surface fluxes. For the bucket model, call the
turbulent_fluxes! function in ClimaLand directly for surface flux calculation.

#### Change default simulation output directory PR[#1653](https://github.com/CliMA/ClimaCoupler.jl/pull/1653)
Changes the default output directory from `experiments/ClimaEarth/output/` to `output`.
The internal structure within the output directory remains the same.

#### Update interfaces to surface-flux calculator to support dynamic roughness PR[#1567](https://github.com/CliMA/ClimaCoupler.jl/pull/1567)
Adds the `roughness_model` argument to the surface-flux input constructor. Default
behaviour is the `SF.ScalarRoughness()` type. Upcoming changes will account for the
value of the `roughness_model` to be inferred from the mask value (e.g. prescribed
aerodynamic roughness over land, and computed as a function of shear over oceans).

#### Remove the `bucket_ic_august.nc` initial condition file PR[#1564](https://github.com/CliMA/ClimaCoupler.jl/pull/1564)
Remove the file `experiments/ClimaEarth/input/bucket_ic_august.nc`, which
was originally added for DYAMOND simulations and is no longer needed.
We still keep the config file option `bucket_initial_condition`.

#### Add an option `restart_cache` PR[#1548](https://github.com/CliMA/ClimaCoupler.jl/pull/1548)
Adds an option so that if restart files are available, we can choose to
restart the state only (`restart_cache` false), or to restart both the state
and the cache (`restart_cache` true). `restart_cache` is true by default.

#### Remove `FluxCalculator.surface_inputs` helper function PR[#1543](https://github.com/CliMA/ClimaCoupler.jl/pull/1543)
We can simplify the flux calculation by calling `SF.ValuesOnly` directly.
Since we now remap all quantities onto the boundary space when we compute
fluxes, there's no need to access the underlying data layouts of ClimaCore Fields.

#### Provide `SW_d`, `LW_d` to surface models instead of `F_radiative` PR[#1518](https://github.com/CliMA/ClimaCoupler.jl/pull/1518)
This allows us to correctly compute radiative flux over surface models by
computing each contribution individually (`SW_u, SW_d, LW_u, LW_d`).

#### Change some closures used in the Oceananigans model PR[#1524](https://github.com/CliMA/ClimaCoupler.jl/pull/1524)
Per ocean team recommendation, we changed some closures to be non-default.

#### Change the behavior of `detect_restart_file` PR[#1515](https://github.com/CliMA/ClimaCoupler.jl/pull/1515)
Users now only need to specify `restart_dir` and `restart_t` to restart a simulation from a specific file, and do
not need to set `detect_restart_file` to true. `detect_restart_file` is used for detecting restart file automatically.

#### Fixes for sea ice PR[#1519](https://github.com/CliMA/ClimaCoupler.jl/pull/1519)
Don't weight fluxes by area fraction when passing them to the surface models,
only when providing them to the atmosphere.
Zero out tendencies of prescribed sea ice and slab ocean where the area
fraction is zero.
For now, use binary area fractions for all surface models, until we correctly
handle area fraction weighting.

#### Use `update_turbulent_fluxes!` instead of `update_field!` for atmosphere PR[#1511](https://github.com/CliMA/ClimaCoupler.jl/pull/1511)
Instead of using an `update_field!` method that dispatches on `::Val{:turbulent_fluxes}`
to update turbulent fluxes in the atmosphere, we switch to using a function `update_turbulent_fluxes!`.
This is consistent with what we do for surface models.

This PR also removes area fraction weighting in the default radiative flux update
in FluxCalculator.jl.

#### Replace `TD.PhaseEquil_ρTq` with `TD.PhaseNonEquil_ρTq` PR[#1506](https://github.com/CliMA/ClimaCoupler.jl/pull/1506)
This should be more physically correct.

#### Use individual surface model temperatures to compute turbulent fluxes PR[#1498](https://github.com/CliMA/ClimaCoupler.jl/pull/1498)
We use a partitioned approach to computing turbulent fluxes, so we should be using
the individual surface temperature (and humidity, air density, etc) for each component
model. Previously we were using a combined surface temperature across all surface models.
This change doesn't affect simulations using integer area fractions, but will affect simulations
with fractional area fractions.

This PR also adds ensures that the area fraction of an OceananigansSimulation
never exceeds the latitude limits when using a LatitudeLongitudeGrid. Note that
this simulation requires running with an ice model, which is used to fill in
the polar regions.

#### Combine LW fluxes to get surface temperature for radiation PR[#1492](https://github.com/CliMA/ClimaCoupler.jl/issues/1492)
Previously, we combined the temperatures of each surface models directly using an
area-weighted sum. Now, we instead compute the longwave flux for each component, compute
the area fraction-weighted sum of those, and then convert the total flux back to get
surface temperature. This temperature is then provided to the atmosphere for radiation.

#### Remove coarse nightly CMIP PR[#1485](https://github.com/CliMA/ClimaCoupler.jl/pull/1485)
To avoid depending on the main branch of too many packages in the nightly pipeline,
we remove the CMIP nightly run and will only test AMIP nightly.

#### Use JuliaFormatter; remove .dev folder PR[#1484](https://github.com/CliMA/ClimaCoupler.jl/pull/1484)
Use JuliaFormatter to format the repo, rather than the previous `.dev/juliaformat.jl`.
The `.dev/` folder is removed, and a documentation page about contributing to
ClimaCoupler.jl is added.

#### Correctly set land domain with `share_surface_space = true` PR[#1464](https://github.com/CliMA/ClimaCoupler.jl/pull/1464)
Previously, we passed but did not use arguments related to the vertical resolution
in `make_land_domain`. Now we correctly use the provided `dz_tuple` and `n_elements_vert`.

#### Add option `detect_restart_files` PR[#1463](https://github.com/CliMA/ClimaCoupler.jl/pull/1463)
Add a CLI option to signal whether restart files should be automatically used
to restart a simulation. This is false by default.

#### Remove ED/EDMF aquaplanet longruns PR[#1461](https://github.com/CliMA/ClimaCoupler.jl/pull/1461)
Removes the ED-only and diag. EDMF aquaplanet longruns.
These can be run manually as needed for debugging, rather than running every week.
Also increases the RMSE limit for `rsutcs` from 7.4 to 10.8.

#### Remove `atmos_config_repo` PR[#1448](https://github.com/CliMA/ClimaCoupler.jl/pull/1448)
ClimaCoupler doesn't depend on configuration files in ClimaAtmos anymore. All the atmosphere
configuration files are now specified in ClimaCoupler.

#### Remove bucket `get_new_cache` PR[#1437](https://github.com/CliMA/ClimaCoupler.jl/pull/1437)

As of ClimaLand v0.16.2, total energy and and water are always stored in the bucket cache.
This PR removes the `get_new_cache`, which allocated space for the energy field,
and instead accesses the cached fields directly.

#### Add option for integrated land spun up IC. PR[#1318](https://github.com/CliMA/ClimaCoupler.jl/pull/1318)

Adds a boolean flag `land_spun_up_ic` that allows the user to request reading integrated land
initial conditions from a saved simulation, rather than setting them with default values.
If this option is true, the land/sea mask must be used, since the spun-up ICs are only defined
over land (i.e. this cannot be used in the `slabplanet_terra` mode). The default is true.

#### Use EN4 dataset for ocean initial conditions and forcing PR[#1425](https://github.com/CliMA/ClimaCoupler.jl/pull/1425)

Use EN4 instead of ECCO for ocean initial conditions and forcing, to avoid
authentication requirements.

#### Don't plot constant fields in debug plots PR[#1417](https://github.com/CliMA/ClimaCoupler.jl/pull/1417)

Heatmaps are no longer generated for spatially constant fields in our debug plots.
This is done to avoid errors when generating the plots.

#### Add aquaplanets; longrun fixes PR[#1411](https://github.com/CliMA/ClimaCoupler.jl/pull/1411)

Add 2 aquaplanet longruns using slab ocean and atmosphere with ED only and diag.
EDMF respectively. Also removes the atmos standalone longrun, and increases
the conservation RSE threshold.

#### Remove intermediate checkpoints PR[#1397](https://github.com/CliMA/ClimaCoupler.jl/pull/1397)

Throughout the simulation, the previous checkpoint is now deleted whenever a new
one is saved. The most recent checkpoint will always be available, so restarting
is still supported. The field `prev_checkpoint_t` in the CoupledSimulation object
is used to remove intermediate checkpoints.

#### Remove `dt_save_state_to_disk` and `dt_save_to_sol` options PR[#1394](https://github.com/CliMA/ClimaCoupler.jl/pull/1394)

`dt_save_state_to_disk` was unused and is removed from all configs in this
commit. Note that ClimaAtmos does have an option with this name, but we
pass `checkpoint_dt` to it. `dt_save_to_sol` is also removed as an option,
in favor of using our more robust checkpointing infrastructure via `checkpoint_dt`.

#### Misc. interface cleanup PR[#1341](https://github.com/CliMA/ClimaCoupler.jl/pull/1341)

Including:
- Remove `ρ_sfc` from surface model caches
- Remove `F_radiative` from `add_coupler_fields` since it's a default
- Remove atmosphere `get_field` method for `thermo_state_int`, since it's now constructed on the boundary space
- Update Interfacer docs
- `comms_ctx` and `boundary_space` are no longer in the `CoupledSimulation` object. Use accessor functions instead.

#### Construct thermo states from exchanged T, q, ρ. PR[#1293](https://github.com/CliMA/ClimaCoupler.jl/pull/1293)

Instead of exchanging atmosphere and surface thermo states, exchange T, q, ρ
and use these to construct the thermo states. This will be necessary as we
run with models on different grids and require real remapping.

#### Allow < 3 surface models. PR[#1286](https://github.com/CliMA/ClimaCoupler.jl/pull/1286)

Previously, land, ocean, and sea ice were all required to be defined, and the
area fraction of a model would be set to 0 to exclude it. Now, a `CoupledSimulation`
can be created with fewer than 4 components in `model_sims`.

#### Change signature for `turbulent_fluxes!`. PR[#1327](https://github.com/CliMA/ClimaCoupler.jl/pull/1327)

`FluxCalculator.turbulent_fluxes!` can now be called in two ways:
```julia
turbulent_fluxes!(coupled_simulation::CoupledSimulation)
turbulent_fluxes!(coupler_fields, model_sims, thermo_params)
```
The previous signature was
```julia
turbulent_fluxes!(model_sims, coupler_fields, boundary_space, thermo_params)
```
The new signature simplifies calling the function, ensures that the mutating
convention is respected (`turbulent_fluxes!` mutates `coupler_fields`), and
removes unnecessary arguments and type restrictions.

Similarly, `get_surface_fluxes!` was renamed to `get_surface_fluxes`.

#### Simplify initial component model exchange. PR[#1305](https://github.com/CliMA/ClimaCoupler.jl/pull/1305)

Surface humidity is now computed from the atmosphere state and surface temperature,
rather than sometimes computing it and sometimes retrieving it from a model cache.
This allows us to simplify the initial component exchange, since we don't need
to `step!` component models to get them to compute humidity. Without `step!`,
we can also remove `reinit!`.

#### Add `get_field` with `remap`. PR[#1314](https://github.com/CliMA/ClimaCoupler.jl/pull/1314)

PR [#1314](https://github.com/CliMA/ClimaCoupler.jl/pull/1314) introduces a new
method for `get_field`: `get_field(sim, Val(), target_space)`. This new method
automatically remaps onto `target_space`. Currently, only trivial remapping is
supported (ie, the same behavior as `dummmy_remap`).

Similarly, `get_field!(target_field, sim, Val())` gets the field in-place.

#### Turbulent energy flux is split into LHF, SHF. PR[#1309](https://github.com/CliMA/ClimaCoupler.jl/pull/1309)
Previously, we have exchanged the combined turbulent energy flux `F_turb_energy`;
This PR splits up the exchanged quantity into `F_lh` and `F_sh`.
This is helpful because some component models store the two fluxes separately,
and those that use the combined quantity can easily be updated with the sum.

#### Rename bucket-specific options. PR[#1310](https://github.com/CliMA/ClimaCoupler.jl/pull/1310)

`land_albedo_type` is now `bucket_albedo_type`, and `land_initial_condition` is now `bucket_initial_condition`.

#### Rename bucket-specific options. PR[#1306](https://github.com/CliMA/ClimaCoupler.jl/pull/1306)

`Interfacer.name` was removed. Use `nameof` instead.

#### Test AMIP with integrated land model. PR[#1254](https://github.com/CliMA/ClimaCoupler.jl/pull/1254)

The integrated ClimaLand model can now be used in coupled simulations.
A short run of AMIP using the full land model is now tested in our regular
Buildkite pipeline, and the restarts test uses the full land model
rather than the bucket.

This PR adds the config option `land_model`, which can be set to either
`bucket` or `integrated` to choose which land model to run with.

#### Remove `nans_to_zero`. PR[#1278](https://github.com/CliMA/ClimaCoupler.jl/pull/1278)

Instead of zeroing out all NaNs in a surface field, we zero out all values
where the area fraction is 0, which is logically what we want to do.

#### Removed `SurfaceScheme`. PR[#1280](https://github.com/CliMA/ClimaCoupler.jl/pull/1280)

The `BulkScheme` option for computing fluxes was removed. Now, fluxes
are always computed with the `MoninObukhovScheme`.

#### Removed `CombinedStateFluxes`. PR[#1276](https://github.com/CliMA/ClimaCoupler.jl/pull/1276)

The `CombinedStateFluxes` option for computing fluxes was removed. Now, fluxes
are always computed with the option formerly known as `PartitionedStateFluxes`
(no longer an option). Related code was also removed, such as the
`atmos_turbulent_fluxes_most!` or `combined_turbulent_fluxes!` functions,
and the `TurbulentFluxPartition` types.

The `partitioned_turbulent_fluxes!` function was renamed to `turbulent_fluxes!`.

As a result of this, the signature of certain functions has changed:
`update_sim`, `update_model_sims!`, `import_atmos_fields!`, and
`import_combined_surface_fields!` no longer take the `turbulent_fluxes` argument.


#### Remove `area_mask`, `binary_mask`, `mono_surface`. [PR#1268](https://github.com/CliMA/ClimaCoupler.jl/pull/1268/files)
Removes all uses of `area_mask`, as multiplying quantities by both `area_fraction`
and `area_mask` will potentially lead to physically inaccurate results.
It also defeats the purpose of maintaining that the sum of surface models'
area fractions is 1 at all points. See issue [#1255](https://github.com/CliMA/ClimaCoupler.jl/issues/1255) for more details.
The function `Utilities.binary_mask` is also removed, as it's now unused.
The option `mono_surface` is no longer supported, as it was rarely exercised
and did not do what it was documented to do. All simulations now have
behavior equivalent to using `mono_surface: false` previously.

#### Remove `EisenmanSeaIce`. PR[#1284](https://github.com/CliMA/ClimaCoupler.jl/pull/1284)

The `EisenmanSeaIce` was removed. The last commit that contains this model is
[a3b32d1](https://github.com/CliMA/ClimaCoupler.jl/commit/a3b32d169137f7dad2edf33fd2f5e29ebd6d5356).
Please refer to this commit if you are interested in running this model. The
function `FluxCalculator.differentiate_turbulent_fluxes!` was no longer needed
and was removed.

#### Removed hierarchy experiments. PR[#1277](https://github.com/CliMA/ClimaCoupler.jl/pull/1277)

The hierarchy experiments have been removed. The last commit that contains them
is
[a6557a3](https://github.com/CliMA/ClimaCoupler.jl/commit/a6557a3bd5853e099429c6f3dda4644c3e28c0d0).

#### Switch to `PartitionedStateFluxes` by default. PR[#1117](https://github.com/CliMA/ClimaCoupler.jl/pull/1117)

Fixed `PartitionedStateFluxes` option. Now `PartitionedStateFluxes` is the
default: instead of combining the surface states and computing fluxes once, we
compute surface fluxes for each component and combine them. Results might be
different.

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
