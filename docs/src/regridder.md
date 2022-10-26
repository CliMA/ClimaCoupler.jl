# Regridder

This module contains functions to regrid information between spaces.
Many of the functions used in this module call TempestRemap functions
via ClimaCoreTempestRemap wrappers.

Information about the TempestRemap library can be found [here](https://github.com/ClimateGlobalChange/tempestremap).
Multiple remapping approaches from TempestRemap have been tested with our
implementation, and information about them is located [here](https://github.com/CliMA/ClimaCoupler.jl/wiki/ClimaCoupler-Lessons-Learned).

## Regridder API

```@docs
ClimaCoupler.Regridder.write_to_hdf5
ClimaCoupler.Regridder.read_from_hdf5
ClimaCoupler.Regridder.dummmy_remap!
ClimaCoupler.Regridder.remap_field_cgll_to_rll
ClimaCoupler.Regridder.land_sea_mask
ClimaCoupler.Regridder.update_masks!
ClimaCoupler.Regridder.combine_surfaces!
```


## Regridder Internal Functions

```@docs
ClimaCoupler.Regridder.reshape_cgll_sparse_to_field!
ClimaCoupler.Regridder.hdwrite_regridfile_rll_to_cgll
ClimaCoupler.Regridder.write_datafile_cc
ClimaCoupler.Regridder.binary_mask
```