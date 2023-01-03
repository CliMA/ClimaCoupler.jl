# BCReader

This module coordinates reading of boundary conditions from NetCDF files,
as well as regridding calls and temporal interpolations from
monthly to daily intervals.

## BCReader API

```@docs
ClimaCoupler.BCReader.BCFileInfo
ClimaCoupler.BCReader.bcfile_info_init
ClimaCoupler.BCReader.update_midmonth_data!
ClimaCoupler.BCReader.next_date_in_file
ClimaCoupler.BCReader.interpolate_midmonth_to_daily
```


## BCReader Internal Functions

```@docs
ClimaCoupler.BCReader.no_scaling
ClimaCoupler.BCReader.interpol
```