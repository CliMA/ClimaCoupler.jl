# Diagnostics

This module contains functions for defining, gathering and outputting model diagnostics from the Coupler. 

Note that "diagnostics" refers to direct model output generated online (i.e., as the model runs), usually on the model grid.

Offline post-processing (i.e., processing model output after the model is run), such as regridded data to lat-lon, will be done elsewhere. 

## Diagnostics API

```@docs
    ClimaCoupler.Diagnostics.AbstractOutputGroup
    ClimaCoupler.Diagnostics.DiagnosticsGroup
    ClimaCoupler.Diagnostics.AbstractDiagnosticsOperations
    ClimaCoupler.Diagnostics.TimeMean
    ClimaCoupler.Diagnostics.get_var
    ClimaCoupler.Diagnostics.accumulate_diagnostics!
    ClimaCoupler.Diagnostics.save_diagnostics
```


## Diagnostics Internal Functions

```@docs
    ClimaCoupler.Diagnostics.collect_diags
    ClimaCoupler.Diagnostics.iterate_operations
    ClimaCoupler.Diagnostics.operation
    ClimaCoupler.Diagnostics.pre_save
    ClimaCoupler.Diagnostics.post_save
    ClimaCoupler.Diagnostics.save_time_format
```