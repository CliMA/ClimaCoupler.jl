# Diagnostics

This module contains functions for defining, gathering and outputting model diagnostics from the Coupler.

Note that `ClimaCoupler.Diagnostics` is deployed online (i.e., as the model runs), working with cached model data (usually) on the model grid. This does not include offline post-processing (i.e., manipulating saved model output after the model is run, such as regridding data to the latitude-longitude grid). See `ClimaCoupler.PostProcessor` for offline model data treatment.

## Diagnostics API

```@docs
    ClimaCoupler.Diagnostics.AbstractOutputGroup
    ClimaCoupler.Diagnostics.DiagnosticsGroup
    ClimaCoupler.Diagnostics.AbstractDiagnosticsOperations
    ClimaCoupler.Diagnostics.TimeMean
    ClimaCoupler.Diagnostics.get_var
    ClimaCoupler.Diagnostics.accumulate_diagnostics!
    ClimaCoupler.Diagnostics.save_diagnostics
    ClimaCoupler.Diagnostics.init_diagnostics
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