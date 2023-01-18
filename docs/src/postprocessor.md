# PostProcessor

This module contains functions for postprocessing model data that was saved during the simulation 
by `ClimaCoupler.Diagnostics`. This module is used for offline regridding, slicing and spatial 
averages. It can also handle data from other sources (e.g., NCEP reanalysis). 

## Diagnostics API

```@docs
ClimaCoupler.PostProcessor.postprocess
ClimaCoupler.PostProcessor.PostProcessedData
ClimaCoupler.PostProcessor.ZLatLonData
ClimaCoupler.PostProcessor.ZLatData
ClimaCoupler.PostProcessor.LatLonData
ClimaCoupler.PostProcessor.LatData
ClimaCoupler.PostProcessor.RawData
ClimaCoupler.PostProcessor.DataPackage

```

## Diagnostics Internal Functions

```@docs
ClimaCoupler.PostProcessor.read_remapped_field
```
