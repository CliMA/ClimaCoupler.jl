# Conservation Checks

If the model is a physically closed system (e.g., in the `slabplanet` configuration with free slip conditions), it should conserve mass (including water), energy and momentum. The conservation checker logs global conservation.

Only energy and water are currently implemented.

Note that kinetic energy is not included in the calculation of the global energy, reflecting the formulation on `ClimaAtmos`, which assumes that kinetic energy is negligible in comparison with the moist static energy components.

## ConservationChecker API

```@docs
ClimaCoupler.ConservationChecker.EnergyConservationCheck
ClimaCoupler.ConservationChecker.WaterConservationCheck
ClimaCoupler.ConservationChecker.check_conservation!
```

## ConservationChecker Internal Functions
```@docs
ClimaCoupler.ConservationChecker.surface_water_gain_from_rates
```
