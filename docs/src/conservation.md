# ConservationChecker

If the model is a physically closed system (e.g., in the `slabplanet` configuration with free slip conditions),
it should conserve mass (including water), energy and momentum. The conservation checker logs global conservation.

Only energy and water are currently implemented, as the conserved quantities
`TotalEnergy` and `TotalWater`. Set `conservation_check: true` in the
configuration file to enable them. Each `ConservationCheck` holds one
timeseries per component simulation, plus others that are tracked to ensure that
the budget closes (e.g., energy lost through the top of the atmosphere). 

Note that kinetic energy is not included in the calculation of the global energy, reflecting the formulation on `ClimaAtmos`,
which assumes that kinetic energy is negligible in comparison with the moist static energy components.

## ConservationChecker API

```@docs
ConservationChecker.ConservationCheck
ConservationChecker.check_conservation!
```
