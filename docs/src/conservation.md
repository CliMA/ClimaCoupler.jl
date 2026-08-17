# ConservationChecker

If the model is a physically closed system (e.g., in the `slabplanet` configuration with free slip conditions),
it should conserve mass (including water), energy and momentum. The conservation checker logs global conservation.

Only energy and water are currently implemented, as the conserved quantities
`TotalEnergy` and `TotalWater`. Set `conservation_check: true` in the
configuration file to enable them.

Note that kinetic energy is not included in the calculation of the global energy, reflecting the formulation on `ClimaAtmos`,
which assumes that kinetic energy is negligible in comparison with the moist static energy components.

Each `ConservationCheck` holds one timeseries per component simulation, plus one
for each part of the budget that no component reports: energy lost through the
top of the atmosphere, water shed as runoff, and water a component holds without
ever having received it as a flux. The total is not stored separately. Every
timeseries carries the sign that makes summing them give the total, so anything
removed from the budget is stored negated.

## ConservationChecker API

```@docs
ConservationChecker.ConservationCheck
ConservationChecker.check_conservation!
```
