# TimeManager

This module contains functions that handle dates and times
in simulations. The functions in this module often call
functions from Julia's [Dates](https://docs.julialang.org/en/v1/stdlib/Dates/) module.

## TimeManager API

```@docs
ClimaCoupler.TimeManager.current_date
ClimaCoupler.TimeManager.trigger_callback
ClimaCoupler.TimeManager.AbstractFrequency
ClimaCoupler.TimeManager.Monthly
ClimaCoupler.TimeManager.EveryTimestep
ClimaCoupler.TimeManager.trigger_callback!
ClimaCoupler.TimeManager.CouplerCallback
ClimaCoupler.TimeManager.HourlyCallback
ClimaCoupler.TimeManager.MonthlyCallback
ClimaCoupler.TimeManager.update_firstdayofmonth!
ClimaCoupler.TimeManager.dt_cb
ClimaCoupler.TimeManager.do_nothing
```
