# Coupler Object

The CouplerMachine defines a type ```CplState``` for a _container_ variable that holds information about the field 
values that are being used to couple between components. Components can use a ```put!``` operation to 
export a set of field values to a ```CplState``` variable. A ```get``` operation is used to retrieve
a set field values from a ```CplState``` variable.

## Coupler Object API

```@docs
    CouplerMachine.CplState
    CouplerMachine.register_cpl_field!
    CouplerMachine.put!
    CouplerMachine.get
```
