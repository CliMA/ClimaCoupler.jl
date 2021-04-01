# CouplerMachine.jl

CouplerMachine.jl provides means to couple climate model components from and within
[ClimateMachine.jl](https://github.com/CliMA/ClimateMachine.jl) and [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl). Functionality includes:
- coupled system time stepping control
- support for mapping import and export boundary information between components.

```@docs
    CouplerMachine
```

```@contents
Pages = [
    "timestepping.md",
    "couplerstate.md",
]
```