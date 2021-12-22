# ClimaCoupler Design Doc

## Introduction

The `CliMACoupler` is designed primarily to support communication between different model components of the CliMA ESM, using a unified interface of ClimaSimulations.jl. However, the design is such that any Julia-based models should be feasible to couple assuming that they are physically compatible. As well as providing the most optimal coupling functionality for the current CliMA needs, our goal is to achieve optimal software interface modularity/generality, simplicity, performance and robustness, backed with with test cases that demonstrate its use for prototyping, debugging and running model simulations.

#### Code: https://github.com/CliMA/ClimaCoupler.jl

## Design sections
- [terminology](terminology.md)
- [interface](interface.md)
- [coupler timestepping](timestepping.md)
- [regridding & masking](regridding.md)
- [coupled boundary fluxes](boundary_fluxes.md)
- [quality control](quality_control.md)
- [test cases](test_cases.md)




