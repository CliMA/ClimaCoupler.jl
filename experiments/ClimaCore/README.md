# ClimaCore Experiments

- prerequisites: Julia 1.6+
- these experiments represent the basic cases for testing physical soundness of coupling methods in [ClimaCore.jl](https://github.com/CliMA/ClimaCore.jl/)
    - 1D finite-difference column test cases
        - Integration tests for individual functionalities
            - Dry heat diffusion + slab (`experiments/ClimaCore/heat-diffusion/`)
                - Purpose:
                    - minimal example for prototype developments and integration tests
                    - flux calculation / accumulation inside the atmos model `rhs!`
                    - one-file tutorial

        - Interface tests
            - Ekman column + Land (+ Oceananigans) (`experiments/ClimaCore/atm-ocn-lnd/`)
                - Purpose:
                    - combine all three interfaces
                    - showcase ClimaCoupler interface functions
                    - let flux calculation be done by coupler
                    - tested with SurfaceFluxes.jl

- documentation and test results can be found in our [Design Doc](https://www.overleaf.com/read/bgfmhgtncpws).
