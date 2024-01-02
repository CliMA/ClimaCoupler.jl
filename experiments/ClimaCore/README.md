# ClimaCore Experiments

- prerequisites: Julia 1.6+
- these experiments represent the basic cases for testing physical soundness of coupling methods in [ClimaCore.jl](https://github.com/CliMA/ClimaCore.jl/)
    - 1D finite-difference column test cases
        - Integration tests for individual functionalities
            - Dry heat diffusive atmos + slab surface (`experiments/ClimaCore/heat-diffusion/`)
                - Purpose:
                    - minimal example for prototype developments and integration tests
                    - flux calculation / accumulation inside the atmos model `rhs!`
                    - one-file tutorial
    - 2D box test cases
        - Integration test - seabreeze experiment
            - Dry heat diffusive atmos + slab ocean + slab land (`experiments/ClimaCore/sea_breeze/`)
                - Purpose:
                    - demonstrate coupling of a simple atmosphere model with multiple surface models
                    - expand domain from 1D to 2D case
                    - flux calculations, exchange, and accumulation

- documentation and test results can be found in our [Design Doc](https://www.overleaf.com/read/bgfmhgtncpws).
