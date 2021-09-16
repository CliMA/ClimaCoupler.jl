# ClimaCore Experiments

- prerequisites: Julia 1.6+
- these experiments represent the basic cases for testing physical soundness of coupling methods in [ClimaCore.jl](https://github.com/CliMA/ClimaCore.jl/)
    - 1D finite-difference column test cases 
        - Integration tests for individual functionalities
            - TC1: Dry heat diffusion + slab 
                - Purpose: 
                    - minimal example for prototype developments and integration tests
                    - flux calculation / accumulation inside the atmos model `rhs!`
                    - one-file tutorial
                - to be maintained by CI
            - TC2: Ekman column + slab
                - uses ClimaAtmos interface
                - not maintained by CI (discontinued upon final clean up of TC3)     
            - TC4: Dry Ekman Column + slab + surface fluxes
                - uses local ClimaAtmos hack
                - Purpose:
                    - for prototyping ClimaCore-like equations with some lightweight ClimaAtmos interface wrappers
                    - all flux calculations / accumulations inside atmos model `rhs!`
                    - compare effect of different `calculate_sfc_fluxes_energy` functions
                - maintained by CI 

        - Interface tests
            - TC3: Ekman column + Land (+ Oceananigans)
                - Purpose:
                    - combine all three interfaces
                    - showcase ClimaCoupler interface functions
                    - let flux calculation be done by coupler
                    - tested with SurfaceFluxes.jl

- documentation and test results can be found in our [Design Doc](https://www.overleaf.com/read/bgfmhgtncpws). 
