# **Coupler stepping**

# Run Sequences (in model time)
![](https://earthsystemmodeling.org/assets/images/nuopc/coupling_leapfrog.png)
![](https://earthsystemmodeling.org/assets/images/nuopc/coupling_explicit.png)

# Modes of execution (in execution time)
There are two modes in which we can run the coupler:
1. **sequential** 
    - models are called in series (e.g. ECMWF IFS)
    - we currently use this with the leap-frog sequencing for prototyping
2. **concurrent** 
    - models are called in parallel (e.g. OASIS-MCT)
    - to be used with the simple explicit sequencing
    - parallelisation tools:
        - multi-node: via MPI (e.g., each model on a different node)
        - single-node: via julia multithreading and shared arrays
    - layout - workload optimisation optimisation
    - lagging - to avoid race conditions

# Plan 
- single-node parallel toy-model prototype *(end Feb)*
    - dependency on `ClimaSimulation.jl` inteface, so each model interrupted at some model time via `DiffEq` callbacks
    - preliminary steps
        1. initialize ICs, BCs and domains as in the **standalone** versions
        2. define for each direction that remapping is required (e.g. model1>model2, model2>model1)
        3. decide on best parallel layout (can this be automated?)
        4. `Simulation(::CoupledSimulation)`
            - define coupler fields and `parameters`
            - override BCs and ICs of stadalone setups
            - `ODEProblem` per model initializing the `rhs!`
            - `model_integ = init(problem, solver, dt, ...)` per model is initialized
            - concurrent loop
                - pre-model ops (e.g., reset flux accumulations, `remap` operators)
                - solve!(model1_integ, callback...)
                - solve!(model2_integ, callback...)
                - post-model ops (barrier until all models done, update coupler fields)
                - ...
    - include performance scaling tests
- single-node SCM / sea breeze case *(March)*
- multi-node parallel prototype (??)
    - dependency on MPI implementation

- future research - [time parallelisation](https://computing.llnl.gov/projects/sundials)

# Ref:
- https://earthsystemmodeling.org/nuopc/
- OASIS: 
    - [craig et al 17](https://gmd.copernicus.org/articles/10/3297/2017/)
    - [OASIS User guide](https://gitlab.com/cerfacs/oasis3-mct/-/blob/OASIS3-MCT_5.0/doc/oasis3mct_UserGuide.pdf)