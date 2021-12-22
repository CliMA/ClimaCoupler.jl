# Coupler stepping

## Run Sequences (in model time)
![](https://earthsystemmodeling.org/assets/images/nuopc/coupling_leapfrog.png)
![](https://earthsystemmodeling.org/assets/images/nuopc/coupling_explicit.png)

## Modes of execution (in execution time)
There are two modes in which we can run the coupler:
- sequential - models are called in series (e.g. ECMWF IFS)
    - we use this with the leap-frog sequencing for prototyping
- concurrent - models are called in parallel (e.g. OASIS-MCT)
    - to be used with the simple explicit sequencing
    - parallelisation tools:
        - multi-node: via MPI (e.g., each model on a different node)
        - single-node: via julia multithreading and shared arrays
    - layout - workload optimisation optimisation
    - lagging - to avoid race conditions


# Ref:
- https://earthsystemmodeling.org/nuopc/